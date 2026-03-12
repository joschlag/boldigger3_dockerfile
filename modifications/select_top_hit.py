import duckdb, datetime, more_itertools, re, time
import pandas as pd
import dask.dataframe as dd
import numpy as np
from tqdm import tqdm
from boldigger3.id_engine import parse_fasta
from string import punctuation, digits


def clean_dataframe(dataframe: object) -> object:
    # replace missing values and empty strings in metadata to pd.NA
    metadata_columns = [
        "processid",
        "sex",
        "life_stage",
        "inst",
        "country/ocean",
        "identified_by",
        "identification_method",
        "coord",
        "nuc",
        "marker_code",
    ]

    # clean na values
    dataframe[metadata_columns] = dataframe[metadata_columns].replace(
        [None, "None", ""], pd.NA
    )
    # remove all punctuation and digits except '-', some species names contain a "-"
    specials = re.escape("".join(c for c in punctuation + digits if c != "-"))
    pattern = f"[{specials}]"

    # Levels to clean
    levels = ["phylum", "class", "order", "family", "genus", "species"]

    # clean Na and None values from levels
    dataframe[levels] = dataframe[levels].replace([None, "None", ""], pd.NA)

    # replace any occurrence of invalid characters with pd.NA
    for level in levels:
        col = dataframe[level].astype("string")
        mask = col.str.contains(pattern, regex=True, na=pd.NA)
        dataframe[level] = col.where(~mask)

    # replace all empty strings in dataframe with pd.NA
    dataframe = dataframe.replace("", pd.NA)

    # make all object columns strings so we do not get unintended behaviour later
    object_columns = dataframe.select_dtypes(include="object").columns
    dataframe[object_columns] = dataframe[object_columns].astype("string")

    try:
        # extract the lat lon values
        dataframe[["lat", "lon"]] = (
            dataframe["coord"].str.strip("[]").str.split(",", expand=True)
        )
        dataframe["lat"], dataframe["lon"] = dataframe["lat"].astype(
            "float"
        ), dataframe["lon"].astype("float")
    except ValueError:
        dataframe["lat"], dataframe["lon"] = np.nan, np.nan

    # drop coord column
    dataframe = dataframe.drop("coord", axis=1)

    return dataframe


def stream_hits_to_excel(id_engine_db_path, project_directory, fasta_dict, fasta_name):
    # chunk the fasta dicts keys to retrieve from duckdb
    chunks = enumerate(more_itertools.chunked(fasta_dict.keys(), n=8_000), start=1)

    # define the output path
    output_path = project_directory.joinpath("boldigger3_data")

    with duckdb.connect(id_engine_db_path) as connection:
        # retrieve one chunk of a maximum of 10_000 process_ids
        for part, chunk in chunks:
            query = f"""SELECT * FROM final_results 
            WHERE id IN ?
            ORDER BY fasta_order ASC, pct_identity DESC"""
            chunk_data = connection.execute(query, [chunk]).df()
            chunk_data = clean_dataframe(chunk_data)

            # drop the fasta order just before saving
            chunk_data = chunk_data.drop("fasta_order", axis=1)

            chunk_data.to_excel(
                output_path.joinpath(f"{fasta_name}_bold_results_part_{part}.xlsx"),
                index=False,
                engine="xlsxwriter",
            )


def get_threshold(hit_for_id: object, thresholds: list) -> object:
    """Function to find a threshold for a given id from the complete dataset.

    Args:
        hit_for_id (object): The hits for the respective id as dataframe.
        thresholds (list): Lists of thresholds to to use for the selection of the top hit.

    Returns:
        object: Single line dataframe containing the top hit
    """
    # find the highest similarity value for the threshold
    threshold = hit_for_id["pct_identity"].max()

    # check for no matches first
    if "no-match" in hit_for_id.astype(str).values:
        return 0, "no-match"
    else:
        # move through the taxonomy if it is no nomatch hit or broken record
        if threshold >= thresholds[0]:
            return thresholds[0], "species"
        elif threshold >= thresholds[1]:
            return thresholds[1], "genus"
        elif threshold >= thresholds[2]:
            return thresholds[2], "family"
        elif threshold >= thresholds[3]:
            return thresholds[3], "order"
        elif threshold >= thresholds[4]:
            return thresholds[4], "class"
        else:
            return thresholds[5], "phylum"


def move_threshold_up(threshold: int, thresholds: list) -> tuple:
    """Function to move the threshold up one taxonomic level safely."""
    levels = ["species", "genus", "family", "order", "class", "phylum"]

    idx = thresholds.index(threshold)
    if idx + 1 >= len(thresholds):
        # Already at highest threshold, return the same
        return threshold, levels[idx]

    return thresholds[idx + 1], levels[idx + 1]


def flag_hits(top_hits: object, final_top_hit: object):
    # initialize the flags
    flags = [""] * 5

    # flag 1: Reverse BIN taxonomy
    id_method = top_hits["identification_method"].dropna()

    if (
        not id_method.empty
        and id_method.str.contains("BOLD|ID|Tree|BIN", regex=False).all()
    ):
        flags[0] = "1"

    # flag 2: top hit ratio < 90%
    #if final_top_hit["records_ratio"].item() < 0.9:
    #if final_top_hit["records_ratio"].iloc[0] < 0.9:
    if not final_top_hit.empty and final_top_hit["records_ratio"].iloc[0] < 0.9:
        flags[1] = "2"

    # flag 3: all of the selected top hits are private
    if top_hits["status"].isin(["private"]).all():
        flags[2] = "3"

    if len(top_hits.index) == 1:
        flags[3] = "4"

    # flag 5: top hit is represented by multiple bins
    #if len(final_top_hit["BIN"].str.split("|").item()) > 1:
    bin_value = final_top_hit["BIN"].iloc[0]
    if pd.notna(bin_value) and len(str(bin_value).split("|")) > 1:
        flags[4] = "5"

    flags = "|".join(flags)

    return flags


def find_top_hit(hits_for_id: object, thresholds: list) -> object:
    """Funtion to find the top hit for a given ID.

    Args:
        hits_for_id (object): Dataframe with the data for a given ID
        thresholds (list): List of thresholds to perform the top hit selection with.

    Returns:
        object: Single line dataframe with the selected top hit
    """

    iteration = 0
    # safety limit to avoid infinite loops
    max_iterations = 20  

    # get the thrshold and taxonomic level
    threshold, level = get_threshold(hits_for_id, thresholds)

    # if a nomatch is found, a no-match can directly be retured
    if level == "no-match":
        return_value = hits_for_id.query("species == 'no-match'").head(1)
        fasta_order = return_value["fasta_order"]
        # columns to return
        return_value = return_value[
            [
                "id",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
                "pct_identity",
                "status",
            ]
        ]

        # fill the missing data with correct types
        data_to_type = {
            "records": 0,
            "selected_level": pd.NA,
            "BIN": pd.NA,
            "flags": "||||",
        }
        for key, value in data_to_type.items():
            return_value[key] = value

        # add the fasta order back in
        return_value["fasta_order"] = fasta_order

        return_value = return_value.astype(
            {
                "selected_level": "string[python]",
                "BIN": "string[python]",
                "flags": "string[python]",
            }
        )

        return return_value

    # define the levels for the groupby. care about the selector string later
    all_levels = ["phylum", "class", "order", "family", "genus", "species"]
    # go through the hits to make the selection
    while iteration < max_iterations:
        iteration += 1
        # copy the hits to perform modifications
        hits_above_similarity = hits_for_id[hits_for_id["pct_identity"] > threshold].copy()
        # -----------------------------------------------------------------
        # Robust handling for mixed taxonomy hits
        # If species-level records exist, ignore genus-only rows.
        # This prevents mixed rows like:
        #   Asaphes vulgaris
        #   Asaphes
        # -----------------------------------------------------------------
        if level == "species":
            species_rows = hits_above_similarity[
                hits_above_similarity["species"].notna()
                & (hits_above_similarity["species"].astype(str).str.strip() != "")
            ]
            if not species_rows.empty:
                hits_above_similarity = species_rows
        # -----------------------------------------------------------------
        levels = all_levels[: all_levels.index(level) + 1]
        # select the hits above similarity and of interest
        hits_above_similarity = hits_above_similarity[levels].copy()
        hits_above_similarity = (
            hits_above_similarity.groupby(by=levels, sort=False, dropna=False)
            .size()
            .reset_index(name="count")
        )
        hits_above_similarity = hits_above_similarity.dropna(subset=[level], axis=0)


        if hits_above_similarity.empty:
            try:
                threshold, level = move_threshold_up(threshold, thresholds)
                continue
            except IndexError:
                fallback = hits_for_id.head(1).copy()
                fallback["records"] = 1
                fallback["records_ratio"] = 1.0
                fallback["selected_level"] = level
                fallback["BIN"] = ""
                fallback["flags"] = ""
                return fallback[
                    [
                        "id",
                        "phylum",
                        "class",
                        "order",
                        "family",
                        "genus",
                        "species",
                        "pct_identity",
                        "status",
                        "records",
                        "records_ratio",
                        "selected_level",
                        "BIN",
                        "flags",
                        "fasta_order",
                    ]
                ]

        top_hit_row = hits_above_similarity.iloc[0]
        top_count = top_hit_row["count"]
        top_ratio = top_count / hits_above_similarity["count"].sum()

        #query_parts = [f"`{col}` == '{str(top_hit_row[col]).replace("'", "''")}'" for col in top_hit_row.index[:-1]]
        query_parts = []
        for col in top_hit_row.index[:-1]:
            value = str(top_hit_row[col]).replace("'", "''")
            query_parts.append(f"`{col}` == '{value}'")
        
        query_string = " and ".join(query_parts)
        top_hits = hits_for_id.query(query_string)

        final_top_hit = top_hits.head(1).copy()
        final_top_hit["records"] = top_count
        final_top_hit["records_ratio"] = top_ratio
        final_top_hit["selected_level"] = level

        # safe BIN collection
        top_hit_bins = top_hits["bin_uri"].dropna().unique() if threshold == thresholds[0] else []
        final_top_hit["BIN"] = "|".join(top_hit_bins) if len(top_hit_bins) > 0 else pd.NA

        if threshold != thresholds[0]:
            idx = all_levels.index(level)
            levels_to_remove = all_levels[idx + 1 :]
            final_top_hit[levels_to_remove] = pd.NA
            final_top_hit[levels_to_remove] = final_top_hit[levels_to_remove].astype("string")

        final_top_hit["flags"] = flag_hits(top_hits, final_top_hit)

        final_top_hit = final_top_hit[
            [
                "id",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
                "pct_identity",
                "status",
                "records",
                "records_ratio",
                "selected_level",
                "BIN",
                "flags",
                "fasta_order",
            ]
        ]
        return final_top_hit

    # fallback if max iterations reached
    fallback = hits_for_id.head(1).copy()
    fallback["records"] = 1
    fallback["records_ratio"] = 1.0
    fallback["selected_level"] = level
    fallback["BIN"] = ""
    fallback["flags"] = ""
    return fallback[
        [
            "id",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
            "pct_identity",
            "status",
            "records",
            "records_ratio",
            "selected_level",
            "BIN",
            "flags",
            "fasta_order",
        ]
    ]


def gather_top_hits(
    fasta_dict, id_engine_db_path, project_directory, fasta_name, thresholds
):
    # store top hits here until n are reached, flush to parquet inbetween
    top_hits_buffer = []
    buffer_counter = 0

    with duckdb.connect(id_engine_db_path) as connection:
        for query_id in tqdm(fasta_dict.keys(), desc="Top hit calculation"):
            sql_query = f"SELECT * FROM final_results WHERE id='{query_id}' ORDER BY fasta_order ASC, pct_identity DESC"
            hits_for_id = clean_dataframe(connection.execute(sql_query).df())

            if hits_for_id.empty:
                print(f"[WARN] No hits found for ID: {query_id}")
                continue

            try:
                result = find_top_hit(hits_for_id, thresholds)
                if result is not None and not result.empty:
                    top_hits_buffer.append(result)

            except Exception as e:
                print(f"[ERROR] Failed processing ID {query_id}: {e}")

                # ---------------- DEBUG BLOCK ----------------
                print("\n----- DEBUG HITS START -----")
                print(f"Query ID: {query_id}")
                print("Number of hits:", len(hits_for_id))

                debug_cols = [
                    "id",
                    "phylum",
                    "class",
                    "order",
                    "family",
                    "genus",
                    "species",
                    "pct_identity",
                    "status",
                    "bin_uri",
                ]

                existing_cols = [c for c in debug_cols if c in hits_for_id.columns]
                print(hits_for_id[existing_cols].head(20))
                print("----- DEBUG HITS END -----\n")
                # ------------------------------------------------

                continue

            # flush buffer
            if len(top_hits_buffer) >= 1_000:
                parquet_output = project_directory.joinpath(
                    "boldigger3_data",
                    f"{fasta_name}_top_hit_buffer_{buffer_counter}.parquet.snappy"
                )
                pd.concat(top_hits_buffer, axis=0).reset_index(drop=True).to_parquet(parquet_output)
                buffer_counter += 1
                top_hits_buffer = []

        # final buffer flush
        if top_hits_buffer:
            parquet_output = project_directory.joinpath(
                "boldigger3_data",
                f"{fasta_name}_top_hit_buffer_{buffer_counter}.parquet.snappy",
            )
            pd.concat(top_hits_buffer, axis=0).reset_index(drop=True).to_parquet(parquet_output)
            #top_hits_buffer = pd.concat(top_hits_buffer, axis=0).reset_index(drop=True)
            #top_hits_buffer.to_parquet(parquet_output)


def save_results(project_directory, fasta_name):
    # load all data into dask dataframe
    data_paths = project_directory.joinpath("boldigger3_data").glob(
        f"{fasta_name}_top_hit_buffer_*.parquet.snappy"
    )
    data_paths = list(data_paths)

    if not data_paths:
        print("[WARN] No top hits generated.")
        return

    all_top_hits = dd.read_parquet([str(f) for f in data_paths])

    # order values globally by fasta order, drop afterwards
    all_top_hits = (
        all_top_hits.set_index("fasta_order", sorted=True)
        .reset_index(drop=True)
        .compute()
    )

    # write the output
    parquet_output = project_directory.joinpath(
        "boldigger3_data", "{}_identification_result.parquet.snappy".format(fasta_name)
    )
    excel_output = project_directory.joinpath(
        "boldigger3_data", "{}_identification_result.xlsx".format(fasta_name)
    )

    # save the data
    all_top_hits.to_excel(excel_output, index=False, engine="xlsxwriter")
    all_top_hits.to_parquet(parquet_output)

    # unlink the buffer files
    for file in project_directory.joinpath("boldigger3_data").glob(
        f"{fasta_name}_top_hit_buffer_*.parquet.snappy"
    ):
        if file.exists():
            file.unlink()


def main(fasta_path: str, thresholds: list):
    tqdm.write(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Removing digits and punctuation from hits."
    )

    # load the fasta data
    fasta_dict, fasta_name, project_directory = parse_fasta(fasta_path)

    # define the id engine database path
    id_engine_db_path = project_directory.joinpath(
        "boldigger3_data", f"{fasta_name}.duckdb"
    )

    tqdm.write(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Streaming all hits to excel."
    )

    # # stream the data from duckdb to excel first
    stream_hits_to_excel(id_engine_db_path, project_directory, fasta_dict, fasta_name)

    tqdm.write(f"{datetime.datetime.now().strftime('%H:%M:%S')}: Calculating top hits.")

    gather_top_hits(
        fasta_dict, id_engine_db_path, project_directory, fasta_name, thresholds
    )
    tqdm.write(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Saving results. This may take a while."
    )

    save_results(project_directory, fasta_name)

    tqdm.write(f"{datetime.datetime.now().strftime('%H:%M:%S')}: Finished.")
