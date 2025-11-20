import datetime, duckdb, sys, pickle, more_itertools, requests_html, json, time, os
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from tqdm import tqdm
from collections import OrderedDict
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from boldigger3.exceptions import DownloadFinished
from json.decoder import JSONDecodeError
from requests.exceptions import ReadTimeout


class BoldIdRequest:
    """A class to represent the data for a BOLD id engine request

    Attributes
    ----------

    Methods
    -------
    """

    def __init__(
        self,
        # base_url: str,
        # params: dict,
        # query_generator: object,
        # timestamp: object,
        # result_url: str,
        # last_checked: object
    ):
        """Constructs the neccessary attribues for the BoldIdRequest object

        Args:
            base_url (str): Represents the base url for the post request
            params (dict): The parameters to send with the post request
            query_generator (object): A generator holding the data to send with the post request
            timestamp (object): Timestamp that is set when the request is sent to BOLD
            result_url (str): The result url to download the data from
            last_checked (object): The last time the download url was checked for updates

        """
        self.base_url = ""
        self.params = {}
        self.query_data = []
        self.result_url = ""
        self.timestamp = None
        self.database = None
        self.operating_mode = None
        self.download_url = ""
        self.last_checked = None


def parse_fasta(fasta_path: str) -> tuple:
    """Function to read a fasta file and parse it into a dictionary.

    Args:
        fasta_path (str): Path to the fasta file to be identified.

    Returns:
        tuple: Data of the fasta file in a dict, the full path to the fasta file, the directory where this fasta file is located.
    """
    # extract the directory from the fasta path
    fasta_path = Path(fasta_path)
    fasta_name = fasta_path.stem
    project_directory = fasta_path.parent

    # use SeqIO to read the data into dict- automatically check fir the type of fasta
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    # trim header to maximum allowed chars of 99. names are preserved in the SeqRecord object
    fasta_dict = {key: value for key, value in fasta_dict.items()}

    # create a set of all valid DNA characters
    valid_chars = {
        "A",
        "C",
        "G",
        "T",
        "M",
        "R",
        "W",
        "S",
        "Y",
        "K",
        "V",
        "H",
        "D",
        "B",
        "X",
        "N",
    }

    # check all sequences for invalid characters
    raise_invalid_fasta = False

    # check if the sequences contain invalid chars
    for key in fasta_dict.keys():
        if not set(fasta_dict[key].seq.upper()).issubset(valid_chars):
            print(
                f"{datetime.datetime.now().strftime('%H:%M:%S')}: Sequence {key} contains invalid characters."
            )
            raise_invalid_fasta = True

    if not raise_invalid_fasta:
        return fasta_dict, fasta_name, project_directory
    else:
        sys.exit()


def already_downloaded(fasta_dict: dict, database_path: str) -> dict:
    """Function to check if any of the requests has been downloaded and stored in the duckdb database.

    Args:
        fasta_dict (dict): The dictionary with the fasta data.
        database_path (str): Path to the duckdb database

    Returns:
        dict: The dictionary with the fasta data with already downloaded sequences removed.
    """
    if database_path.is_file():
        # connect to the database
        database = duckdb.connect(database_path)

        # collect all unique values from the id column
        downloaded_ids = database.execute(
            "SELECT DISTINCT id FROM id_engine_results"
        ).fetchall()

        # close the connection
        database.close()

        # flatten the downloaded ids into a set for faster lookup
        downloaded_ids = {row[0] for row in downloaded_ids}

        # filter the fasta_dict
        fasta_dict = {
            id: seq for (id, seq) in fasta_dict.items() if id not in downloaded_ids
        }

        # return the updated dict
        return fasta_dict
    else:
        # return the fasta dict unchanged
        return fasta_dict


def build_url_params(database: int, operating_mode: int) -> tuple:
    """Function that generates a base URL and the params for the POST request to the ID engine.

    Args:
        database (int): Between 1 and 7 referring to the database, see readme for details.
        operating_mode (int): Between 1 and 3 referring to the operating mode, see readme for details

    Returns:
        tuple: Contains the base URL as str and the params as dict
    """

    # the database int is translated here
    idx_to_database = {
        1: "public.tax-derep",
        2: "species",
        3: "all.tax-derep",
        4: "DS-CANREF22",
        5: "public.plants",
        6: "public.fungi",
        7: "all.animal-alt",
        8: "DS-IUCNPUB",
    }

    # the operating mode is translated here
    idx_to_operating_mode = {
        1: {"mi": 0.94, "maxh": 25},
        2: {"mi": 0.9, "maxh": 50},
        3: {"mi": 0.75, "maxh": 100},
    }

    # params can be calculated from the database and operating mode
    params = {
        "db": idx_to_database[database],
        "mi": idx_to_operating_mode[operating_mode]["mi"],
        "mo": 100,
        "maxh": idx_to_operating_mode[operating_mode]["maxh"],
        "order": 3,
    }

    # format the base url
    base_url = f"https://id.boldsystems.org/submission?db={params['db']}&mi={params['mi']}&mo={params['mo']}&maxh={params['maxh']}&order={params['order']}"

    return base_url, params


def build_download_queue(fasta_dict: dict, database: int, operating_mode: int) -> dict:
    """Function to build the download queue.

    Args:
        fasta_dict (dict): Dict that holds the data in the fasta file.
        download_queue_name (str): String that holds the path where the download queue is saved.
        database (int): Between 1 and 7 referring to the database, see readme for details.
        operating_mode (int): Between 1 and 3 referring to the operating mode, see readme for details

    Returns:
        dict: The dictionary with the downloaded queue

    """
    # initialize the download queue
    download_queue = {"waiting": OrderedDict(), "active": dict()}

    # build the base url and the params
    base_url, params = build_url_params(database, operating_mode)

    # determine the query size from the params
    query_size_dict = {0.94: 1000, 0.9: 200, 0.75: 100}
    query_size = query_size_dict[params["mi"]]

    # split the fasta dict in query sized chunks
    query_data = more_itertools.chunked(fasta_dict.keys(), query_size)

    # produce a generator that holds all sequence and key data to loop over for the post requests
    query_generators = (
        [f">{key}\n{fasta_dict[key].seq}\n" for key in query_subset]
        for query_subset in query_data
    )

    for idx, query_generator in enumerate(query_generators, start=1):
        # initialize the bold id engine request
        bold_request = BoldIdRequest()
        bold_request.base_url = base_url
        bold_request.params = params
        bold_request.query_data = query_generator
        bold_request.database = database
        bold_request.operating_mode = operating_mode
        download_queue["waiting"][idx] = bold_request

    return download_queue


def build_post_request(BoldIdRequest: object) -> object:
    """Function to send the POST request for the dataset to the BOLD id engine.

    Args:
        bold_id_request (object): A BoldIdRequest object that holds all the information needed to send the request

    Returns:
        object: Returns the BoldIdRequest object with an added result url
    """
    # send the post requests
    with requests_html.HTMLSession() as session:

        # build a retry strategy for the html session
        session.headers.update(
            {
                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/98.0.4758.82 Safari/537.36"
            }
        )
        retry_strategy = Retry(total=10, backoff_factor=1)
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("https://", adapter)

        data = "".join(BoldIdRequest.query_data)

        # generate the files to send via the id engine
        files = {"fasta_file": ("submitted.fas", data, "text/plain")}

        while True:
            try:
                # submit the post request
                response = session.post(
                    BoldIdRequest.base_url,
                    params=BoldIdRequest.params,
                    files=files,
                    timeout=10000,
                )

                # fetch the result
                result = json.loads(response.text)
                break
            except (JSONDecodeError, ReadTimeout):
                # user output
                tqdm.write(
                    f"{datetime.datetime.now().strftime('%H:%M:%S')}: Building the request failed. Waiting 60 seconds for repeat."
                )
                # wait 60 seconds
                time.sleep(60)

        result_url = f"https://id.boldsystems.org/submission/results/{result['sub_id']}"

        # append the resulting url
        BoldIdRequest.result_url = result_url
        BoldIdRequest.timestamp = datetime.datetime.now()

        return BoldIdRequest


def add_no_match(result: object, BoldIdRequest: object, fasta_order: dict) -> dict:
    """Function to add a no match in case BOLD does not return any results

    Args:
        result (object): JSON returned from BOLD.
        BoldIdRequest (object): Id Request object.
        fasta_order (dict): Order of the original fasta file

    Returns:
        dict: A no match row for the dataframe
    """
    seq_id = result["seqid"]

    row = {
        "id": seq_id,
        "phylum": "no-match",
        "class": "no-match",
        "order": "no-match",
        "family": "no-match",
        "subfamily": "no-match",
        "genus": "no-match",
        "species": "no-match",
        "taxid_count": 0.0,
        "pct_identity": 0.0,
        "process_id": "",
        "bin_uri": "",
        "request_date": pd.Timestamp.now().strftime("%Y-%m-%d %X"),
        "database": BoldIdRequest.database,
        "operating_mode": BoldIdRequest.operating_mode,
        "status": "",
        "fasta_order": fasta_order[seq_id],
    }

    return row


def safe_status(record_key, index, default=""):
    return record_key[index] if len(record_key) > index else default


def parse_and_save_data(
    BoldIdRequest: object,
    response: object,
    fasta_order: dict,
    request_id: int,
    database_path: str,
    fasta_name: str,
):
    """Function to parse the JSON returned by BOLD and save it as parquet.

    Args:
        BoldIdRequest (object): BoldIdRequest object.
        response (object): http response to parse.
        fasta_order (dict): Order of the original fasta file, can be used to order the table after metadata addition.
        request_id (int): Request id, used to save the file.
        database_path (str): Path to the database to write to.
    """
    # parse out all objects from the response
    json_objects = [json.loads(line) for line in response.iter_lines()]

    # store the resulting tables rows here when parsing the jsons
    rows = []

    for result in json_objects:
        # handle no-matches here
        if not result["results"]:
            # create a no-match hit
            row = add_no_match(result, BoldIdRequest, fasta_order)
            rows.append(row)
            continue

        # parse the json
        seq_id = result["seqid"]

        for record_key, record_data in result["results"].items():
            record_key = record_key.split("|")

            row = {
                "id": seq_id,
                "process_id": record_key[0],
                "bin_uri": record_key[2],
                "status": safe_status(record_key, 4),
                "pct_identity": 100.0 - record_data.get("pdist"),
            }  # definition changed, now 100 - pdist

            taxonomy = record_data.get("taxonomy", {})
            row.update(taxonomy)
            rows.append(row)

    # transform to dataframe, drop unused labels, add ordering column, use correct type annotations
    id_engine_result = pd.DataFrame(rows)
    id_engine_result = id_engine_result.drop(
        labels=["subfamily", "taxid_count"], axis=1
    )
    id_engine_result["pct_identity"] = id_engine_result["pct_identity"].astype(
        "float64"
    )
    id_engine_result["request_date"] = pd.Timestamp.now().strftime("%Y-%m-%d %X")
    id_engine_result["database"] = BoldIdRequest.database
    id_engine_result["operating_mode"] = BoldIdRequest.operating_mode
    id_engine_result["fasta_order"] = id_engine_result["id"].map(fasta_order)

    # reorder the columns
    id_engine_result = id_engine_result[
        [
            "id",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
            "pct_identity",
            "process_id",
            "bin_uri",
            "request_date",
            "database",
            "operating_mode",
            "status",
            "fasta_order",
        ]
    ]

    # finally stream to parquet to later load into duckdb, update the active queue
    output_file = database_path.joinpath(
        "boldigger3_data", f"request_id_{request_id}_{fasta_name}.parquet.snappy"
    )

    id_engine_result.to_parquet(output_file)


# def download_json(
#     active_queue: dict, fasta_order: dict, project_directory: str, fasta_name: str,
# ):
#     """Function to download the JSON results from the id engine and store them in temporary parquet files.

#     Args:
#         active_queue (dict): Queue with active BOLD requests.
#         fasta_order (dict): Dict that can be appended to the results. Needed for creating a sorted duckdb table later.
#     """
#     # initialize the last check's for the active queue
#     for key in active_queue.keys():
#         if not active_queue[key].last_checked:
#             active_queue[key].last_checked = datetime.datetime.now()

#     # check the requests every ten seconds
#     with requests_html.HTMLSession() as session:
#         # loop over the active queue until any request is finished
#         while active_queue:
#             for key in active_queue.keys():
#                 now = datetime.datetime.now()
#                 # check the timestamp of the key, if it is older than 10 minutes, pop it from the active
#                 # queue to fetch it in a later run
#                 if now - active_queue[key].timestamp > datetime.timedelta(minutes=15):
#                     tqdm.write(
#                         f"{datetime.datetime.now().strftime('%H:%M:%S')}: Request ID {key} has timed out. Will be requeued."
#                     )
#                     active_queue.pop(key)
#                     return active_queue

#                 # if the url has not been checked in the last 10 seconds, check again, else skip the url
#                 if now - active_queue[key].last_checked > datetime.timedelta(
#                     seconds=15
#                 ):
#                     response = session.get(active_queue[key].result_url)
#                 else:
#                     continue
#                 # if there's no data in the response yet, continue
#                 if response.status_code == 404:
#                     active_queue[key].last_checked = datetime.datetime.now()
#                     continue
#                 # if timing is sufficient AND data in the response, parse the response
#                 else:
#                     # parse the response here and save, update the active queue
#                     parse_and_save_data(
#                         active_queue[key],
#                         response,
#                         fasta_order,
#                         key,
#                         project_directory,
#                         fasta_name,
#                     )
#                     active_queue.pop(key)

#                     # give user output
#                     tqdm.write(
#                         f"{datetime.datetime.now().strftime('%H:%M:%S')}: Request ID {key} has successfully been downloaded."
#                     )

#                     # return the active queue
#                     return active_queue

def download_json(
        active_queue: dict,
        fasta_order: dict,
        project_directory: str,
        fasta_name: str,
        retry_queue: dict,
):
    """
    Download JSON results for all active jobs.
    Adds robust retry handling, backoff, and permanent-failure detection.
    """

    # Set missing last_checked timestamps
    for key, req in active_queue.items():
        if req.last_checked is None:
            req.last_checked = datetime.datetime.now()

        # ---- NEW: Ensure retry counters always exist ----
        if not hasattr(req, "retry_count"):
            req.retry_count = 0                                    # ← NEW
        if not hasattr(req, "max_retries"):
            req.max_retries = 5                                    # ← NEW

    with requests_html.HTMLSession() as session:
        completed = []

        for key, req in list(active_queue.items()):
            now = datetime.datetime.now()

            # ===============================
            # 1. TIMEOUT HANDLING
            # ===============================
            if now - req.timestamp > datetime.timedelta(minutes=15):

                tqdm.write(f"{now:%H:%M:%S}: Request {key} timed out.")

                req.retry_count += 1

                if req.retry_count <= req.max_retries:
                    # ---- NEW: Exponential backoff ----
                    req.next_attempt = now + datetime.timedelta(
                        seconds=2 ** req.retry_count               # ← NEW
                    )
                    tqdm.write(
                        f"{now:%H:%M:%S}: Retrying {key} "
                        f"({req.retry_count}/{req.max_retries}) – "
                        f"waiting {2 ** req.retry_count}s."        # ← NEW
                    )
                    retry_queue[key] = req
                else:
                    tqdm.write(
                        f"{now:%H:%M:%S}: Request {key} failed permanently."
                    )

                del active_queue[key]
                continue

            # -------- NEW: Backoff check --------
            if hasattr(req, "next_attempt") and now < req.next_attempt:  # ← NEW
                continue                                                  # ← NEW

            # ===============================
            # 2. Respect minimum poll interval
            # ===============================
            if now - req.last_checked < datetime.timedelta(seconds=15):
                continue

            # Try requesting JSON
            try:
                response = session.get(req.result_url)
            except Exception as e:
                tqdm.write(f"{now:%H:%M:%S}: Network error for {key}: {e}")
                # Treat as transient
                req.retry_count += 1
                if req.retry_count <= req.max_retries:
                    req.next_attempt = now + datetime.timedelta(seconds=2 ** req.retry_count)
                    retry_queue[key] = req
                else:
                    tqdm.write(f"{now:%H:%M:%S}: Network failure permanent for {key}.")
                del active_queue[key]
                continue

            req.last_checked = now

            # ===============================
            # 3. Handle BOLD server responses
            # ===============================

            # ---- Not ready yet (normal) ----
            if response.status_code == 404:
                continue

            # ---- NEW: Permanent errors ----
            if response.status_code in {400, 403, 410}:  # ← NEW
                tqdm.write(
                    f"{now:%H:%M:%S}: Permanent failure for {key} "
                    f"(HTTP {response.status_code}), discarding."    # ← NEW
                )
                del active_queue[key]
                continue

            # ---- NEW: Transient server errors ----
            if response.status_code >= 500:  # ← NEW
                req.retry_count += 1         # ← NEW
                if req.retry_count <= req.max_retries:
                    req.next_attempt = now + datetime.timedelta(seconds=2 ** req.retry_count)
                    tqdm.write(
                        f"{now:%H:%M:%S}: Server error {response.status_code} "
                        f"for {key}, retry {req.retry_count}/{req.max_retries} "
                        f"(waiting {2 ** req.retry_count}s)."        # ← NEW
                    )
                    retry_queue[key] = req
                else:
                    tqdm.write(
                        f"{now:%H:%M:%S}: Server error permanent for {key}."
                    )
                del active_queue[key]
                continue

            # ---- SUCCESS ----
            parse_and_save_data(
                req,
                response,
                fasta_order,
                key,
                project_directory,
                fasta_name,
            )
            tqdm.write(f"{now:%H:%M:%S}: Downloaded request {key}.")
            completed.append(key)

        # Remove completed from active queue
        for key in completed:
            active_queue.pop(key, None)

    return active_queue

def parquet_to_duckdb(project_directory, database_path):
    """Function to stream the parquet output to duckdb.

    Args:
        project_directory (str): Project directory to work in.
        database_path (str): Path to the database.
    """
    # check if duck db already exists
    db_exists = database_path.exists()

    # fetch the parquet path
    parquet_path = project_directory.joinpath(
        "boldigger3_data", "request_id_*.parquet.snappy"
    )

    # connect to the database
    database = duckdb.connect(database_path)

    try:
        # insert is db exists already, create table if not
        if not db_exists:
            # Create the table
            database.execute(
                f"""
                CREATE TABLE id_engine_results AS
                SELECT * FROM read_parquet('{parquet_path}')
            """
            )
        else:
            database.execute(
                f"""
            INSERT INTO id_engine_results
            SELECT * FROM read_parquet('{parquet_path}')
            """
            )
    except duckdb.IOException:
        # if no parquet files are there dues to killed downloads, just continue
        pass

    # remove all parquet files if finished successfully
    for file in project_directory.joinpath("boldigger3_data").glob("*.parquet.snappy"):
        if file.is_file():
            file.unlink()

def main(fasta_path: str, database: int, operating_mode: int) -> None:
    """Main function to run the BOLD identification engine."""

    tqdm.write(f"{datetime.datetime.now():%H:%M:%S}: Reading input fasta.")

    fasta_dict, fasta_name, project_directory = parse_fasta(fasta_path)

    fasta_dict_order = {key: idx for idx, key in enumerate(fasta_dict.keys())}

    database_path = project_directory.joinpath("boldigger3_data", f"{fasta_name}.duckdb")

    download_queue_name = project_directory.joinpath(
        "boldigger3_data", f"{fasta_name}_download_queue.pkl"
    )

    data_dir = project_directory.joinpath("boldigger3_data")
    data_dir.mkdir(exist_ok=True)

    # Check for prior downloads
    fasta_dict = already_downloaded(fasta_dict, database_path)

    if not fasta_dict:
        tqdm.write(f"{datetime.datetime.now():%H:%M:%S}: All data has already been downloaded.")
        return None

    # Load existing queue
    try:
        with open(download_queue_name, "rb") as f:
            download_queue = pickle.load(f)
        tqdm.write(f"{datetime.datetime.now():%H:%M:%S}: Found unfinished downloads. Continuing.")

        # ---- FIX 1 — Ensure retry queue exists ----
        if "retry" not in download_queue:
            download_queue["retry"] = OrderedDict()

    except FileNotFoundError:
        tqdm.write(f"{datetime.datetime.now():%H:%M:%S}: Building the download queue.")
        download_queue = build_download_queue(fasta_dict, database, operating_mode)

        # ---- FIX 2 — Ensure retry exists on new queue ----
        download_queue["retry"] = OrderedDict()

        with open(download_queue_name, "wb") as f:
            pickle.dump(download_queue, f)

        tqdm.write(
            f"{datetime.datetime.now():%H:%M:%S}: Added {len(download_queue['waiting'])} "
            f"requests to the download queue."
        )

    total_downloads = len(download_queue["waiting"]) + len(download_queue["active"])

    # ------------------- MAIN DOWNLOAD LOOP -------------------
    with tqdm(total=total_downloads, desc="Finished downloads") as pbar:

        while True:
            try:
                if download_queue["waiting"] or download_queue["active"] or download_queue["retry"]:

                    # ---- FIX 3 — Move retry items safely back to waiting ----
                    if not download_queue["waiting"] and download_queue["retry"]:
                        tqdm.write(
                            f"{datetime.datetime.now():%H:%M:%S}: "
                            f"Moving {len(download_queue['retry'])} retry items back to waiting."
                        )
                        for key, req in list(download_queue["retry"].items()):
                            download_queue["waiting"][key] = req
                            del download_queue["retry"][key]

                    # Fill active queue up to four jobs
                    if len(download_queue["active"]) < 4 and download_queue["waiting"]:
                        request_id, req_obj = download_queue["waiting"].popitem(last=False)

                        tqdm.write(
                            f"{datetime.datetime.now():%H:%M:%S}: "
                            f"Request ID {request_id} moved to active downloads."
                        )

                        download_queue["active"][request_id] = build_post_request(req_obj)

                    else:
                        # ---- FIX 4 — Safe active queue update ----
                        before = len(download_queue["active"])
                        download_queue["active"] = download_json(
                            download_queue["active"],
                            fasta_dict_order,
                            project_directory,
                            fasta_name,
                            download_queue["retry"],
                        )
                        after = len(download_queue["active"])

                        # Only advance the bar for *actual* completed downloads
                        finished_now = before - after
                        if finished_now > 0:
                            pbar.update(finished_now)

                    # Persist queue every loop
                    with open(download_queue_name, "wb") as out_stream:
                        pickle.dump(download_queue, out_stream)

                else:
                    parquet_to_duckdb(project_directory, database_path)
                    raise DownloadFinished

            except DownloadFinished:
                fasta_dict = already_downloaded(fasta_dict, database_path)

                if fasta_dict:
                    tqdm.write(f"{datetime.datetime.now():%H:%M:%S}: Requeuing incomplete downloads.")

                    download_queue = build_download_queue(fasta_dict, database, operating_mode)
                    download_queue["retry"] = OrderedDict()

                    total_downloads = len(download_queue["active"]) + len(download_queue["waiting"])
                    pbar.reset()
                    pbar.total = total_downloads
                    pbar.refresh()

                else:
                    tqdm.write(f"{datetime.datetime.now():%H:%M:%S}: All downloads finished successfully.")
                    os.remove(download_queue_name)
                    break



