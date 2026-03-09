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


def log(level: str, message: str):
    """Simple logging function that writes to tqdm output with timestamp."""
    tqdm.write(f"{datetime.datetime.now():%H:%M:%S} [{level}] {message}")

MIN_BATCH_MINUTES = 30

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
        """Constructs the neccessary attributes for the BoldIdRequest object

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
        
        # ================================================
        # NEW FIELDS REQUIRED BY PATCH 1 + PATCH 2
        # ================================================

        # Count how many times a request has failed
        self.retry_count = 0              # ← Needed for timeouts & 5xx retries

        # Maximum number of retries before permanent failure
        self.max_retries = 3              # ← Set default here

        # Used for exponential backoff (download_json checks this)
        self.next_attempt = None          # ← Needed to avoid hammering the API

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

    # use SeqIO to read the data into dict- automatically check for the type of fasta
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
        raise ValueError("One or more sequences contain invalid characters in the FASTA file.")


def already_downloaded(fasta_dict: dict, database_path: str) -> dict:
    """Function to check if any of the requests has been downloaded and stored in the duckdb database.

    Args:
        fasta_dict (dict): The dictionary with the fasta data.
        database_path (str): Path to the duckdb database

    Returns:
        dict: The dictionary with the fasta data with already downloaded sequences removed.
    """
    # If DB does not exist → nothing downloaded yet
    if not database_path.is_file():
        return fasta_dict

    # Connect to DuckDB
    database = duckdb.connect(database_path)

    try:
        # Try to read the stored IDs
        downloaded_ids = database.execute(
            "SELECT DISTINCT id FROM id_engine_results"
        ).fetchall()

    except duckdb.CatalogException:
        # Table does not exist → treat as empty database
        database.close()
        return fasta_dict

    # Close DB
    database.close()

    # Convert rows to a set of IDs
    downloaded_ids = {row[0] for row in downloaded_ids}

    # Filter out already downloaded IDs
    return {
        seq_id: seq
        for seq_id, seq in fasta_dict.items()
        if seq_id not in downloaded_ids
    }


def build_url_params(database: int, operating_mode: int) -> tuple:
    """
    Returns a clean base URL (without query parameters)
    and a params dict for the POST request.
    """
    # database translation
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

    # operating mode translation
    idx_to_operating_mode = {
        1: {"mi": 0.94, "maxh": 25},
        2: {"mi": 0.9, "maxh": 50},
        3: {"mi": 0.75, "maxh": 100},
    }

    base_url = "https://id.boldsystems.org/submission"

    params = {
        "db": idx_to_database[database],
        "mi": idx_to_operating_mode[operating_mode]["mi"],
        "mo": 100,
        "maxh": idx_to_operating_mode[operating_mode]["maxh"],
        "order": 3,
    }

    return base_url, params


def build_download_queue(fasta_dict: dict, database: int, operating_mode: int) -> dict:
    """Function to build the download queue.

    Args:
        fasta_dict (dict): Dict that holds the data in the fasta file.
        database (int): Between 1 and 7 referring to the database, see readme for details.
        operating_mode (int): Between 1 and 3 referring to the operating mode, see readme for details

    Returns:
        dict: The dictionary with the download queue
    """

    # Initialize queue including retry container
    download_queue = {
        "waiting": OrderedDict(),
        "active": dict(),
        "retry": dict(),
    }

    # Base URL + generic params
    base_url, base_params = build_url_params(database, operating_mode)

    # Determine chunk size
    query_size_dict = {0.94: 250, 0.9: 150, 0.75: 75}
    query_size = query_size_dict[base_params["mi"]]

    # Split fasta sequences into batches
    query_data = more_itertools.chunked(fasta_dict.keys(), query_size)

    # Create the actual payload strings
    query_generators = (
        [f">{key}\n{fasta_dict[key].seq}\n" for key in query_subset]
        for query_subset in query_data
    )

    # Build request objects
    for idx, query_generator in enumerate(query_generators, start=1):
        bold_request = BoldIdRequest()

        bold_request.base_url = base_url

        # ⚠️ IMPORTANT: copy params so each request has its own dict
        bold_request.params = base_params.copy()

        bold_request.query_data = query_generator
        bold_request.database = database
        bold_request.operating_mode = operating_mode

        # No need to reassign retry/backoff fields — already set in __init__()
        # # Retry / backoff fields
        # bold_request.retry_count = 0
        # bold_request.max_retries = 5
        # bold_request.next_attempt = None
        # bold_request.last_checked = None
        # bold_request.timestamp = None

        download_queue["waiting"][idx] = bold_request

    return download_queue

def build_post_request(BoldIdRequest: object) -> object:
    """
    Sends a POST request to the BOLD ID engine.
    Preserves all fields needed for the download loop.

    Changes:
    - Pre-flight DNS/TCP/TLS probe (non-fatal) and logging.             # ← NEW
    - Payload/header debug logging just before POST.                   # ← NEW
    - Attach fasta_header and sequence to request object for logging.  # ← NEW
    """

    # --- ensure retry/backoff attributes exist and have sensible types ---
    for attr in ["retry_count", "max_retries", "next_attempt", "timestamp", "last_checked"]:
        if not hasattr(BoldIdRequest, attr):
            # integer counters default to 0, others default to None
            if attr in ["retry_count", "max_retries"]:
                setattr(BoldIdRequest, attr, 0 if attr == "retry_count" else 5)
            else:
                setattr(BoldIdRequest, attr, None)

    # --- Prepare payload and attach FASTA header/sequence for logging/troubleshooting ---  # ← NEW
    data = "".join(BoldIdRequest.query_data)


    # -------- NEW: avoid duplicate POST if job already submitted --------
    if getattr(BoldIdRequest, "result_url", ""):
        return BoldIdRequest


    # Try to extract the first FASTA header and the first sequence for logging
    fasta_header = "UNKNOWN_HEADER"
    fasta_sequence = ""
    try:
        # split by lines, find first line that starts with '>'
        lines = data.splitlines()
        for i, ln in enumerate(lines):
            if ln.startswith(">"):
                fasta_header = ln[1:].strip()
                # next non-header lines until next '>' or end -> sequence
                seq_lines = []
                for later in lines[i + 1 :]:
                    if later.startswith(">"):
                        break
                    seq_lines.append(later.strip())
                fasta_sequence = "".join(seq_lines)
                break
    except Exception:
        fasta_header = "UNKNOWN_HEADER"
        fasta_sequence = ""

    # attach them so log_failed_request can access them later  # ← NEW
    try:
        BoldIdRequest.fasta_header = fasta_header
        BoldIdRequest.sequence = fasta_sequence
    except Exception:
        # non-fatal - keep going
        pass

    files = {"fasta_file": ("submitted.fas", data, "text/plain")}

    # --- Quick validation of base_url & params ---  # ← NEW
    base_url = getattr(BoldIdRequest, "base_url", "") or ""
    params = getattr(BoldIdRequest, "params", {}) or {}

    # If base_url looks empty or invalid, mark for retry and return
    if not base_url or not isinstance(base_url, str):
        BoldIdRequest.retry_count = getattr(BoldIdRequest, "retry_count", 0) + 1
        if BoldIdRequest.retry_count <= getattr(BoldIdRequest, "max_retries", 5):
            wait_seconds = 2 ** BoldIdRequest.retry_count
            BoldIdRequest.next_attempt = datetime.datetime.now() + datetime.timedelta(seconds=wait_seconds)
            tqdm.write(
                f"{datetime.datetime.now():%H:%M:%S}: POST aborted — invalid base_url. "
                f"Retry {BoldIdRequest.retry_count}/{BoldIdRequest.max_retries}; waiting {wait_seconds}s."
            )
        else:
            tqdm.write(f"{datetime.datetime.now():%H:%M:%S}: POST permanently failed — invalid base_url.")
            BoldIdRequest.next_attempt = None
        return BoldIdRequest

    # --- Pre-flight probe: DNS + TCP + TLS handshake (non-fatal) ---  # ← NEW
    try:
        import socket, ssl
        from urllib.parse import urlparse

        parsed = urlparse(base_url)
        host = parsed.hostname or parsed.path  # fallback if base_url is just host-like
        port = parsed.port or 443

        # DNS lookup
        try:
            socket.getaddrinfo(host, port)
            dns_ok = True
        except Exception as dns_e:
            dns_ok = False
            tqdm.write(f"{datetime.datetime.now():%H:%M:%S}: Preflight DNS lookup failed for {host}: {dns_e}")

        # TCP + optional TLS handshake (quick, short timeout)
        if dns_ok:
            try:
                sock = socket.create_connection((host, port), timeout=5)
                try:
                    # TLS handshake to ensure server accepts TLS
                    ctx = ssl.create_default_context()
                    ssock = ctx.wrap_socket(sock, server_hostname=host)
                    ssock.close()
                    tcp_tls_ok = True
                except Exception as tls_e:
                    tcp_tls_ok = False
                    tqdm.write(f"{datetime.datetime.now():%H:%M:%S}: Preflight TLS handshake failed for {host}:{port}: {tls_e}")
                    try:
                        sock.close()
                    except Exception:
                        pass
            except Exception as conn_e:
                tcp_tls_ok = False
                tqdm.write(f"{datetime.datetime.now():%H:%M:%S}: Preflight TCP connect failed for {host}:{port}: {conn_e}")
        else:
            tcp_tls_ok = False
    except Exception:
        # If preflight tooling isn't available, don't block the request — just continue
        tcp_tls_ok = False

    # --- Debug log: show what we're about to POST (size + first header) ---  # ← NEW
    try:
        payload_len = len(data.encode("utf-8"))
    except Exception:
        payload_len = len(data)

    log("POST", f"Submitting job (payload={payload_len} bytes, first_seq={fasta_header})")

    # --- Make the POST request with a prepared session and conservative retry strategy ---
    with requests_html.HTMLSession() as session:
        session.headers.update({
            "User-Agent": "BOLDigger3/ID-Engine",
            "Accept": "application/json,text/html;q=0.9,*/*;q=0.8",
            # Note: files= will set multipart/form-data automatically, Content-Type header not needed here
        })

        # keep urllib3 internal retries modest (we do our own retry/backoff)
        retry_strategy = Retry(total=3, backoff_factor=1)
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("https://", adapter)

        try:
            response = session.post(
                base_url,               # clean base URL (no query string embedded)
                params=params,         # params in the query string
                files=files,
                timeout=30,
            )
            # parse JSON result (may raise)
            result = json.loads(response.text)

        except Exception as e:
            # treat as transient by default and schedule retry/backoff
            BoldIdRequest.retry_count = getattr(BoldIdRequest, "retry_count", 0) + 1
            if BoldIdRequest.retry_count <= getattr(BoldIdRequest, "max_retries", 5):
                wait_seconds = 2 ** BoldIdRequest.retry_count
                BoldIdRequest.next_attempt = datetime.datetime.now() + datetime.timedelta(seconds=wait_seconds)
                log("RETRY", f"POST failed ({type(e).__name__}) retry {BoldIdRequest.retry_count}/{BoldIdRequest.max_retries} wait={wait_seconds}s")    
            else:
                tqdm.write(
                    f"{datetime.datetime.now():%H:%M:%S}: POST permanently failed ({type(e).__name__}: {e})."
                )
                BoldIdRequest.next_attempt = None
            return BoldIdRequest

    # --- Validate result structure ---
    if not isinstance(result, dict) or "sub_id" not in result:
        BoldIdRequest.retry_count = getattr(BoldIdRequest, "retry_count", 0) + 1
        if BoldIdRequest.retry_count <= getattr(BoldIdRequest, "max_retries", 5):
            wait_seconds = 2 ** BoldIdRequest.retry_count
            BoldIdRequest.next_attempt = datetime.datetime.now() + datetime.timedelta(seconds=wait_seconds)
            tqdm.write(
                f"{datetime.datetime.now():%H:%M:%S}: Invalid BOLD response. "
                f"Retry {BoldIdRequest.retry_count}/{BoldIdRequest.max_retries}; waiting {wait_seconds}s."
            )
        else:
            tqdm.write(f"{datetime.datetime.now():%H:%M:%S}: Invalid response permanent failure.")
            BoldIdRequest.next_attempt = None
        return BoldIdRequest

    # ---------------- SUCCESS ----------------
    BoldIdRequest.result_url = f"https://id.boldsystems.org/submission/results/{result['sub_id']}"
    BoldIdRequest.timestamp = datetime.datetime.now()
    BoldIdRequest.last_checked = None
    BoldIdRequest.next_attempt = None

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
    raw_lines = list(response.iter_lines())

    # Case 1 — No content returned at all
    if not raw_lines:
        tqdm.write(f"⚠ Warning: request {request_id} returned EMPTY response; saving empty file.")
        id_engine_result = pd.DataFrame()
        output_file = database_path.joinpath(
            "boldigger3_data", f"request_id_{request_id}_{fasta_name}.parquet.snappy"
        )
        id_engine_result.to_parquet(output_file)
        return
    # Case 2 — Try to decode JSON lines
    json_objects = []
    for line in raw_lines:
        try:
            json_objects.append(json.loads(line))
        except json.JSONDecodeError as e:
            tqdm.write(
                f"⚠ JSON decode error in request {request_id}: {e}. "
                f"Corrupt or partial data; saving empty result."
            )
            id_engine_result = pd.DataFrame()
            output_file = database_path.joinpath(
                "boldigger3_data", f"request_id_{request_id}_{fasta_name}.parquet.snappy"
            )
            id_engine_result.to_parquet(output_file)
            return


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


#                     return active_queue
def log_failed_request(req, key, fasta_order, project_directory, reason):
    """
    Log any request that permanently failed (POST or GET).
    Always safe: never references undefined attributes.

    Produces:
        failed_requests.log  (append mode)
        failed_fasta/        (FASTA sequences of failed entries)
    """

    log_path = project_directory.joinpath("failed_requests.log")
    failed_fasta_dir = project_directory.joinpath("failed_fasta")
    failed_fasta_dir.mkdir(exist_ok=True)

    # --- Extract available URLs safely ---
    base_url = getattr(req, "base_url", "N/A")
    result_url = getattr(req, "result_url", "N/A")

    # --- Extract sequence safely ---
    seq = None
    if hasattr(req, "query_data") and req.query_data:
        # req.query_data is a list of FASTA lines
        seq = "".join(req.query_data).strip()
    else:
        seq = "NO_SEQUENCE_FOUND"

    # --- Extract original FASTA header if possible ---
    fasta_header = fasta_order.get(key, "UNKNOWN_FASTA_ENTRY")

    # --- Build a failure classification line ---
    classification = "UNKNOWN"
    if "POST" in reason.upper():
        classification = "POST_FAILURE"
    elif "HTTP" in reason.upper():
        classification = "HTTP_ERROR"
    elif "TIMEOUT" in reason.upper():
        classification = "TIMEOUT"
    elif "JSON" in reason.upper():
        classification = "INVALID_JSON"
    elif "EMPTY" in reason.upper():
        classification = "EMPTY_RESPONSE"

    # === Write to failed_requests.log ===
    with open(log_path, "a", encoding="utf-8") as f:
        f.write("\n" + "="*80 + "\n")
        f.write(f"Timestamp: {datetime.datetime.now()}\n")
        f.write(f"Request ID: {key}\n")
        f.write(f"FASTA Header: {fasta_header}\n")
        f.write(f"Failure Reason: {reason}\n")
        f.write(f"Classification: {classification}\n")

        f.write("\n--- URLs ---\n")
        f.write(f"Submission URL (POST): {base_url}\n")
        f.write(f"Result URL (GET):      {result_url}\n")

        f.write("\n--- Retry Status ---\n")
        f.write(f"Retry Count: {getattr(req, 'retry_count', 'N/A')}\n")
        f.write(f"Max Retries: {getattr(req, 'max_retries', 'N/A')}\n")
        f.write(f"Next Attempt: {getattr(req, 'next_attempt', 'N/A')}\n")
        f.write(f"Timestamp First Seen: {getattr(req, 'timestamp', 'N/A')}\n")

        f.write("\n--- Query Data (FASTA) ---\n")
        f.write(seq + "\n")

    # Also save the failed FASTA sequence for re-submission
    failed_fasta_path = failed_fasta_dir.joinpath(f"{key}.fas")
    with open(failed_fasta_path, "w", encoding="utf-8") as f:
        f.write(seq)

    return

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

        # ---- Ensure retry/backoff fields exist ----
        if not hasattr(req, "retry_count"):
            req.retry_count = 0
        if not hasattr(req, "max_retries"):
            req.max_retries = 3
        if not hasattr(req, "next_attempt"):
            req.next_attempt = None
        if not hasattr(req, "timestamp"):
            req.timestamp = None

    with requests_html.HTMLSession() as session:
        completed = []

        for key, req in list(active_queue.items()):
            now = datetime.datetime.now()

            # ===============================
            # 1. TIMEOUT HANDLING
            # ===============================
            if req.timestamp is None:
                req.timestamp = now

            if now - req.timestamp > datetime.timedelta(minutes=15):
                log("TIMEOUT", f"Request {key} exceeded 15 min runtime")

                req.retry_count += 1

                if req.retry_count <= req.max_retries:
                    wait = 2 ** req.retry_count
                    req.next_attempt = now + datetime.timedelta(seconds=wait)
                    tqdm.write(
                        f"{now:%H:%M:%S}: Retrying {key} "
                        f"({req.retry_count}/{req.max_retries}) – waiting {wait}s."
                    )
                    retry_queue[key] = req
                else:
                    tqdm.write(f"{now:%H:%M:%S}: Request {key} failed permanently.")
                    log_failed_request(req, key, fasta_order, project_directory, "Timeout exceeded")

                active_queue.pop(key, None)
                continue

            # -------- Backoff check --------
            if req.next_attempt and now < req.next_attempt:
                continue

            # ===============================
            # 2. Respect minimum poll interval
            # ===============================
            if now - req.last_checked < datetime.timedelta(seconds=15):
                continue

            # Try requesting JSON
            try:
                response = session.get(req.result_url, timeout=30)
            except Exception as e:
                log("ERROR", f"Network error request {key}: {e}")
                req.retry_count += 1

                if req.retry_count <= req.max_retries:
                    wait = 2 ** req.retry_count
                    req.next_attempt = now + datetime.timedelta(seconds=wait)
                    retry_queue[key] = req
                    tqdm.write(f"{now:%H:%M:%S}: Network retry {req.retry_count}/{req.max_retries}, waiting {wait}s.")
                else:
                    tqdm.write(f"{now:%H:%M:%S}: Network failure permanent for {key}.")
                    log_failed_request(req, key, fasta_order, project_directory, f"Network error: {e}")
                    #retry_queue.pop(key, None)

                active_queue.pop(key, None)
                continue

            req.last_checked = now

            # ===============================
            # 3. Handle BOLD server responses
            # ===============================
            if response.status_code == 404:
                continue

            # Permanent client errors -> discard
            if response.status_code in {400, 403, 410}:
                tqdm.write(
                    f"{now:%H:%M:%S}: Permanent failure for {key} "
                    f"(HTTP {response.status_code}), discarding."
                )
                log_failed_request(req, key, fasta_order, project_directory, f"HTTP {response.status_code}")

                active_queue.pop(key, None)
                continue

            # Server errors -> retry with backoff
            if response.status_code >= 500:
                req.retry_count += 1

                if req.retry_count <= req.max_retries:
                    wait = 2 ** req.retry_count
                    req.next_attempt = now + datetime.timedelta(seconds=wait)

                    log("RETRY", f"HTTP {response.status_code} request {key} retry {req.retry_count}/{req.max_retries} wait={wait}s")

                    retry_queue[key] = req
                else:
                    log("ERROR", f"Server error permanent request {key}")
                    log_failed_request(req, key, fasta_order, project_directory,
                                       f"HTTP {response.status_code} (server error)")

                active_queue.pop(key, None)
                continue

            # ---- SUCCESS OR PARSE FAILURE ----
            if response.text.strip():
                try:
                    parse_and_save_data(
                        req,
                        response,
                        fasta_order,
                        key,
                        project_directory,
                        fasta_name,
                    )
                    log("GET", f"Request {key} completed")
                    completed.append(key)

                except Exception as e:
                
                    log("RETRY", f"JSON parse failure request {key}: {e}")
                    req.retry_count += 1

                    if req.retry_count <= req.max_retries:
                        wait = 2 ** req.retry_count
                        req.next_attempt = now + datetime.timedelta(seconds=wait)
                        retry_queue[key] = req
                        tqdm.write(
                            f"{now:%H:%M:%S}: Parse retry {req.retry_count}/{req.max_retries}, "
                            f"waiting {wait}s."
                        )
                    else:
                        tqdm.write(
                            f"{now:%H:%M:%S}: Request {key} failed permanently due to repeated invalid JSON."
                        )
                        log_failed_request(req, key, fasta_order, project_directory, "Repeated invalid JSON")

                    active_queue.pop(key, None)
                    continue

            else:
                log("RETRY", f"Empty response request {key}")

                req.retry_count += 1

                if req.retry_count <= req.max_retries:
                    wait = 2 ** req.retry_count
                    req.next_attempt = now + datetime.timedelta(seconds=wait)
                    retry_queue[key] = req
                else:
                    tqdm.write(
                        f"{now:%H:%M:%S}: Request {key} failed permanently due to empty response."
                    )
                    log_failed_request(req, key, fasta_order, project_directory, "Repeated empty response")

                active_queue.pop(key, None)
                continue


    # return active_queue
        # Remove completed jobs
        for key in completed:
            active_queue.pop(key, None)

        # -------- NEW: immediately requeue retry jobs if active queue has space --------
        while len(active_queue) < 4 and retry_queue:
            retry_key, retry_req = retry_queue.popitem()
            active_queue[retry_key] = retry_req

    return active_queue

def parquet_to_duckdb(project_directory, database_path):
    data_dir = project_directory.joinpath("boldigger3_data")
    parquet_files = list(data_dir.glob("request_id_*.parquet.snappy"))

    if not parquet_files:
        log("DB", "No parquet files to import")
        return

    con = duckdb.connect(database_path)

    # --- NEW: detect whether the table exists ---
    table_exists = False
    try:
        result = con.execute(
            """
            SELECT table_name 
            FROM information_schema.tables 
            WHERE table_name = 'id_engine_results'
            """
        ).fetchone()
        table_exists = result is not None
    except Exception:
        table_exists = False
    # --------------------------------------------

    try:
        if not table_exists:
            # First insert: create table from parquet files
            #parquet_list = ", ".join(f"'{str(f)}'" for f in parquet_files)
            parquet_paths = [str(f) for f in parquet_files]
            con.execute(
                f"""
                CREATE TABLE id_engine_results AS
                SELECT * FROM read_parquet([{parquet_list}])
                """
            )
        else:
            # Insert into existing table
            for file in parquet_files:
                try:
                    con.execute(
                        f"""
                        INSERT INTO id_engine_results
                        SELECT * FROM read_parquet('{file}')
                        """
                    )
                except duckdb.IOException as e:
                    log("DB", f"Skipped parquet {file.name} ({e})")

    finally:
        con.close()

    # Remove successfully ingested parquet files
    for file in parquet_files:
        if file.is_file():
            file.unlink()


def split_into_batches(fasta_dict, batch_size=5000):
    """Split fasta_dict into batches of max batch_size sequences."""
    items = list(fasta_dict.items())
    for i in range(0, len(items), batch_size):
        yield dict(items[i:i + batch_size])


def main(fasta_path: str, database: int, operating_mode: int) -> None:
    """Main function to run the BOLD identification engine."""

    log("INFO", "Reading input FASTA")

    fasta_dict, fasta_name, project_directory = parse_fasta(fasta_path)
    fasta_dict_order = {key: idx for idx, key in enumerate(fasta_dict.keys())}

    database_path = project_directory.joinpath("boldigger3_data", f"{fasta_name}.duckdb")
    download_queue_name = project_directory.joinpath(
        "boldigger3_data", f"{fasta_name}_download_queue.pkl"
    )

    data_dir = project_directory.joinpath("boldigger3_data")
    data_dir.mkdir(exist_ok=True)

    # REQUIRED BY BOLD API guidelines (6 seconds)
    submission_delay = 6

    # Check for prior downloads
    fasta_dict = already_downloaded(fasta_dict, database_path)

    if not fasta_dict:
        tqdm.write(f"{datetime.datetime.now():%H:%M:%S}: All data has already been downloaded.")
        return None

    # -------- NEW: split sequences into 5000-sequence batches --------
    batches = list(split_into_batches(fasta_dict, 5000))

    log("INFO", f"{len(batches)} batches detected (≤5000 sequences each)")


    # -------- NEW: process batches sequentially --------
    for batch_index, batch_dict in enumerate(batches, start=1):

        log("BATCH", f"Starting batch {batch_index}/{len(batches)} ({len(batch_dict)} sequences)")

        batch_start_time = datetime.datetime.now()

        download_queue = build_download_queue(batch_dict, database, operating_mode)
        download_queue["retry"] = OrderedDict()

        total_downloads = len(download_queue["waiting"])

        # ------------------- MAIN DOWNLOAD LOOP -------------------
        with tqdm(total=total_downloads, desc="Finished downloads") as pbar:
            while True:
                try:
                    if download_queue["waiting"] or download_queue["active"] or download_queue["retry"]:

                        # Move retry items back to waiting if nothing is waiting
                        if not download_queue["waiting"] and download_queue["retry"]:

                            log("QUEUE", f"Requeueing {len(download_queue['retry'])} retry jobs")

                            for key, req in list(download_queue["retry"].items()):
                                download_queue["waiting"][key] = req
                                del download_queue["retry"][key]

                        # Fill active queue up to four jobs
                        if len(download_queue["active"]) < 4 and download_queue["waiting"]:
                            request_id, req_obj = download_queue["waiting"].popitem(last=False)

                            log("QUEUE", f"Request {request_id} → ACTIVE")

                            download_queue["active"][request_id] = build_post_request(req_obj)

                            # REQUIRED BY BOLD API guidelines
                            time.sleep(submission_delay)

                        else:
                            # Update active queue and remove completed requests
                            before = len(download_queue["active"])
                            download_queue["active"] = download_json(
                                download_queue["active"],
                                fasta_dict_order,
                                project_directory,
                                fasta_name,
                                download_queue["retry"],
                            )
                            after = len(download_queue["active"])

                            # Only advance the progress bar for actual completed downloads
                            finished_now = before - after

                            if finished_now > 0:
                                pbar.update(finished_now)

                            if finished_now == 0:
                                submission_delay = min(submission_delay + 2, 30)
                            else:
                                submission_delay = max(6, submission_delay - 1)
                            # NEW: avoid spamming GET requests
                            time.sleep(1)    


                    else:
                        parquet_to_duckdb(project_directory, database_path)
                        raise DownloadFinished

                except DownloadFinished:
                    fasta_dict = already_downloaded(fasta_dict, database_path)


                    # ---------------- Minimum batch duration enforcement ----------------
                    min_batch_duration = datetime.timedelta(minutes=MIN_BATCH_MINUTES)
                    batch_elapsed = datetime.datetime.now() - batch_start_time
                    if batch_elapsed < min_batch_duration:
                        wait_time = (min_batch_duration - batch_elapsed).total_seconds()
                        log("BATCH", f"Batch {batch_index} finished early ({batch_elapsed}), waiting {int(wait_time)}s to enforce 30 min minimum")
                        time.sleep(wait_time)
                    else:
                        log("BATCH", f"Batch {batch_index} took {batch_elapsed}, no wait needed")
                    # -------------------------------------------------------------------------

                    if fasta_dict:
                        tqdm.write(f"{datetime.datetime.now():%H:%M:%S}: Requeuing incomplete downloads.")

                        download_queue = build_download_queue(fasta_dict, database, operating_mode)
                        download_queue["retry"] = OrderedDict()

                        total_downloads = len(download_queue["active"]) + len(download_queue["waiting"])
                        pbar.reset()
                        pbar.total = total_downloads
                        pbar.refresh()

                    else:
                        log("DONE", "All downloads finished successfully")

                        # os.remove(download_queue_name)
                        if os.path.exists(download_queue_name):
                            os.remove(download_queue_name)
                        break

