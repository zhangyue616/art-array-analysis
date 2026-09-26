#!/usr/bin/env python3
"""Submit one recoverable multi-query BLASTP job through NCBI URL API.

The script is intentionally single-use.  If either the raw submission response
or the submission receipt already exists, it exits before contacting NCBI.
It submits but never polls; a later controller heartbeat may retrieve by RID.
"""

from __future__ import annotations

import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
QUERY_PATH = ROOT / "data" / "raw" / "p6_protein" / "SA1_MarsHill_PB50_rt_queries.faa"
RAW_RESPONSE = ROOT / "data" / "raw" / "p6_protein" / "ncbi_urlapi_put_response.txt"
RECEIPT = ROOT / "data" / "source_metadata" / "p6_protein" / "ncbi_urlapi_submission.json"
STATUS = ROOT / "data" / "source_metadata" / "p6_protein" / "ncbi_urlapi_status.json"
ENDPOINT = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
OFFICIAL_URL_API_DOC = "https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html"
OFFICIAL_USAGE_DOC = "https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html"


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def iso(dt: datetime) -> str:
    return dt.isoformat().replace("+00:00", "Z")


def fasta_summary(text: str) -> list[dict[str, object]]:
    records = []
    header = None
    chunks: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append({"header": header, "length_aa": len("".join(chunks))})
            header = line[1:]
            chunks = []
        else:
            if header is None:
                raise ValueError("Sequence encountered before FASTA header")
            chunks.append(line)
    if header is not None:
        records.append({"header": header, "length_aa": len("".join(chunks))})
    return records


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def sanitize_response_text(text: str) -> str:
    """Remove ephemeral NCBI session/account identifiers but retain RID/RTOE."""
    text = re.sub(
        r'(<meta\s+name="ncbi_(?:sessionid|phid)"\s+content=")[^"]+',
        r"\1[redacted]",
        text,
        flags=re.IGNORECASE,
    )
    text = re.sub(r"(MYNCBI%5FUSER%3D)\d+", r"\1REDACTED", text, flags=re.IGNORECASE)
    text = re.sub(
        r'(name="MYNCBI_USER"[^>]*value=")[^"]*',
        r"\1[redacted]",
        text,
        flags=re.IGNORECASE,
    )
    return text


def main() -> int:
    if RECEIPT.exists() or RAW_RESPONSE.exists() or STATUS.exists():
        print("Submission artifact already exists; refusing to contact NCBI again.", file=sys.stderr)
        return 3

    query = QUERY_PATH.read_text(encoding="ascii")
    query_records = fasta_summary(query)
    if len(query_records) < 1 or len(query_records) > 3:
        raise ValueError(f"Expected 1-3 FASTA records, found {len(query_records)}")

    # Only currently documented Common URL API Put parameters are sent.  The
    # current parameter table does not list ENTREZ_QUERY, so this request is nr
    # without a server-side virus restriction; returned hits require post-filtering.
    parameters = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": "nr",
        "QUERY": query,
        "EXPECT": "1e-5",
        "HITLIST_SIZE": "200",
        "WORD_SIZE": "3",
        "GAPCOSTS": "11 2",
        "MATRIX": "BLOSUM62",
        "COMPOSITION_BASED_STATISTICS": "2",
        "FILTER": "F",
        "FORMAT_TYPE": "XML2",
    }
    submitted = utc_now()
    encoded = urllib.parse.urlencode(parameters).encode("ascii")
    request = urllib.request.Request(
        ENDPOINT,
        data=encoded,
        method="POST",
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "User-Agent": "ART-P6-recoverable-URLAPI/1.0",
        },
    )

    response_body = b""
    http_status = None
    response_headers: dict[str, str] = {}
    error = None
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            response_body = response.read()
            http_status = response.status
            response_headers = dict(response.headers.items())
    except urllib.error.HTTPError as exc:
        http_status = exc.code
        response_body = exc.read()
        response_headers = dict(exc.headers.items()) if exc.headers else {}
        error = f"HTTPError: {exc}"
    except Exception as exc:  # one attempt only; preserve the exact failure
        error = f"{type(exc).__name__}: {exc}"

    original_response_bytes = len(response_body)
    response_text = sanitize_response_text(response_body.decode("utf-8", errors="replace"))
    RAW_RESPONSE.parent.mkdir(parents=True, exist_ok=True)
    RAW_RESPONSE.write_text(response_text, encoding="utf-8")
    rid_match = re.search(r"^\s*RID\s*=\s*([A-Z0-9-]+)\s*$", response_text, flags=re.MULTILINE)
    rtoe_match = re.search(r"^\s*RTOE\s*=\s*(\d+)\s*$", response_text, flags=re.MULTILINE)
    rid = rid_match.group(1) if rid_match else None
    rtoe_seconds = int(rtoe_match.group(1)) if rtoe_match else None

    # NCBI's current developer guidance says not to poll an RID more often than
    # once per minute.  RTOE is the server estimate, so use the later bound.
    wait_seconds = max(60, rtoe_seconds or 0)
    next_check = submitted + timedelta(seconds=wait_seconds)
    status_name = "SUBMITTED_NOT_POLLED" if rid else "SUBMISSION_FAILED_NO_RETRY"

    public_parameters = {k: v for k, v in parameters.items() if k != "QUERY"}
    receipt = {
        "submitted_utc": iso(submitted),
        "attempt_count": 1,
        "resubmitted": False,
        "status": status_name,
        "http_status": http_status,
        "error": error,
        "rid": rid,
        "rtoe_seconds": rtoe_seconds,
        "next_allowed_check_utc": iso(next_check) if rid else None,
        "minimum_poll_interval_seconds": 60,
        "endpoint": ENDPOINT,
        "method": "POST",
        "parameters": public_parameters,
        "query_path": str(QUERY_PATH.relative_to(ROOT)).replace("\\", "/"),
        "query_records": query_records,
        "raw_response_path": str(RAW_RESPONSE.relative_to(ROOT)).replace("\\", "/"),
        "raw_response_bytes_original": original_response_bytes,
        "raw_response_bytes_saved": len(response_text.encode("utf-8")),
        "response_body_sanitized_fields": ["ncbi_sessionid", "ncbi_phid", "MYNCBI_USER"],
        "response_headers": {
            key: value
            for key, value in response_headers.items()
            if key.lower() in {"date", "server", "content-type", "content-length", "transfer-encoding"}
        },
        "official_documentation": [OFFICIAL_URL_API_DOC, OFFICIAL_USAGE_DOC],
        "server_side_taxon_filter": None,
        "taxon_filter_note": "Current Common URL API supported-parameter table does not list ENTREZ_QUERY. The request therefore uses unfiltered nr; any returned hit must be source-resolved and post-filtered before the unchanged array-blind RT/partner/neighborhood rules are applied.",
        "relationship_to_old_cli_job": "new authorized URL API submission; not recovery or retry of the prior CLI process that produced no RID",
        "retrieval": {
            "default_get_url": f"{ENDPOINT}?" + urllib.parse.urlencode({"CMD": "Get", "RID": rid}) if rid else None,
            "xml2_get_url": f"{ENDPOINT}?" + urllib.parse.urlencode({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML2"}) if rid else None,
            "allowed_action": "After next_allowed_check_utc, issue one Get; if not ready, wait at least 60 seconds before any later Get.",
        },
    }
    write_json(RECEIPT, receipt)
    write_json(
        STATUS,
        {
            "status": status_name,
            "rid": rid,
            "rtoe_seconds": rtoe_seconds,
            "submitted_utc": iso(submitted),
            "next_allowed_check_utc": iso(next_check) if rid else None,
            "minimum_poll_interval_seconds": 60,
            "get_url": receipt["retrieval"]["xml2_get_url"],
            "source_receipt": str(RECEIPT.relative_to(ROOT)).replace("\\", "/"),
            "result_path_when_ready": "data/raw/p6_protein/ncbi_urlapi_result.xml2",
            "postfilter_required": "unfiltered nr; resolve source taxonomy, then apply frozen RT+direct-partner+synteny+boundary rules without inspecting arrays",
        },
    )
    print(json.dumps({"status": status_name, "rid": rid, "rtoe_seconds": rtoe_seconds, "next_allowed_check_utc": receipt["next_allowed_check_utc"]}, indent=2))
    return 0 if rid else 2


if __name__ == "__main__":
    raise SystemExit(main())
