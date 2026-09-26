#!/usr/bin/env python3
"""Retrieve the already-submitted P6 NCBI BLAST URL-API RID once.

This recovery helper can only issue ``CMD=Get`` for the RID already recorded in
``ncbi_urlapi_status.json``.  It never submits or resubmits a BLAST job.  Every
HTTP response is written to the raw-data directory before its service status is
interpreted, and the compact status/receipt records the next legal poll time.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta, timezone
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
RAW_ROOT = PROJECT / "data" / "raw" / "p6_protein"
META_ROOT = PROJECT / "data" / "source_metadata" / "p6_protein"
STATUS_PATH = META_ROOT / "ncbi_urlapi_status.json"
SUBMISSION_PATH = META_ROOT / "ncbi_urlapi_submission.json"
GET_RECEIPT_PATH = META_ROOT / "ncbi_urlapi_get_receipt.json"
READY_RESULT_PATH = RAW_ROOT / "ncbi_urlapi_result.xml2"
EXPECTED_ENDPOINT = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
USER_AGENT = "ART-P6-recoverable-URLAPI/1.0"


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def iso(value: datetime) -> str:
    return value.isoformat(timespec="microseconds").replace("+00:00", "Z")


def parse_iso(value: str | None) -> datetime | None:
    if not value:
        return None
    return datetime.fromisoformat(value.replace("Z", "+00:00"))


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def sanitize_response_payload(payload: bytes) -> tuple[bytes, list[str]]:
    """Remove ephemeral NCBI session/account values from textual output."""

    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError:
        return payload, []
    fields = []
    patterns = [
        (
            "ncbi_sessionid_or_phid",
            r'(<meta\s+name="ncbi_(?:sessionid|phid)"\s+content=")[^"]+',
            r"\1[redacted]",
        ),
        ("myncbi_url_user", r"(MYNCBI%5FUSER%3D)\d+", r"\1REDACTED"),
        (
            "myncbi_hidden_user",
            r'(name="MYNCBI_USER"[^>]*value=")[^"]*',
            r"\1[redacted]",
        ),
    ]
    for name, pattern, replacement in patterns:
        text, count = re.subn(pattern, replacement, text, flags=re.IGNORECASE)
        if count:
            fields.append(name)
    return text.encode("utf-8"), fields


def safe_response_headers(headers: dict[str, str]) -> dict[str, str]:
    """Retain useful transport metadata without cookies or session identifiers."""

    blocked = {"set-cookie", "ncbi-phid", "ncbi-sid"}
    return {key: value for key, value in headers.items() if key.lower() not in blocked}


def classify(payload: bytes, http_status: int) -> tuple[str, str | None]:
    """Return a compact service state and any parsed XML root tag."""

    text = payload.decode("utf-8", errors="replace")
    status_match = re.search(r"\bStatus\s*=\s*(WAITING|READY|FAILED|UNKNOWN)\b", text, re.I)
    if status_match:
        qblast_status = status_match.group(1).upper()
        if qblast_status == "WAITING":
            return "WAITING", None
        if qblast_status == "FAILED":
            return "FAILED", None
        if qblast_status == "UNKNOWN":
            return "EXPIRED_OR_UNKNOWN", None
        # READY without result content is followed by a later Get according to
        # the legacy QBlast contract; preserve it as ready-without-payload.
        if qblast_status == "READY" and "ThereAreHits=yes" in text:
            return "READY_METADATA_ONLY", None
        if qblast_status == "READY" and "ThereAreHits=no" in text:
            return "READY_NO_HITS_METADATA_ONLY", None

    lower = text.lower()
    if "rid" in lower and ("expired" in lower or "not found" in lower or "unknown" in lower):
        return "EXPIRED_OR_UNKNOWN", None
    if http_status >= 400:
        return "HTTP_ERROR", None
    if not payload.strip():
        return "EMPTY_RESPONSE", None

    try:
        root = ET.fromstring(payload)
    except ET.ParseError:
        return "UNRECOGNIZED_RESPONSE", None

    root_tag = root.tag.rsplit("}", 1)[-1]
    if root_tag in {"BlastXML2", "BlastOutput2"} or root.find(".//BlastOutput2") is not None:
        return "READY", root_tag
    if root_tag.lower() in {"error", "errors"} or root.find(".//error") is not None:
        return "FAILED", root_tag
    return "UNRECOGNIZED_XML", root_tag


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--resume-authorized",
        action="store_true",
        help="Record that the user explicitly resumed P6 after the local pause.",
    )
    args = parser.parse_args()

    status = json.loads(STATUS_PATH.read_text(encoding="utf-8"))
    submission = json.loads(SUBMISSION_PATH.read_text(encoding="utf-8"))
    rid = str(status.get("rid") or "").strip()
    if not rid or rid != str(submission.get("rid") or "").strip():
        raise RuntimeError("status/submission RID mismatch; refusing retrieval")
    if status.get("terminal") or status.get("polling_complete"):
        raise RuntimeError(
            f"RID has terminal status {status.get('status')}; refusing any further Get"
        )
    if status.get("polling_paused") and not args.resume_authorized:
        raise RuntimeError("local P6 polling is paused; --resume-authorized is required")

    get_url = str(status.get("get_url") or "")
    parsed = urllib.parse.urlparse(get_url)
    query = urllib.parse.parse_qs(parsed.query)
    if (
        parsed.scheme != "https"
        or parsed.netloc != "blast.ncbi.nlm.nih.gov"
        or parsed.path != "/Blast.cgi"
        or query.get("CMD") != ["Get"]
        or query.get("RID") != [rid]
        or query.get("FORMAT_TYPE") != ["XML2"]
        or get_url.split("?", 1)[0] != EXPECTED_ENDPOINT
    ):
        raise RuntimeError("recorded Get URL is outside the frozen RID/format contract")

    now = utc_now()
    next_allowed = parse_iso(status.get("next_allowed_check_utc"))
    next_scheduled = parse_iso(status.get("next_scheduled_poll_utc"))
    guards = [value for value in (next_allowed, next_scheduled) if value is not None]
    effective_next = max(guards) if guards else None
    if effective_next and now < effective_next:
        remaining = (effective_next - now).total_seconds()
        raise RuntimeError(f"poll interval guard: wait another {remaining:.1f} seconds")

    requested_utc = now
    request = urllib.request.Request(get_url, headers={"User-Agent": USER_AGENT})
    response_error = None
    try:
        with urllib.request.urlopen(request, timeout=180) as response:
            payload = response.read()
            http_status = int(response.status)
            headers = dict(response.headers.items())
    except urllib.error.HTTPError as exc:
        payload = exc.read()
        http_status = int(exc.code)
        headers = dict(exc.headers.items()) if exc.headers else {}
        response_error = f"HTTPError: {exc}"
    except Exception as exc:
        payload = b""
        http_status = 0
        headers = {}
        response_error = f"{type(exc).__name__}: {exc}"

    received_utc = utc_now()
    original_response_bytes = len(payload)
    payload, sanitized_fields = sanitize_response_payload(payload)
    headers = safe_response_headers(headers)
    stamp = requested_utc.strftime("%Y%m%dT%H%M%S%fZ")
    raw_path = RAW_ROOT / f"ncbi_urlapi_get_{stamp}.response"
    poll_receipt_path = META_ROOT / f"ncbi_urlapi_get_{stamp}_receipt.json"
    raw_path.write_bytes(payload)

    service_status, xml_root = classify(payload, http_status)
    if response_error and service_status == "EMPTY_RESPONSE":
        service_status = "NETWORK_ERROR"
    if service_status == "READY":
        READY_RESULT_PATH.write_bytes(payload)

    minimum_interval = int(status.get("minimum_poll_interval_seconds") or 60)
    next_poll = requested_utc + timedelta(seconds=max(60, minimum_interval))
    scheduled_interval = int(status.get("scheduled_poll_interval_seconds") or 0)
    next_scheduled_poll = (
        requested_utc + timedelta(seconds=max(60, scheduled_interval))
        if service_status == "WAITING" and scheduled_interval
        else None
    )
    receipt = {
        "rid": rid,
        "request_method": "GET",
        "get_url": get_url,
        "requested_utc": iso(requested_utc),
        "received_utc": iso(received_utc),
        "http_status": http_status,
        "response_headers": headers,
        "response_bytes_original": original_response_bytes,
        "response_bytes_saved": len(payload),
        "response_body_sanitized_fields": sanitized_fields,
        "raw_response_path": raw_path.relative_to(PROJECT).as_posix(),
        "poll_receipt_path": poll_receipt_path.relative_to(PROJECT).as_posix(),
        "service_status": service_status,
        "xml_root": xml_root,
        "error": response_error,
        "minimum_poll_interval_seconds": minimum_interval,
        "next_allowed_check_utc": iso(next_poll),
        "scheduled_poll_interval_seconds": scheduled_interval or None,
        "next_scheduled_poll_utc": iso(next_scheduled_poll) if next_scheduled_poll else None,
        "submitted_utc": submission.get("submitted_utc"),
        "resubmitted": False,
        "relationship_to_old_cli_job": submission.get("relationship_to_old_cli_job"),
        "result_path_when_ready": (
            READY_RESULT_PATH.relative_to(PROJECT).as_posix() if service_status == "READY" else None
        ),
    }
    write_json(poll_receipt_path, receipt)
    write_json(GET_RECEIPT_PATH, receipt)

    status.update(
        {
            "status": service_status,
            "last_poll_utc": iso(requested_utc),
            "last_response_utc": iso(received_utc),
            "last_http_status": http_status,
            "last_response_bytes": len(payload),
            "last_raw_response_path": raw_path.relative_to(PROJECT).as_posix(),
            "next_allowed_check_utc": iso(next_poll),
            "next_scheduled_poll_utc": iso(next_scheduled_poll) if next_scheduled_poll else None,
            "get_receipt": GET_RECEIPT_PATH.relative_to(PROJECT).as_posix(),
            "last_poll_receipt": poll_receipt_path.relative_to(PROJECT).as_posix(),
            "poll_count_local": int(status.get("poll_count_local") or 0) + 1,
            "polling_paused": False,
            "resume_requires_explicit_user_instruction": False,
            "local_work_status": "P6_RESUMED_BY_USER",
            "resume_authorized_utc": iso(requested_utc) if args.resume_authorized else status.get("resume_authorized_utc"),
            "result_path_when_ready": (
                READY_RESULT_PATH.relative_to(PROJECT).as_posix()
                if service_status == "READY"
                else status.get("result_path_when_ready")
            ),
        }
    )
    write_json(STATUS_PATH, status)

    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    if service_status in {
        "READY",
        "WAITING",
        "READY_METADATA_ONLY",
        "READY_NO_HITS_METADATA_ONLY",
    }:
        return 0
    if service_status in {"EXPIRED_OR_UNKNOWN", "FAILED"}:
        return 20
    return 21


if __name__ == "__main__":
    sys.exit(main())
