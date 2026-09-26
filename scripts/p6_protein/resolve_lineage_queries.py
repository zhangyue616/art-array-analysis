#!/usr/bin/env python3
"""Resolve cached bounded NCBI lineage ESearch IDs to versioned accessions.

The four ESearch responses were frozen before any upstream-array inspection.
This script makes at most one small official ESummary request for the union of
their nuccore IDs, caches the response, and otherwise works only from cache.
"""

from __future__ import annotations

import csv
import json
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw" / "p6_protein"
META_DIR = ROOT / "data" / "source_metadata" / "p6_protein"
OUT_DIR = ROOT / "data" / "processed" / "p6_extension"
INPUTS = {
    "LENTUS_COMPLETE_GENOME": RAW_DIR / "LENTUS_nuccore_esearch.json",
    "MADAWASKA_COMPLETE_GENOME": RAW_DIR / "MADAWASKA_nuccore_esearch.json",
    "MACHIAS_COMPLETE_GENOME": RAW_DIR / "MACHIAS_nuccore_esearch.json",
    "WALLMARK_COMPLETE_GENOME": RAW_DIR / "WALLMARK_nuccore_esearch.json",
}
SUMMARY_RAW = RAW_DIR / "lineage_nuccore_esummary.json"
OUTPUT_TSV = OUT_DIR / "lineage_nuccore_results.tsv"
RECEIPT = META_DIR / "lineage_nuccore_receipt.json"
FIXED = ROOT / "data" / "processed" / "p2_loci" / "fixed_10_loci.tsv"


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    query_data = {}
    all_ids = set()
    for query_id, path in INPUTS.items():
        payload = json.loads(path.read_text(encoding="utf-8"))
        result = payload["esearchresult"]
        ids = result["idlist"]
        query_data[query_id] = {
            "query_translation": result["querytranslation"],
            "reported_count": int(result["count"]),
            "ids": ids,
            "raw_path": str(path.relative_to(ROOT)).replace("\\", "/"),
        }
        all_ids.update(ids)

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(
        {"db": "nuccore", "id": ",".join(sorted(all_ids)), "retmode": "json"}
    )
    request_status = "cache_reused"
    if not SUMMARY_RAW.exists() or SUMMARY_RAW.stat().st_size == 0:
        req = urllib.request.Request(url, headers={"User-Agent": "ART-P6-bounded-lineage-resolution/1.0"})
        with urllib.request.urlopen(req, timeout=45) as response:
            data = response.read()
            http_status = response.status
        SUMMARY_RAW.write_bytes(data)
        request_status = "retrieved"
    else:
        http_status = 200

    summary_payload = json.loads(SUMMARY_RAW.read_text(encoding="utf-8"))
    result = summary_payload["result"]
    fixed_by_accession = {row["accession"]: row["label"] for row in read_tsv(FIXED)}
    query_ids_by_uid = {}
    for query_id, record in query_data.items():
        for uid in record["ids"]:
            query_ids_by_uid.setdefault(uid, []).append(query_id)

    rows = []
    for uid in sorted(all_ids, key=int):
        item = result[uid]
        accession = item.get("accessionversion", "")
        rows.append(
            {
                "nuccore_uid": uid,
                "accession_version": accession,
                "title": item.get("title", ""),
                "organism": item.get("organism", ""),
                "length_nt": item.get("slen", ""),
                "source_query_ids": ";".join(sorted(query_ids_by_uid[uid])),
                "is_fixed10_accession": str(accession in fixed_by_accession).lower(),
                "fixed10_label": fixed_by_accession.get(accession, ""),
                "source_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession or uid}",
            }
        )
    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_TSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    receipt = {
        "generated_utc": now_utc(),
        "request_status": request_status,
        "http_status": http_status,
        "request_url": url,
        "raw_esummary_path": str(SUMMARY_RAW.relative_to(ROOT)).replace("\\", "/"),
        "queries": query_data,
        "union_uid_count": len(all_ids),
        "resolved_accession_count": len(rows),
        "fixed10_accession_count": sum(row["is_fixed10_accession"] == "true" for row in rows),
        "new_accession_count": sum(row["is_fixed10_accession"] != "true" for row in rows),
        "output": str(OUTPUT_TSV.relative_to(ROOT)).replace("\\", "/"),
        "interpretation_limit": "Taxonomic complete-genome index scope only; no absence claim outside these four bounded queries.",
    }
    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: receipt[k] for k in ("union_uid_count", "resolved_accession_count", "fixed10_accession_count", "new_accession_count")}, indent=2))


if __name__ == "__main__":
    main()
