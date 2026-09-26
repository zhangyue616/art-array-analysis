#!/usr/bin/env python3
"""Fetch versioned GenPept records for the fixed ten RTs and partners."""

from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
TABLE = PROJECT / "data" / "processed" / "p2_loci" / "fixed_10_loci.tsv"
RAW = PROJECT / "data" / "raw" / "p6_protein" / "fixed_proteins"
META = PROJECT / "data" / "source_metadata" / "p6_protein" / "fixed_protein_record_receipt.json"


def now():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(TABLE.open("r", encoding="utf-8"), delimiter="\t"))
    proteins = []
    for row in rows:
        proteins.extend(
            [
                (row["label"], "RT", row["rt_protein_id"]),
                (row["label"], "partner", row["partner_protein_id"]),
            ]
        )
    receipts = []
    for label, role, accession in proteins:
        params = urllib.parse.urlencode(
            {"db": "protein", "id": accession, "rettype": "gp", "retmode": "text"}
        )
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + params
        request = urllib.request.Request(url, headers={"User-Agent": "ART-P6-fixed-protein/1.0"})
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = response.read()
            status = response.status
        path = RAW / (accession + ".gp")
        path.write_bytes(payload)
        receipts.append(
            {
                "label": label,
                "role": role,
                "accession": accession,
                "url": url,
                "http_status": status,
                "bytes": len(payload),
                "raw_path": path.relative_to(PROJECT).as_posix(),
            }
        )
        time.sleep(0.36)
    META.parent.mkdir(parents=True, exist_ok=True)
    META.write_text(
        json.dumps(
            {
                "retrieved_utc": now(),
                "source": "NCBI Protein EFetch GenPept",
                "records": receipts,
                "record_count": len(receipts),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print("retrieved", len(receipts))


if __name__ == "__main__":
    main()
