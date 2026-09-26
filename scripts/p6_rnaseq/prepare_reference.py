#!/usr/bin/env python3
"""Build the joint reference and frozen SA1 feature/window tables."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "python312_packages"))

from Bio import SeqIO  # type: ignore  # noqa: E402


PHAGE_GB = ROOT / "data" / "raw" / "ncbi" / "MW218148.1.gb"
HOST_GB = ROOT / "data" / "raw" / "p6_rnaseq" / "NZ_CP059679.1.gb"
OUT_DIR = ROOT / "data" / "processed" / "p6_rnaseq"
META_DIR = ROOT / "data" / "source_metadata" / "p6_rnaseq"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean(value: object) -> str:
    return str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    META_DIR.mkdir(parents=True, exist_ok=True)

    phage = SeqIO.read(PHAGE_GB, "genbank")
    host = SeqIO.read(HOST_GB, "genbank")
    if phage.id != "MW218148.1" or len(phage) != 260727:
        raise RuntimeError(f"Unexpected SA1 reference: {phage.id}, {len(phage)} nt")
    if host.id != "NZ_CP059679.1":
        raise RuntimeError(f"Unexpected host reference: {host.id}")

    reference_path = OUT_DIR / "MW218148.1_NZ_CP059679.1.fna"
    with reference_path.open("w", encoding="ascii", newline="\n") as handle:
        for record in (phage, host):
            handle.write(f">{record.id}\n")
            sequence = str(record.seq).upper()
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")

    features = []
    for ordinal, feature in enumerate((f for f in phage.features if f.type == "CDS"), 1):
        start0 = int(feature.location.start)
        end0 = int(feature.location.end)
        strand = int(feature.location.strand or 0)
        qualifiers = feature.qualifiers
        protein_id = qualifiers.get("protein_id", [""])[0]
        locus_tag = qualifiers.get("locus_tag", [""])[0]
        product = qualifiers.get("product", [""])[0]
        feature_id = protein_id or locus_tag or f"SA1_CDS_{ordinal:03d}"
        features.append(
            {
                "feature_ordinal": ordinal,
                "feature_id": clean(feature_id),
                "protein_id": clean(protein_id),
                "locus_tag": clean(locus_tag),
                "product": clean(product),
                "start0": start0,
                "end0": end0,
                "length_nt": end0 - start0,
                "strand": "+" if strand == 1 else "-" if strand == -1 else ".",
            }
        )
    if len(features) != 258:
        raise RuntimeError(f"Expected 258 SA1 CDS features, found {len(features)}")

    feature_path = OUT_DIR / "sa1_features.tsv"
    with feature_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(features[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(features)

    windows = [
        ("array_whole", "whole_array", 9447, 10646, "+", "published_SA1_array"),
        ("array_distal_residual", "array_partition", 9447, 9597, "+", "array_whole"),
        ("U1_core_to_core", "core_to_core", 9597, 9809, "+", "array_whole"),
        ("U2_core_to_core", "core_to_core", 9809, 9980, "+", "array_whole"),
        ("U3_core_to_core", "core_to_core", 9980, 10181, "+", "array_whole"),
        ("U4_core_to_core", "core_to_core", 10181, 10357, "+", "array_whole"),
        ("array_terminal_residual", "array_partition", 10357, 10646, "+", "array_whole"),
        ("U1_spacer", "repeat_excluded_spacer", 9615, 9805, "+", "U1_core_to_core"),
        ("U2_spacer", "repeat_excluded_spacer", 9827, 9976, "+", "U2_core_to_core"),
        ("U3_spacer", "repeat_excluded_spacer", 9998, 10177, "+", "U3_core_to_core"),
        ("U4_spacer", "repeat_excluded_spacer", 10199, 10353, "+", "U4_core_to_core"),
        ("repeat01_mask", "repeat_mask", 9593, 9615, "+", "SA1_copy01"),
        ("repeat02_mask", "repeat_mask", 9805, 9827, "+", "SA1_copy02"),
        ("repeat03_mask", "repeat_mask", 9976, 9998, "+", "SA1_copy03"),
        ("repeat04_mask", "repeat_mask", 10177, 10199, "+", "SA1_copy04"),
        ("repeat05_mask", "repeat_mask", 10353, 10375, "+", "SA1_copy05"),
        ("RT_CDS", "context_CDS", 10646, 12152, "+", "QPI16926.1"),
        ("partner_CDS", "context_CDS", 12153, 12951, "+", "QPI16927.1"),
    ]
    partition = [windows[index] for index in (1, 2, 3, 4, 5, 6)]
    if sum(end - start for _, _, start, end, _, _ in partition) != 1199:
        raise RuntimeError("Whole-array partition does not sum to 1,199 nt")
    if any(partition[index][3] != partition[index + 1][2] for index in range(5)):
        raise RuntimeError("Whole-array partition is not contiguous")

    window_path = OUT_DIR / "sa1_windows.tsv"
    with window_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = ["window_id", "window_class", "start0", "end0", "length_nt", "strand", "parent"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for name, kind, start, end, strand, parent in windows:
            writer.writerow(
                {
                    "window_id": name,
                    "window_class": kind,
                    "start0": start,
                    "end0": end,
                    "length_nt": end - start,
                    "strand": strand,
                    "parent": parent,
                }
            )

    # Independent interface check against the RNA-structure lane.
    interface_path = ROOT / "data" / "processed" / "p6_rna" / "unit_sequence_objects.tsv"
    interface_check = {"path": str(interface_path.relative_to(ROOT)), "status": "not_present"}
    if interface_path.exists():
        with interface_path.open(encoding="utf-8", newline="") as handle:
            rows = [
                row
                for row in csv.DictReader(handle, delimiter="\t")
                if row["label"] == "SA1" and row["context"] == "one_repeat"
            ]
        expected = {
            "SA1|unit01": (-1053, -841, 190),
            "SA1|unit02": (-841, -670, 149),
            "SA1|unit03": (-670, -469, 179),
            "SA1|unit04": (-469, -293, 154),
        }
        observed = {
            row["unit_id"]: (
                int(row["object_start_relative_to_RT"]),
                int(row["object_end_exclusive_relative_to_RT"]),
                int(row["spacer_length_nt"]),
            )
            for row in rows
        }
        interface_check["status"] = "match" if observed == expected else "mismatch"
        interface_check["expected"] = expected
        interface_check["observed"] = observed
        if observed != expected:
            raise RuntimeError(f"P6 RNA interface mismatch: {observed!r}")

    receipt = {
        "reference_ids": [phage.id, host.id],
        "reference_lengths_nt": {phage.id: len(phage), host.id: len(host)},
        "phage_cds_count": len(features),
        "source_files": {
            str(PHAGE_GB.relative_to(ROOT)): {"sha256": sha256(PHAGE_GB), "bytes": PHAGE_GB.stat().st_size},
            str(HOST_GB.relative_to(ROOT)): {"sha256": sha256(HOST_GB), "bytes": HOST_GB.stat().st_size},
        },
        "derived_files": {
            str(reference_path.relative_to(ROOT)): {"sha256": sha256(reference_path), "bytes": reference_path.stat().st_size},
            str(feature_path.relative_to(ROOT)): {"sha256": sha256(feature_path), "bytes": feature_path.stat().st_size},
            str(window_path.relative_to(ROOT)): {"sha256": sha256(window_path), "bytes": window_path.stat().st_size},
        },
        "p6_rna_interface_check": interface_check,
    }
    with (META_DIR / "reference_receipt.json").open("w", encoding="utf-8") as handle:
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
