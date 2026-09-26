#!/usr/bin/env python3
"""Prepare alignment-column BPP data for the representative P6 figure.

This script does not refold RNA or rescore any comparison.  It maps the saved
ViennaRNA BPP edges for two fixed, already-supported unit pairs onto the same
affine global alignment used by the P2/P6 scoring code.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LOCAL_PACKAGES = ROOT / "tools" / "python312_packages"
if str(LOCAL_PACKAGES) not in sys.path:
    sys.path.insert(0, str(LOCAL_PACKAGES))

from Bio.Align import PairwiseAligner  # noqa: E402


INPUT_DIR = ROOT / "data" / "processed" / "p6_rna"
OUTPUT_DATA = INPUT_DIR / "representative_alignment_bpp.tsv"
OUTPUT_META = INPUT_DIR / "representative_alignment_metadata.tsv"

REPRESENTATIVES = [
    ("PALS_2|unit01", "PB50|unit01", "U1"),
    ("PALS_2|unit02", "PB50|unit02", "U2"),
]
CONTEXT = "two_repeat"
DISPLAY_PROBABILITY_FLOOR = 0.01


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, str(path))
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sequence_position_to_alignment_column(gapped: str):
    mapping = {}
    sequence_position = 0
    for column, base in enumerate(gapped):
        if base != "-":
            mapping[sequence_position] = column
            sequence_position += 1
    return mapping


def main():
    p2 = load_module(ROOT / "scripts" / "p2_analysis" / "run_unified_locus_analysis.py", "p2_plot_library")
    aligner = p2.configure_aligner(PairwiseAligner)

    objects = {
        (row["unit_id"], row["context"]): row
        for row in read_tsv(INPUT_DIR / "unit_sequence_objects.tsv")
    }
    pair_scores = {
        (row["p0_unit_id"], row["pb_unit_id"], row["context"]): row
        for row in read_tsv(INPUT_DIR / "unit_pair_structure_similarity.tsv")
        if row["supported_correspondence"] == "true"
    }

    edges_by_object = {}
    for row in read_tsv(INPUT_DIR / "unit_bpp_edges.tsv"):
        key = (row["unit_id"], row["context"])
        edges_by_object.setdefault(key, []).append(row)

    plot_rows = []
    metadata_rows = []
    for left_id, right_id, ordinal in REPRESENTATIVES:
        left = objects[(left_id, CONTEXT)]
        right = objects[(right_id, CONTEXT)]
        evidence = p2.alignment_evidence(aligner, left["sequence_dna"], right["sequence_dna"])
        left_map = sequence_position_to_alignment_column(evidence["target_gapped"])
        right_map = sequence_position_to_alignment_column(evidence["query_gapped"])
        score = pair_scores[(left_id, right_id, CONTEXT)]

        for side, object_id, gapped, position_map in [
            ("PALS_2", left_id, evidence["target_gapped"], left_map),
            ("PB50", right_id, evidence["query_gapped"], right_map),
        ]:
            for edge in edges_by_object[(object_id, CONTEXT)]:
                probability = float(edge["probability"])
                if probability < DISPLAY_PROBABILITY_FLOOR:
                    continue
                left_position = int(edge["left_position_1based"]) - 1
                right_position = int(edge["right_position_1based"]) - 1
                plot_rows.append({
                    "ordinal": ordinal,
                    "side": side,
                    "unit_id": object_id,
                    "alignment_length": len(gapped),
                    "left_alignment_column_1based": position_map[left_position] + 1,
                    "right_alignment_column_1based": position_map[right_position] + 1,
                    "pair_class": edge["pair_class"],
                    "probability": "{:.8f}".format(probability),
                })

        metadata_rows.append({
            "ordinal": ordinal,
            "left_unit_id": left_id,
            "right_unit_id": right_id,
            "context": CONTEXT,
            "alignment_length": len(evidence["target_gapped"]),
            "left_sequence_length": len(left["sequence_dna"]),
            "right_sequence_length": len(right["sequence_dna"]),
            "left_gapped_sequence": evidence["target_gapped"],
            "right_gapped_sequence": evidence["query_gapped"],
            "left_repeat_columns_1based": "{}-{};{}-{}".format(
                left_map[int(left["left_repeat_start_0based"])] + 1,
                left_map[int(left["left_repeat_end_exclusive_0based"]) - 1] + 1,
                left_map[int(left["right_repeat_start_0based"])] + 1,
                left_map[int(left["right_repeat_end_exclusive_0based"]) - 1] + 1,
            ),
            "right_repeat_columns_1based": "{}-{};{}-{}".format(
                right_map[int(right["left_repeat_start_0based"])] + 1,
                right_map[int(right["left_repeat_end_exclusive_0based"]) - 1] + 1,
                right_map[int(right["right_repeat_start_0based"])] + 1,
                right_map[int(right["right_repeat_end_exclusive_0based"]) - 1] + 1,
            ),
            "context_global_identity": score["context_global_identity"],
            "saved_trimmed_spacer_global_identity": score["saved_trimmed_spacer_global_identity"],
            "spacer_spacer_bpp_overlap": score["spacer_spacer_bpp_overlap"],
            "repeat_repeat_bpp_overlap": score["repeat_repeat_bpp_overlap"],
            "mixed_bpp_overlap": score["mixed_bpp_overlap"],
            "display_probability_floor": DISPLAY_PROBABILITY_FLOOR,
        })

    write_tsv(
        OUTPUT_DATA,
        plot_rows,
        [
            "ordinal", "side", "unit_id", "alignment_length",
            "left_alignment_column_1based", "right_alignment_column_1based",
            "pair_class", "probability",
        ],
    )
    write_tsv(
        OUTPUT_META,
        metadata_rows,
        [
            "ordinal", "left_unit_id", "right_unit_id", "context", "alignment_length",
            "left_sequence_length", "right_sequence_length", "left_gapped_sequence",
            "right_gapped_sequence", "left_repeat_columns_1based", "right_repeat_columns_1based",
            "context_global_identity", "saved_trimmed_spacer_global_identity",
            "spacer_spacer_bpp_overlap", "repeat_repeat_bpp_overlap", "mixed_bpp_overlap",
            "display_probability_floor",
        ],
    )
    print(json.dumps({
        "representatives": len(metadata_rows),
        "display_edges": len(plot_rows),
        "data": str(OUTPUT_DATA),
        "metadata": str(OUTPUT_META),
    }, indent=2))


if __name__ == "__main__":
    main()
