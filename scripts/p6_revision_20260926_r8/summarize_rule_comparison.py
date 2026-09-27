#!/usr/bin/env python3
"""Summarize fixed component variants of the saved P2 correspondence rule.

This bounded diagnostic reads saved P2 tables.  It does not realign sequences,
change thresholds, or treat any rule variant as a truth or accuracy benchmark.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
P2 = ROOT / "data" / "processed" / "p2_comparison"
OUT = ROOT / "data" / "processed" / "p6_revision_20260926_r8"

BEST_HITS = P2 / "unit_free_best_hits.tsv"
UNITS = P2 / "unified_units.tsv"
COPIES = P2 / "unified_repeat_copies.tsv"
SUPPORTED = P2 / "supported_unit_correspondences.tsv"

TIE_TOLERANCE = 1e-12
BOUNDARY_TOLERANCE_NT = 10
P0_GROUP = "P0_seed7"
PB_GROUP = "PB50_near3"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: encode(row.get(field)) for field in fields})


def encode(value: object) -> object:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return value


def as_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def optional_int(value: str) -> int | None:
    return None if value == "" else int(value)


def relative(path: Path) -> str:
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return path.name


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="Repository root containing data/processed.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Destination for rule_comparison.tsv, rule_comparison_pairs.tsv, "
            "and rule_comparison_receipt.json. The default is the saved-data directory."
        ),
    )
    return parser.parse_args()


def canonical_pair(unit_a: dict[str, str], unit_b: dict[str, str]) -> tuple[dict[str, str], dict[str, str]]:
    if unit_a["analysis_subgroup"] == P0_GROUP and unit_b["analysis_subgroup"] == PB_GROUP:
        return unit_a, unit_b
    if unit_b["analysis_subgroup"] == P0_GROUP and unit_a["analysis_subgroup"] == PB_GROUP:
        return unit_b, unit_a
    raise AssertionError("Expected one primary unit from each comparison group")


def main() -> None:
    global ROOT, P2, OUT, BEST_HITS, UNITS, COPIES, SUPPORTED
    args = parse_args()
    ROOT = args.project_root.resolve()
    P2 = ROOT / "data" / "processed" / "p2_comparison"
    OUT = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else ROOT / "data" / "processed" / "p6_revision_20260926_r8"
    )
    BEST_HITS = P2 / "unit_free_best_hits.tsv"
    UNITS = P2 / "unified_units.tsv"
    COPIES = P2 / "unified_repeat_copies.tsv"
    SUPPORTED = P2 / "supported_unit_correspondences.tsv"

    copy_rows = read_tsv(COPIES)
    primary_copy_ids = {
        row["copy_id"] for row in copy_rows if as_bool(row["participates_in_seed_chain"])
    }
    unit_rows = read_tsv(UNITS)
    primary_units = {
        row["unit_id"]: row
        for row in unit_rows
        if row["left_copy_id"] in primary_copy_ids and row["right_copy_id"] in primary_copy_ids
    }
    if len(primary_copy_ids) != 47 or len(primary_units) != 37:
        raise AssertionError(
            f"Primary object count drift: copies={len(primary_copy_ids)}, units={len(primary_units)}"
        )

    directed = []
    for raw in read_tsv(BEST_HITS):
        if raw["query_unit_id"] not in primary_units or raw["best_target_unit_id"] not in primary_units:
            continue
        if raw["query_group"] == raw["target_group"]:
            continue
        row = dict(raw)
        row["best_is_unique_bool"] = as_bool(raw["best_is_unique"])
        row["reciprocal_unique_best_bool"] = as_bool(raw["reciprocal_unique_best"])
        row["left_offset"] = optional_int(raw["left_boundary_offset_to_best_target_nt"])
        row["right_offset"] = optional_int(raw["right_boundary_offset_to_best_target_nt"])
        row["sequence_only_pass"] = row["best_is_unique_bool"] and row["reciprocal_unique_best_bool"]
        row["left_boundary_pass"] = (
            row["sequence_only_pass"]
            and row["left_offset"] is not None
            and abs(row["left_offset"]) <= BOUNDARY_TOLERANCE_NT
        )
        row["right_boundary_pass"] = (
            row["sequence_only_pass"]
            and row["right_offset"] is not None
            and abs(row["right_offset"]) <= BOUNDARY_TOLERANCE_NT
        )
        row["both_boundaries_pass"] = row["left_boundary_pass"] and row["right_boundary_pass"]
        directed.append(row)

    if len(directed) != 147:
        raise AssertionError(f"Expected 147 primary cross-group directed selections; observed {len(directed)}")
    locus_pairs_screened = {
        tuple(sorted((primary_units[row["query_unit_id"]]["candidate_id"], row["target_candidate_id"])))
        for row in directed
    }
    if len(locus_pairs_screened) != 21:
        raise AssertionError(f"Expected 21 cross-group locus pairs; observed {len(locus_pairs_screened)}")

    by_pair: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in directed:
        key = tuple(sorted((str(row["query_unit_id"]), str(row["best_target_unit_id"]))))
        by_pair[key].append(row)

    pair_rows: list[dict[str, object]] = []
    for pair_key, group in sorted(by_pair.items()):
        sequence_rows = [row for row in group if row["sequence_only_pass"]]
        if not sequence_rows:
            continue
        if len(sequence_rows) != 2:
            raise AssertionError(
                f"Reciprocal pair {pair_key} should have two reciprocal directed rows; observed {len(sequence_rows)}"
            )
        unit_a = primary_units[pair_key[0]]
        unit_b = primary_units[pair_key[1]]
        p0, pb = canonical_pair(unit_a, unit_b)
        direction_evidence = []
        for row in sorted(sequence_rows, key=lambda item: str(item["query_unit_id"])):
            direction_evidence.append(
                "{query}->{target}|left={left}|right={right}|L={left_pass}|R={right_pass}|LR={both_pass}".format(
                    query=row["query_unit_id"],
                    target=row["best_target_unit_id"],
                    left="NA" if row["left_offset"] is None else row["left_offset"],
                    right="NA" if row["right_offset"] is None else row["right_offset"],
                    left_pass=str(bool(row["left_boundary_pass"])).lower(),
                    right_pass=str(bool(row["right_boundary_pass"])).lower(),
                    both_pass=str(bool(row["both_boundaries_pass"])).lower(),
                )
            )
        pair_rows.append(
            {
                "pair_id": "__".join((p0["unit_id"], pb["unit_id"])),
                "locus_pair_id": "__".join((p0["label"], pb["label"])),
                "p0_unit_id": p0["unit_id"],
                "p0_label": p0["label"],
                "p0_ordinal": int(p0["unit_ordinal_distal_to_proximal"]),
                "pb_unit_id": pb["unit_id"],
                "pb_label": pb["label"],
                "pb_ordinal": int(pb["unit_ordinal_distal_to_proximal"]),
                "reciprocal_directed_rows_n": len(sequence_rows),
                "sequence_only_pass": True,
                "sequence_plus_left_pass": any(bool(row["left_boundary_pass"]) for row in sequence_rows),
                "sequence_plus_right_pass": any(bool(row["right_boundary_pass"]) for row in sequence_rows),
                "sequence_plus_both_pass": any(bool(row["both_boundaries_pass"]) for row in sequence_rows),
                "component_calls_same_in_both_directions": len(
                    {
                        (
                            bool(row["left_boundary_pass"]),
                            bool(row["right_boundary_pass"]),
                            bool(row["both_boundaries_pass"]),
                        )
                        for row in sequence_rows
                    }
                )
                == 1,
                "direction_evidence": ";".join(direction_evidence),
            }
        )

    rules = (
        ("unique_reciprocal_sequence", "Unique reciprocal sequence best only", "sequence_only_pass"),
        ("sequence_plus_left_boundary", "Unique reciprocal sequence best plus left boundary", "sequence_plus_left_pass"),
        ("sequence_plus_right_boundary", "Unique reciprocal sequence best plus right boundary", "sequence_plus_right_pass"),
        ("sequence_plus_both_boundaries", "Unique reciprocal sequence best plus both boundaries", "sequence_plus_both_pass"),
    )
    summary_rows: list[dict[str, object]] = []
    for rule_order, (rule_id, label, field) in enumerate(rules, 1):
        passing = [row for row in pair_rows if row[field]]
        summary_rows.append(
            {
                "rule_order": rule_order,
                "rule_id": rule_id,
                "rule_label": label,
                "tie_tolerance": TIE_TOLERANCE,
                "boundary_tolerance_nt": BOUNDARY_TOLERANCE_NT,
                "nonredundant_pair_count": len(passing),
                "locus_pair_coverage_count": len({row["locus_pair_id"] for row in passing}),
                "primary_units": len(primary_units),
                "cross_group_locus_pairs_screened": len(locus_pairs_screened),
                "directed_free_best_rows_screened": len(directed),
            }
        )

    saved_full = {
        tuple(sorted((row["left_unit_id"], row["right_unit_id"])))
        for row in read_tsv(SUPPORTED)
        if row["relationship"] == "cross_group"
        and row["left_unit_id"] in primary_units
        and row["right_unit_id"] in primary_units
    }
    recomputed_full = {
        tuple(sorted((row["p0_unit_id"], row["pb_unit_id"])))
        for row in pair_rows
        if row["sequence_plus_both_pass"]
    }
    full_locus_coverage = {
        row["locus_pair_id"] for row in pair_rows if row["sequence_plus_both_pass"]
    }
    if recomputed_full != saved_full:
        missing = sorted(saved_full - recomputed_full)
        extra = sorted(recomputed_full - saved_full)
        raise AssertionError(f"Full-rule pair-set regression failed; missing={missing}, extra={extra}")
    if len(recomputed_full) != 35 or len(full_locus_coverage) != 18:
        raise AssertionError(
            f"Full-rule count regression failed: pairs={len(recomputed_full)}, loci={len(full_locus_coverage)}"
        )

    summary_path = OUT / "rule_comparison.tsv"
    pair_path = OUT / "rule_comparison_pairs.tsv"
    receipt_path = OUT / "rule_comparison_receipt.json"
    write_tsv(
        summary_path,
        summary_rows,
        [
            "rule_order",
            "rule_id",
            "rule_label",
            "tie_tolerance",
            "boundary_tolerance_nt",
            "nonredundant_pair_count",
            "locus_pair_coverage_count",
            "primary_units",
            "cross_group_locus_pairs_screened",
            "directed_free_best_rows_screened",
        ],
    )
    write_tsv(
        pair_path,
        pair_rows,
        [
            "pair_id",
            "locus_pair_id",
            "p0_unit_id",
            "p0_label",
            "p0_ordinal",
            "pb_unit_id",
            "pb_label",
            "pb_ordinal",
            "reciprocal_directed_rows_n",
            "sequence_only_pass",
            "sequence_plus_left_pass",
            "sequence_plus_right_pass",
            "sequence_plus_both_pass",
            "component_calls_same_in_both_directions",
            "direction_evidence",
        ],
    )
    boundary_contingency = {
        "left_pass_right_pass": sum(
            bool(row["sequence_plus_left_pass"]) and bool(row["sequence_plus_right_pass"])
            for row in pair_rows
        ),
        "left_pass_right_fail": sum(
            bool(row["sequence_plus_left_pass"]) and not bool(row["sequence_plus_right_pass"])
            for row in pair_rows
        ),
        "left_fail_right_pass": sum(
            not bool(row["sequence_plus_left_pass"]) and bool(row["sequence_plus_right_pass"])
            for row in pair_rows
        ),
        "left_fail_right_fail": sum(
            not bool(row["sequence_plus_left_pass"]) and not bool(row["sequence_plus_right_pass"])
            for row in pair_rows
        ),
    }
    receipt = {
        "analysis": "saved-data component comparison of the P2 sequence-and-boundary rule",
        "scope": "retrospective descriptive rule comparison; not a truth, accuracy, sensitivity, or superiority analysis",
        "inputs": [relative(path) for path in (BEST_HITS, UNITS, COPIES, SUPPORTED)],
        "parameters": {
            "primary_chain_copies": len(primary_copy_ids),
            "primary_chain_units": len(primary_units),
            "cross_group_locus_pairs": len(locus_pairs_screened),
            "directed_free_best_rows": len(directed),
            "tie_tolerance": TIE_TOLERANCE,
            "boundary_tolerance_nt": BOUNDARY_TOLERANCE_NT,
            "pair_collapse": "unordered unit pair; a component rule passes when either saved reciprocal direction passes, matching the original rowwise-then-deduplicate implementation",
            "new_alignment_or_threshold_search": False,
        },
        "rules": summary_rows,
        "boundary_component_contingency_among_61_sequence_pairs": boundary_contingency,
        "pairs_with_directionally_discordant_component_calls": sum(
            not bool(row["component_calls_same_in_both_directions"]) for row in pair_rows
        ),
        "full_rule_regression": {
            "expected_nonredundant_pairs": 35,
            "observed_nonredundant_pairs": len(recomputed_full),
            "expected_locus_pair_coverage": 18,
            "observed_locus_pair_coverage": len(full_locus_coverage),
            "saved_pair_set_exact_match": recomputed_full == saved_full,
            "status": "pass",
        },
        "outputs": [relative(summary_path), relative(pair_path), relative(receipt_path)],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    with receipt_path.open("w", encoding="utf-8") as handle:
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")

    for row in summary_rows:
        print(
            f"{row['rule_id']}\tpairs={row['nonredundant_pair_count']}\t"
            f"locus_pairs={row['locus_pair_coverage_count']}"
        )
    print("full_rule_regression\tpass")


if __name__ == "__main__":
    main()
