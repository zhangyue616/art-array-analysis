#!/usr/bin/env python3
"""Post hoc exploratory diagnostics from the saved scientific tables.

This script does not rerun RNA folding, read processing, alignment, or the
original correspondence selection.  It only summarizes saved rows and extracts
a small, fixed coverage window around the SA1 array/RT boundary.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data" / "processed" / "p6_revision_20260925"

PAIR_PATH = ROOT / "data" / "processed" / "p6_rna" / "unit_pair_structure_similarity.tsv"
SUPPORTED_PATH = ROOT / "data" / "processed" / "p2_comparison" / "supported_unit_correspondences.tsv"
CONTRAST_PATH = ROOT / "data" / "processed" / "p6_rnaseq" / "aggregate" / "paired_contrasts.tsv"
PARTITION_PATH = ROOT / "data" / "processed" / "p6_rnaseq" / "aggregate" / "whole_array_partitions.tsv"
LIBRARY_PATH = ROOT / "data" / "processed" / "p6_rnaseq" / "aggregate" / "library_metrics.tsv"
WINDOW_PATH = ROOT / "data" / "processed" / "p6_rnaseq" / "sa1_windows.tsv"
FEATURE_PATH = ROOT / "data" / "processed" / "p6_rnaseq" / "sa1_features.tsv"
ENA_PATH = ROOT / "data" / "source_metadata" / "p6_rnaseq" / "ena_PRJNA836150_read_run.json"

IDENTITY_FIELDS = (
    "context_global_identity",
    "saved_core_to_core_global_identity",
    "saved_trimmed_spacer_global_identity",
)
BPP_FIELD = "spacer_spacer_bpp_overlap"
ABSOLUTE_SHARE_THRESHOLD = 0.05
BOUNDARY_1BASED = 10647
LOCAL_WINDOW_NT = 100


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        if not rows:
            raise ValueError(f"Cannot infer fields for empty table: {path}")
        fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def display_path(path: Path) -> str:
    """Return a stable project-relative path when possible."""
    try:
        return str(path.resolve().relative_to(ROOT.resolve())).replace("\\", "/")
    except ValueError:
        return str(path.resolve()).replace("\\", "/")


def read_and_validate_local_coverage(
    path: Path,
    partition_rows: list[dict],
) -> list[dict]:
    """Read the bounded 12-library coverage export used by the light package.

    The export must contain exactly the same 200 fixed positions for each of the
    12 libraries.  Its ``source_coverage_file`` field preserves provenance to
    the full saved vector; using this input does not imply that vector was read.
    """
    required_fields = {
        "run_accession",
        "sample_alias",
        "time_min",
        "culture",
        "boundary_1based",
        "position_1based",
        "position_relative_to_rt_start0",
        "region",
        "plus_fragment_depth",
        "plus_depth_per_million_retained_phage_fragments",
        "source_coverage_file",
    }
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        observed_fields = set(reader.fieldnames or [])
        missing_fields = sorted(required_fields - observed_fields)
        if missing_fields:
            raise ValueError(f"Local coverage input is missing fields: {missing_fields}")
        raw_rows = list(reader)

    expected_by_run = {row["run_accession"]: row for row in partition_rows}
    expected_positions = set(range(
        BOUNDARY_1BASED - LOCAL_WINDOW_NT,
        BOUNDARY_1BASED + LOCAL_WINDOW_NT,
    ))
    expected_total = len(expected_by_run) * len(expected_positions)
    if len(expected_by_run) != 12:
        raise AssertionError(f"Expected 12 libraries before local coverage validation; observed {len(expected_by_run)}")
    if len(raw_rows) != expected_total:
        raise AssertionError(f"Expected {expected_total} local coverage rows; observed {len(raw_rows)}")

    normalized_rows = []
    positions_by_run: dict[str, set[int]] = defaultdict(set)
    for raw in raw_rows:
        run_accession = raw["run_accession"]
        if run_accession not in expected_by_run:
            raise AssertionError(f"Unexpected library in local coverage input: {run_accession}")
        partition_row = expected_by_run[run_accession]
        position = int(raw["position_1based"])
        depth = int(raw["plus_fragment_depth"])
        normalized_depth = float(raw["plus_depth_per_million_retained_phage_fragments"])
        if depth < 0 or not math.isfinite(normalized_depth) or normalized_depth < 0:
            raise ValueError(f"Invalid local coverage depth for {run_accession} at {position}")
        if position not in expected_positions:
            raise AssertionError(f"Out-of-window position for {run_accession}: {position}")
        if position in positions_by_run[run_accession]:
            raise AssertionError(f"Duplicate local coverage position for {run_accession}: {position}")
        positions_by_run[run_accession].add(position)

        expected_region = "terminal_last_100nt" if position < BOUNDARY_1BASED else "RT_first_100nt"
        expected_relative = position - BOUNDARY_1BASED
        checks = {
            "sample_alias": (raw["sample_alias"], str(partition_row["sample_alias"])),
            "time_min": (int(raw["time_min"]), int(partition_row["time_min"])),
            "culture": (int(raw["culture"]), int(partition_row["culture"])),
            "boundary_1based": (int(raw["boundary_1based"]), BOUNDARY_1BASED),
            "position_relative_to_rt_start0": (
                int(raw["position_relative_to_rt_start0"]),
                expected_relative,
            ),
            "region": (raw["region"], expected_region),
        }
        for field, (observed, expected) in checks.items():
            if observed != expected:
                raise AssertionError(
                    f"Local coverage {field} mismatch for {run_accession} at {position}: "
                    f"observed {observed!r}, expected {expected!r}"
                )
        if not raw["source_coverage_file"].strip():
            raise ValueError(f"Missing source_coverage_file for {run_accession} at {position}")

        normalized_rows.append({
            "run_accession": run_accession,
            "sample_alias": partition_row["sample_alias"],
            "time_min": int(partition_row["time_min"]),
            "culture": int(partition_row["culture"]),
            "boundary_1based": BOUNDARY_1BASED,
            "position_1based": position,
            "position_relative_to_rt_start0": expected_relative,
            "region": expected_region,
            "plus_fragment_depth": depth,
            "plus_depth_per_million_retained_phage_fragments": normalized_depth,
            "source_coverage_file": raw["source_coverage_file"],
        })

    observed_runs = set(positions_by_run)
    if observed_runs != set(expected_by_run):
        missing = sorted(set(expected_by_run) - observed_runs)
        extra = sorted(observed_runs - set(expected_by_run))
        raise AssertionError(f"Local coverage library mismatch; missing={missing}, extra={extra}")
    for run_accession, observed_positions in positions_by_run.items():
        if observed_positions != expected_positions:
            missing = sorted(expected_positions - observed_positions)
            extra = sorted(observed_positions - expected_positions)
            raise AssertionError(
                f"Local coverage coordinate mismatch for {run_accession}; missing={missing}, extra={extra}"
            )

    order = {
        row["run_accession"]: (int(row["time_min"]), int(row["culture"]))
        for row in partition_rows
    }
    normalized_rows.sort(key=lambda row: (order[row["run_accession"]], row["position_1based"]))
    return normalized_rows


def as_bool(value: object) -> bool:
    return str(value).strip().lower() == "true"


def median(values) -> float | None:
    data = list(values)
    return statistics.median(data) if data else None


def mean(values) -> float | None:
    data = list(values)
    return statistics.mean(data) if data else None


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    cursor = 0
    while cursor < len(order):
        end = cursor + 1
        while end < len(order) and values[order[end]] == values[order[cursor]]:
            end += 1
        rank = ((cursor + 1) + end) / 2.0
        for position in range(cursor, end):
            ranks[order[position]] = rank
        cursor = end
    return ranks


def pearson(left: list[float], right: list[float]) -> float | None:
    if len(left) != len(right) or len(left) < 2:
        return None
    left_mean = statistics.mean(left)
    right_mean = statistics.mean(right)
    numerator = sum((a - left_mean) * (b - right_mean) for a, b in zip(left, right))
    left_ss = sum((a - left_mean) ** 2 for a in left)
    right_ss = sum((b - right_mean) ** 2 for b in right)
    denominator = math.sqrt(left_ss * right_ss)
    return numerator / denominator if denominator else None


def spearman(left: list[float], right: list[float]) -> float | None:
    return pearson(average_ranks(left), average_ranks(right))


def encode(value: object) -> object:
    if isinstance(value, float):
        return round(value, 12)
    return value


def s1_identity_diagnostics() -> dict:
    raw_rows = [row for row in read_tsv(PAIR_PATH) if row["context"] == "two_repeat"]
    rows = []
    for source_row, raw in enumerate(raw_rows, 1):
        row = dict(raw)
        row["source_row_two_repeat"] = source_row
        row["supported"] = as_bool(raw["supported_correspondence"])
        for field in (*IDENTITY_FIELDS, BPP_FIELD):
            row[field] = float(raw[field])
        rows.append(row)

    supported = [row for row in rows if row["supported"]]
    other = [row for row in rows if not row["supported"]]
    if (len(rows), len(supported), len(other)) != (252, 35, 217):
        raise AssertionError("Expected two_repeat 252 = 35 supported + 217 other")

    by_locus: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        by_locus[row["locus_pair_id"]].append(row)

    global_ranges = {
        field: (min(row[field] for row in supported), max(row[field] for row in supported))
        for field in IDENTITY_FIELDS
    }
    local_ranges: dict[str, dict[str, tuple[float, float] | None]] = {}
    for locus_pair_id, group in by_locus.items():
        locus_supported = [row for row in group if row["supported"]]
        local_ranges[locus_pair_id] = {
            field: (
                (min(row[field] for row in locus_supported), max(row[field] for row in locus_supported))
                if locus_supported else None
            )
            for field in IDENTITY_FIELDS
        }

    pair_audit_rows = []
    for row in rows:
        audit = {
            "source_row_two_repeat": row["source_row_two_repeat"],
            "locus_pair_id": row["locus_pair_id"],
            "p0_unit_id": row["p0_unit_id"],
            "p0_label": row["p0_label"],
            "p0_ordinal": row["p0_ordinal"],
            "pb_unit_id": row["pb_unit_id"],
            "pb_label": row["pb_label"],
            "pb_ordinal": row["pb_ordinal"],
            "supported_correspondence": str(row["supported"]).lower(),
            "same_ordinal": row["same_ordinal"],
            BPP_FIELD: row[BPP_FIELD],
        }
        for field in IDENTITY_FIELDS:
            global_low, global_high = global_ranges[field]
            local = local_ranges[row["locus_pair_id"]][field]
            audit[field] = row[field]
            audit[f"other_in_global_supported_range__{field}"] = str(
                (not row["supported"]) and global_low <= row[field] <= global_high
            ).lower()
            audit[f"other_in_locus_local_supported_range__{field}"] = (
                ""
                if local is None
                else str((not row["supported"]) and local[0] <= row[field] <= local[1]).lower()
            )
        pair_audit_rows.append(audit)

    pair_fields = [
        "source_row_two_repeat", "locus_pair_id", "p0_unit_id", "p0_label", "p0_ordinal",
        "pb_unit_id", "pb_label", "pb_ordinal", "supported_correspondence", "same_ordinal",
        *IDENTITY_FIELDS, BPP_FIELD,
        *[f"other_in_global_supported_range__{field}" for field in IDENTITY_FIELDS],
        *[f"other_in_locus_local_supported_range__{field}" for field in IDENTITY_FIELDS],
    ]
    write_tsv(OUT / "s1_two_repeat_pair_audit.tsv", pair_audit_rows, pair_fields)

    locus_rows = []
    summary_rows = []
    for field in IDENTITY_FIELDS:
        for range_scope in ("global_supported_min_max", "locus_local_supported_min_max"):
            locus_effects = []
            unrestricted_locus_effects = []
            range_other_all_loci = []
            range_other_supported_loci = []
            supported_in_evaluable_loci = []
            excluded_supported_loci = []
            for locus_pair_id, group in sorted(by_locus.items()):
                locus_supported = [row for row in group if row["supported"]]
                locus_all_other = [row for row in group if not row["supported"]]
                interval = global_ranges[field] if range_scope == "global_supported_min_max" else local_ranges[locus_pair_id][field]
                if interval is None:
                    locus_other = []
                else:
                    locus_other = [
                        row for row in group
                        if (not row["supported"]) and interval[0] <= row[field] <= interval[1]
                    ]
                range_other_all_loci.extend(locus_other)
                if locus_supported:
                    range_other_supported_loci.extend(locus_other)
                evaluable = bool(locus_supported and locus_other)
                effect = (
                    median(row[BPP_FIELD] for row in locus_supported)
                    - median(row[BPP_FIELD] for row in locus_other)
                    if evaluable else None
                )
                if evaluable:
                    locus_effects.append(effect)
                    supported_in_evaluable_loci.extend(locus_supported)
                elif locus_supported:
                    excluded_supported_loci.append(locus_pair_id)
                classification = (
                    "no_supported_correspondence"
                    if not locus_supported
                    else "evaluable"
                    if evaluable
                    else "supported_present_but_no_range_restricted_other"
                )
                unrestricted_effect = (
                    median(row[BPP_FIELD] for row in locus_supported)
                    - median(row[BPP_FIELD] for row in locus_all_other)
                    if locus_supported else None
                )
                if unrestricted_effect is not None:
                    unrestricted_locus_effects.append(unrestricted_effect)
                locus_rows.append({
                    "identity_field": field,
                    "range_scope": range_scope,
                    "locus_pair_id": locus_pair_id,
                    "p0_label": group[0]["p0_label"],
                    "pb_label": group[0]["pb_label"],
                    "supported_n": len(locus_supported),
                    "all_other_n": len(locus_all_other),
                    "range_restricted_other_n": len(locus_other),
                    "range_low": "" if interval is None else interval[0],
                    "range_high": "" if interval is None else interval[1],
                    "range_is_degenerate_single_value": "" if interval is None else str(interval[0] == interval[1]).lower(),
                    "classification": classification,
                    "supported_bpp_median": "" if not locus_supported else median(row[BPP_FIELD] for row in locus_supported),
                    "all_other_bpp_median": "" if not locus_all_other else median(row[BPP_FIELD] for row in locus_all_other),
                    "supported_minus_all_other_bpp_median": "" if unrestricted_effect is None else unrestricted_effect,
                    "range_restricted_other_bpp_median": "" if not locus_other else median(row[BPP_FIELD] for row in locus_other),
                    "supported_minus_range_other_bpp_median": "" if effect is None else effect,
                    "effect_positive": "" if effect is None else str(effect > 0).lower(),
                })

            global_low, global_high = global_ranges[field]
            summary_rows.append({
                "identity_field": field,
                "range_scope": range_scope,
                "spearman_all_252": spearman([row[field] for row in rows], [row[BPP_FIELD] for row in rows]),
                "spearman_supported_35": spearman([row[field] for row in supported], [row[BPP_FIELD] for row in supported]),
                "global_supported_min": global_low,
                "global_supported_max": global_high,
                "supported_n_all": len(supported),
                "supported_identity_median_all": median(row[field] for row in supported),
                "supported_bpp_median_all": median(row[BPP_FIELD] for row in supported),
                "range_restricted_other_n_all_loci": len(range_other_all_loci),
                "range_restricted_other_identity_median_all_loci": median(row[field] for row in range_other_all_loci),
                "range_restricted_other_bpp_median_all_loci": median(row[BPP_FIELD] for row in range_other_all_loci),
                "range_restricted_other_n_at_supported_loci": len(range_other_supported_loci),
                "supported_loci_n": sum(any(row["supported"] for row in group) for group in by_locus.values()),
                "evaluable_supported_loci_n": len(locus_effects),
                "excluded_supported_loci_n": len(excluded_supported_loci),
                "excluded_supported_loci": ";".join(excluded_supported_loci),
                "supported_pairs_in_evaluable_loci_n": len(supported_in_evaluable_loci),
                "supported_identity_median_evaluable_loci": median(row[field] for row in supported_in_evaluable_loci),
                "supported_bpp_median_evaluable_loci": median(row[BPP_FIELD] for row in supported_in_evaluable_loci),
                "range_restricted_other_identity_median_at_supported_loci": median(
                    row[field] for row in range_other_supported_loci
                ),
                "range_restricted_other_bpp_median_at_supported_loci": median(
                    row[BPP_FIELD] for row in range_other_supported_loci
                ),
                "unrestricted_positive_locus_effect_n": sum(effect > 0 for effect in unrestricted_locus_effects),
                "unrestricted_locus_effect_median": median(unrestricted_locus_effects),
                "positive_locus_effect_n": sum(effect > 0 for effect in locus_effects),
                "locus_effect_median": median(locus_effects),
            })

    write_tsv(OUT / "s1_identity_range_summary.tsv", summary_rows)
    write_tsv(OUT / "s1_locus_pair_identity_range_diagnostic.tsv", locus_rows)

    identity_correlation_rows = []
    for index, left in enumerate(IDENTITY_FIELDS):
        for right in IDENTITY_FIELDS[index + 1:]:
            identity_correlation_rows.append({
                "identity_field_left": left,
                "identity_field_right": right,
                "spearman_all_252": spearman([row[left] for row in rows], [row[right] for row in rows]),
            })
    write_tsv(OUT / "s1_identity_metric_correlations.tsv", identity_correlation_rows)

    supported_rows = [row for row in read_tsv(SUPPORTED_PATH) if row["relationship"] == "cross_group"]
    if len(supported_rows) != 35:
        raise AssertionError("Expected 35 cross-group supported pairs")
    write_tsv(OUT / "s3_cross_group_supported_pairs.tsv", supported_rows)

    supported_locus_counts = Counter()
    p0_counts = Counter()
    for row in supported_rows:
        p0_label = row["left_label"] if row["left_group"] == "P0_seed7" else row["right_label"]
        pb_label = row["right_label"] if row["right_group"] == "PB50_near3" else row["left_label"]
        supported_locus_counts[(p0_label, pb_label)] += 1
        p0_counts[p0_label] += 1
    locus_count_rows = []
    for locus_pair_id, group in sorted(by_locus.items()):
        p0_label = group[0]["p0_label"]
        pb_label = group[0]["pb_label"]
        count = supported_locus_counts[(p0_label, pb_label)]
        locus_count_rows.append({
            "locus_pair_id": locus_pair_id,
            "p0_label": p0_label,
            "pb_label": pb_label,
            "supported_pair_count": count,
            "has_supported_pair": str(count > 0).lower(),
        })
    write_tsv(OUT / "s3_cross_group_supported_by_locus_pair.tsv", locus_count_rows)

    return {
        "population": {
            "context": "two_repeat",
            "all_combinations_n": len(rows),
            "supported_n": len(supported),
            "other_n": len(other),
            "locus_pairs_n": len(by_locus),
            "unit_combinations_per_locus_pair": sorted({len(group) for group in by_locus.values()}),
        },
        "identity_range_summary": [{key: encode(value) for key, value in row.items()} for row in summary_rows],
        "identity_metric_spearman": [{key: encode(value) for key, value in row.items()} for row in identity_correlation_rows],
        "locus_local_range_sensitivity": {
            "global_supported_min_max_context_result": "9/18 positive; no supported locus pair excluded",
            "locus_local_supported_min_max_context_result": (
                "4/15 positive; MarsHill__AH12, MarsHill__Machias, and MarsHill__PB50 excluded because "
                "each has one supported row, so the local interval is a single identity value with no other row at exactly that value"
            ),
        },
        "selection_boundary": (
            "The saved supported set was selected by reciprocal sequence-best and boundary criteria, not a single fixed "
            "global identity cutoff. The observed supported identity range is a post hoc exploratory diagnostic, not the selection rule."
        ),
        "dependence_boundary": (
            "Rows share loci and units, and the three identity fields are correlated measurements on the same 252 combinations. "
            "They are not independent validations or biological replicates."
        ),
        "supported_cross_group": {
            "n": len(supported_rows),
            "p0_locus_counts": dict(sorted(p0_counts.items())),
            "zero_supported_locus_pairs": [row["locus_pair_id"] for row in locus_count_rows if row["supported_pair_count"] == 0],
        },
    }


def s2_temporal_diagnostics() -> dict:
    rows = [
        row for row in read_tsv(CONTRAST_PATH)
        if row["endpoint"] == "spacer_midpoint_density_share"
    ]
    if len(rows) != 36:
        raise AssertionError("Expected 36 primary-endpoint culture contrasts")
    audit_rows = []
    for row in rows:
        start = float(row["start_value"])
        end = float(row["end_value"])
        delta = float(row["delta"])
        audit_rows.append({
            "comparison": row["comparison"],
            "comparison_role": row["comparison_role"],
            "unit": row["unit"],
            "culture": int(row["culture"]),
            "start_time_min": int(row["start_time_min"]),
            "end_time_min": int(row["end_time_min"]),
            "start_share": start,
            "end_share": end,
            "start_percentage_points": 100 * start,
            "end_percentage_points": 100 * end,
            "delta_share": delta,
            "delta_percentage_points": 100 * delta,
            "relative_change_percent": 100 * delta / start if start else "",
            "absolute_delta_percentage_points": 100 * abs(delta),
            "fixed_threshold_share": ABSOLUTE_SHARE_THRESHOLD,
            "fixed_threshold_percentage_points": 100 * ABSOLUTE_SHARE_THRESHOLD,
            "decrease_branch_reachable_from_start": str(start >= ABSOLUTE_SHARE_THRESHOLD).lower(),
            "increase_branch_reachable_from_start": str(start <= 1 - ABSOLUTE_SHARE_THRESHOLD).lower(),
            "relative_decrease_required_for_5pp_percent": 100 * ABSOLUTE_SHARE_THRESHOLD / start if start else "",
            "relative_increase_required_for_5pp_percent": 100 * ABSOLUTE_SHARE_THRESHOLD / start if start else "",
            "both_libraries_estimable": row["both_libraries_estimable"],
            "pair_status": row["pair_status"],
        })
    audit_rows.sort(key=lambda row: (
        {"5_to_55": 0, "5_to_15": 1, "15_to_55": 2}[row["comparison"]],
        row["unit"], row["culture"],
    ))
    write_tsv(OUT / "s2_temporal_absolute_relative_changes.tsv", audit_rows)

    summary_rows = []
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in audit_rows:
        grouped[(row["comparison"], row["unit"])].append(row)
    for (comparison, unit), group in sorted(
        grouped.items(), key=lambda item: ({"5_to_55": 0, "5_to_15": 1, "15_to_55": 2}[item[0][0]], item[0][1])
    ):
        deltas = [row["delta_share"] for row in group]
        signs = [1 if value > 0 else -1 if value < 0 else 0 for value in deltas]
        same_nonzero_direction = len(set(signs)) == 1 and signs[0] != 0
        median_abs = median(abs(value) for value in deltas)
        summary_rows.append({
            "comparison": comparison,
            "comparison_role": group[0]["comparison_role"],
            "unit": unit,
            "culture_delta_percentage_points": ";".join(f"{row['delta_percentage_points']:.9f}" for row in group),
            "culture_relative_change_percent": ";".join(f"{row['relative_change_percent']:.9f}" for row in group),
            "all_three_same_nonzero_direction": str(same_nonzero_direction).lower(),
            "direction": "increase" if same_nonzero_direction and signs[0] > 0 else "decrease" if same_nonzero_direction else "mixed",
            "median_absolute_delta_share": median_abs,
            "median_absolute_delta_percentage_points": 100 * median_abs,
            "fixed_threshold_share": ABSOLUTE_SHARE_THRESHOLD,
            "passes_fixed_absolute_threshold": str(median_abs >= ABSOLUTE_SHARE_THRESHOLD).lower(),
            "passes_original_material_change_rule": str(
                same_nonzero_direction
                and median_abs >= ABSOLUTE_SHARE_THRESHOLD
                and all(as_bool(row["both_libraries_estimable"]) for row in group)
            ).lower(),
            "decrease_branch_reachable_all_cultures": str(
                all(row["start_share"] >= ABSOLUTE_SHARE_THRESHOLD for row in group)
            ).lower(),
        })
    write_tsv(OUT / "s2_temporal_contrast_summary.tsv", summary_rows)

    unreachable_u1_u4 = [
        row for row in audit_rows
        if row["unit"] in {"U1", "U4"} and row["decrease_branch_reachable_from_start"] != "false"
    ]
    if unreachable_u1_u4:
        raise AssertionError("Expected the 5-percentage-point decrease branch to be unreachable for every U1/U4 contrast")

    return {
        "endpoint": "repeat-excluded spacer midpoint-density share",
        "original_rule_unchanged": True,
        "absolute_threshold_share": ABSOLUTE_SHARE_THRESHOLD,
        "absolute_threshold_percentage_points": 100 * ABSOLUTE_SHARE_THRESHOLD,
        "rows_n": len(audit_rows),
        "contrast_summary": [{key: encode(value) for key, value in row.items()} for row in summary_rows],
        "u1_u4_decrease_branch": (
            "Unreachable from every culture-specific start value in all three prespecified contrasts because each start share is below 0.05."
        ),
        "interpretation": (
            "Relative percent change is reported beside, and does not replace, the prespecified absolute percentage-point endpoint."
        ),
    }


def s4_partition_and_boundary_diagnostics(local_coverage_input: Path | None = None) -> dict:
    partition_rows = read_tsv(PARTITION_PATH)
    library_rows = {row["run_accession"]: row for row in read_tsv(LIBRARY_PATH)}
    metadata_rows = {row["run_accession"]: row for row in json.loads(ENA_PATH.read_text(encoding="utf-8-sig"))}
    windows = {row["window_id"]: row for row in read_tsv(WINDOW_PATH)}
    features = {row["feature_id"]: row for row in read_tsv(FEATURE_PATH)}

    expected_coordinates = {
        "array_distal_residual": (9447, 9597),
        "array_terminal_residual": (10357, 10646),
        "RT_CDS": (10646, 12152),
    }
    for window_id in ("array_distal_residual", "array_terminal_residual"):
        observed = (int(windows[window_id]["start0"]), int(windows[window_id]["end0"]))
        if observed != expected_coordinates[window_id]:
            raise AssertionError(f"Coordinate mismatch for {window_id}: {observed}")
    rt = features["QPI16926.1"]
    if (int(rt["start0"]), int(rt["end0"])) != expected_coordinates["RT_CDS"]:
        raise AssertionError("RT coordinate mismatch")

    by_run: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in partition_rows:
        by_run[row["run_accession"]][row["partition"]] = row
    if len(by_run) != 12:
        raise AssertionError("Expected 12 whole-array partition libraries")

    partition_audit_rows = []
    for run_accession, group in by_run.items():
        first = next(iter(group.values()))
        total = sum(int(row["sense_midpoint_fragments"]) for row in group.values())
        library = library_rows[run_accession]
        if total != int(library["whole_array_sense_midpoint_fragments"]):
            raise AssertionError(f"Whole-array denominator mismatch for {run_accession}")
        unit_total = sum(int(group[unit]["sense_midpoint_fragments"]) for unit in ("U1", "U2", "U3", "U4"))
        partition_audit_rows.append({
            "run_accession": run_accession,
            "sample_alias": first["sample_alias"],
            "time_min": int(first["time_min"]),
            "culture": int(first["culture"]),
            "whole_array_midpoint_denominator": total,
            "spacer_midpoint_count_separate_denominator": int(library["spacer_sense_midpoint_fragments"]),
            "distal_start0": 9447,
            "distal_end0": 9597,
            "distal_midpoint_fragments": int(group["distal_residual"]["sense_midpoint_fragments"]),
            "distal_whole_array_midpoint_share": float(group["distal_residual"]["whole_array_midpoint_share"]),
            "u1_u4_midpoint_fragments": unit_total,
            "u1_u4_whole_array_midpoint_share": unit_total / total if total else "",
            "terminal_start0": 10357,
            "terminal_end0": 10646,
            "terminal_midpoint_fragments": int(group["terminal_residual"]["sense_midpoint_fragments"]),
            "terminal_whole_array_midpoint_share": float(group["terminal_residual"]["whole_array_midpoint_share"]),
            "terminal_fraction_of_residual_midpoint_fragments": (
                int(group["terminal_residual"]["sense_midpoint_fragments"])
                / (
                    int(group["distal_residual"]["sense_midpoint_fragments"])
                    + int(group["terminal_residual"]["sense_midpoint_fragments"])
                )
                if (
                    int(group["distal_residual"]["sense_midpoint_fragments"])
                    + int(group["terminal_residual"]["sense_midpoint_fragments"])
                ) else ""
            ),
        })
    partition_audit_rows.sort(key=lambda row: (row["time_min"], row["culture"]))
    write_tsv(OUT / "s4_whole_array_partitions_by_library.tsv", partition_audit_rows)

    local_start = BOUNDARY_1BASED - LOCAL_WINDOW_NT
    local_end = BOUNDARY_1BASED + LOCAL_WINDOW_NT - 1
    if local_coverage_input is None:
        coverage_input_mode = "saved_full_per_base_vectors"
        coverage_input_display = "data/processed/p6_rnaseq/runs/<run_accession>/locus_coverage.tsv.gz"
        coverage_rows = []
        for partition_row in partition_audit_rows:
            run_accession = partition_row["run_accession"]
            path = ROOT / "data" / "processed" / "p6_rnaseq" / "runs" / run_accession / "locus_coverage.tsv.gz"
            if not path.exists():
                raise FileNotFoundError(path)
            with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
                for raw in csv.DictReader(handle, delimiter="\t"):
                    position = int(raw["position_1based"])
                    if local_start <= position <= local_end:
                        coverage_rows.append({
                            "run_accession": run_accession,
                            "sample_alias": partition_row["sample_alias"],
                            "time_min": partition_row["time_min"],
                            "culture": partition_row["culture"],
                            "boundary_1based": BOUNDARY_1BASED,
                            "position_1based": position,
                            "position_relative_to_rt_start0": position - BOUNDARY_1BASED,
                            "region": "terminal_last_100nt" if position < BOUNDARY_1BASED else "RT_first_100nt",
                            "plus_fragment_depth": int(raw["plus_fragment_depth"]),
                            "plus_depth_per_million_retained_phage_fragments": float(
                                raw["plus_depth_per_million_retained_phage_fragments"]
                            ),
                            "source_coverage_file": display_path(path),
                        })
    else:
        local_coverage_input = local_coverage_input.resolve()
        if not local_coverage_input.exists():
            raise FileNotFoundError(local_coverage_input)
        coverage_input_mode = "validated_local_coverage_table"
        coverage_input_display = display_path(local_coverage_input)
        coverage_rows = read_and_validate_local_coverage(local_coverage_input, partition_audit_rows)

    coverage_by_run: dict[str, list[dict]] = defaultdict(list)
    for row in coverage_rows:
        coverage_by_run[row["run_accession"]].append(row)

    continuity_rows = []
    expected_positions = set(range(local_start, local_end + 1))
    for partition_row in partition_audit_rows:
        run_accession = partition_row["run_accession"]
        local = coverage_by_run[run_accession]
        if len(local) != 2 * LOCAL_WINDOW_NT:
            raise AssertionError(f"Expected {2 * LOCAL_WINDOW_NT} local positions for {run_accession}")
        by_position = {row["position_1based"]: row for row in local}
        if set(by_position) != expected_positions or len(by_position) != len(local):
            raise AssertionError(f"Expected one row for every fixed local position for {run_accession}")
        terminal = [row for row in local if row["region"] == "terminal_last_100nt"]
        rt_start = [row for row in local if row["region"] == "RT_first_100nt"]
        terminal_median = median(row["plus_fragment_depth"] for row in terminal)
        rt_median = median(row["plus_fragment_depth"] for row in rt_start)
        left_boundary_depth = by_position[BOUNDARY_1BASED - 1]["plus_fragment_depth"]
        right_boundary_depth = by_position[BOUNDARY_1BASED]["plus_fragment_depth"]
        if partition_row["time_min"] == 0:
            evidence_status = "zero_min_low_signal_background_only"
        elif left_boundary_depth > 0 and right_boundary_depth > 0:
            evidence_status = "adjacent_boundary_depth_present_descriptive_only"
        else:
            evidence_status = "boundary_depth_absent_or_one_sided"
        continuity_rows.append({
            "run_accession": run_accession,
            "sample_alias": partition_row["sample_alias"],
            "time_min": partition_row["time_min"],
            "culture": partition_row["culture"],
            "boundary_1based": BOUNDARY_1BASED,
            "terminal_window_1based": f"{local_start}-{BOUNDARY_1BASED - 1}",
            "rt_window_1based": f"{BOUNDARY_1BASED}-{local_end}",
            "terminal_last100_plus_depth_median": terminal_median,
            "rt_first100_plus_depth_median": rt_median,
            "rt_to_terminal_median_depth_ratio": rt_median / terminal_median if terminal_median else "",
            "terminal_last100_plus_depth_mean": mean(row["plus_fragment_depth"] for row in terminal),
            "rt_first100_plus_depth_mean": mean(row["plus_fragment_depth"] for row in rt_start),
            "terminal_end_position_depth": left_boundary_depth,
            "rt_start_position_depth": right_boundary_depth,
            "boundary_step_depth": right_boundary_depth - left_boundary_depth,
            "boundary_step_relative_to_larger_depth": (
                (right_boundary_depth - left_boundary_depth) / max(left_boundary_depth, right_boundary_depth)
                if max(left_boundary_depth, right_boundary_depth) else ""
            ),
            "terminal_last100_coverage_breadth_gt0": sum(row["plus_fragment_depth"] > 0 for row in terminal) / LOCAL_WINDOW_NT,
            "rt_first100_coverage_breadth_gt0": sum(row["plus_fragment_depth"] > 0 for row in rt_start) / LOCAL_WINDOW_NT,
            "evidence_status": evidence_status,
        })
    write_tsv(OUT / "s4_terminal_rt_local_coverage_100nt.tsv", coverage_rows)
    write_tsv(OUT / "s4_terminal_rt_continuity_summary.tsv", continuity_rows)

    zero_metadata_rows = []
    for row in partition_audit_rows:
        if row["time_min"] != 0:
            continue
        metadata = metadata_rows[row["run_accession"]]
        zero_metadata_rows.append({
            "run_accession": row["run_accession"],
            "sample_accession": metadata["sample_accession"],
            "sample_alias": metadata["sample_alias"],
            "library_name": metadata["library_name"],
            "experiment_title": metadata["experiment_title"],
            "whole_array_midpoint_denominator": row["whole_array_midpoint_denominator"],
            "spacer_midpoint_count_separate_denominator": row["spacer_midpoint_count_separate_denominator"],
            "interpretation": "0-min baseline/background; B alias alone does not establish pre-phage sampling",
        })
    write_tsv(OUT / "s4_zero_min_metadata_and_denominators.tsv", zero_metadata_rows)

    postinfection_continuity = [row for row in continuity_rows if row["time_min"] > 0]
    return {
        "coordinates_0based_half_open": {
            "distal_residual": [9447, 9597],
            "terminal_residual": [10357, 10646],
            "rt_cds": [10646, 12152],
        },
        "whole_array": {
            "libraries_n": len(partition_audit_rows),
            "distal_midpoint_fragment_range": [
                min(row["distal_midpoint_fragments"] for row in partition_audit_rows),
                max(row["distal_midpoint_fragments"] for row in partition_audit_rows),
            ],
            "terminal_midpoint_share_range_all_libraries": [
                min(row["terminal_whole_array_midpoint_share"] for row in partition_audit_rows),
                max(row["terminal_whole_array_midpoint_share"] for row in partition_audit_rows),
            ],
            "zero_min_whole_array_midpoint_denominators": [
                row["whole_array_midpoint_denominator"] for row in partition_audit_rows if row["time_min"] == 0
            ],
            "zero_min_spacer_counts_separate": [
                row["spacer_midpoint_count_separate_denominator"] for row in partition_audit_rows if row["time_min"] == 0
            ],
        },
        "terminal_rt_local_coverage": {
            "input_mode": coverage_input_mode,
            "input_path": coverage_input_display,
            "source": (
                "saved per-base fragment-span coverage; no alignment or read processing rerun"
                if local_coverage_input is None
                else (
                    "validated bounded local-coverage export; source_coverage_file retains provenance to the saved "
                    "per-base vectors, which were not opened in this mode"
                )
            ),
            "fixed_windows_1based": {
                "terminal_last_100nt": [local_start, BOUNDARY_1BASED - 1],
                "rt_first_100nt": [BOUNDARY_1BASED, local_end],
            },
            "postinfection_libraries_with_depth_on_both_adjacent_boundary_bases": sum(
                row["terminal_end_position_depth"] > 0 and row["rt_start_position_depth"] > 0
                for row in postinfection_continuity
            ),
            "postinfection_libraries_n": len(postinfection_continuity),
            "boundary": (
                "The coverage vector is fragment-span depth. Adjacent coverage is descriptive and does not count "
                "boundary-spanning molecules, establish readthrough, or identify a mature RNA endpoint."
            ),
        },
        "zero_min_metadata": {
            "aliases": [row["sample_alias"] for row in zero_metadata_rows],
            "saved_ena_fields_do_not_state_prephage": True,
            "original_report_design_wording": "0, 5, 15, and 55 min after infection; three cultures per time point",
            "allowed_term": "0-min baseline/background",
            "disallowed_inference": "Do not infer pre-phage solely from the B_0 alias.",
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproduce the bounded P6 science diagnostics from saved tables.",
    )
    parser.add_argument(
        "--local-coverage-input",
        type=Path,
        help=(
            "Use the exported 12-library x 200-position local coverage TSV instead of opening the full "
            "per-library locus_coverage.tsv.gz vectors. Relative paths are resolved from the project root."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    local_coverage_input = args.local_coverage_input
    if local_coverage_input is not None and not local_coverage_input.is_absolute():
        local_coverage_input = ROOT / local_coverage_input

    OUT.mkdir(parents=True, exist_ok=True)
    s4_summary = s4_partition_and_boundary_diagnostics(local_coverage_input)
    summary = {
        "analysis": "Post hoc exploratory identity, temporal and interval diagnostics",
        "date": "2026-09-25",
        "scope": (
            (
                "Saved-table reanalysis plus a fixed 100-nt local extraction from saved per-base coverage; "
                if local_coverage_input is None
                else "Saved-table reanalysis using the validated fixed 100-nt local-coverage export; "
            )
            + "no ViennaRNA rerun, no RNA-seq rerun, no new download, and no original output overwritten."
        ),
        "s1_s3": s1_identity_diagnostics(),
        "s2": s2_temporal_diagnostics(),
        "s4": s4_summary,
    }
    summary_path = OUT / "science_diagnostics_summary.json"
    summary_text = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    summary_path.write_text(summary_text, encoding="utf-8")
    report_json = ROOT / "reports" / "p6_revision_20260925" / "SCIENCE_DIAGNOSTICS.json"
    report_json.parent.mkdir(parents=True, exist_ok=True)
    report_json.write_text(summary_text, encoding="utf-8")
    print(json.dumps({
        "status": "complete",
        "summary": str(summary_path.relative_to(ROOT)).replace("\\", "/"),
        "two_repeat": summary["s1_s3"]["population"],
        "temporal_rows": summary["s2"]["rows_n"],
        "partition_libraries": summary["s4"]["whole_array"]["libraries_n"],
    }, indent=2))


if __name__ == "__main__":
    main()
