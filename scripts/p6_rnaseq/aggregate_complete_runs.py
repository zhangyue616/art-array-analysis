#!/usr/bin/env python3
"""Aggregate complete P6 RNA-seq libraries and render the two frozen result figures."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed" / "p6_rnaseq"
RUNS_DIR = PROCESSED / "runs"
FINAL_DIR = PROCESSED / "aggregate"
PROVISIONAL_DIR = PROCESSED / "provisional"
FIGURE_DIR = ROOT / "figures" / "p6_rnaseq"
LOCAL_PACKAGES = ROOT / "tools" / "p6_rnaseq" / "python_packages"
ENA_META = ROOT / "data" / "source_metadata" / "p6_rnaseq" / "ena_PRJNA836150_read_run.json"

EXPECTED_LIBRARIES = [
    ("SRR19152335", 0, 1),
    ("SRR19152334", 0, 2),
    ("SRR19152331", 0, 3),
    ("SRR19152330", 5, 1),
    ("SRR19152329", 5, 2),
    ("SRR19152328", 5, 3),
    ("SRR19152327", 15, 1),
    ("SRR19152326", 15, 2),
    ("SRR19152325", 15, 3),
    ("SRR19152324", 55, 1),
    ("SRR19152333", 55, 2),
    ("SRR19152332", 55, 3),
]
EXPECTED_BY_RUN = {run: (time_min, culture) for run, time_min, culture in EXPECTED_LIBRARIES}
UNITS = ("U1", "U2", "U3", "U4")
TIMEPOINTS = (0, 5, 15, 55)
POSTINFECTION_TIMEPOINTS = (5, 15, 55)
COMPARISONS = (
    ("5_to_55", 5, 55, "primary"),
    ("5_to_15", 5, 15, "secondary"),
    ("15_to_55", 15, 55, "secondary"),
)
ENDPOINTS = (
    "spacer_midpoint_density_share",
    "spacer_depth_share",
    "core_midpoint_density_share",
    "core_depth_share",
)
PARTITIONS = (
    ("distal_residual", "array_distal_residual"),
    ("U1", "U1_core_to_core"),
    ("U2", "U2_core_to_core"),
    ("U3", "U3_core_to_core"),
    ("U4", "U4_core_to_core"),
    ("terminal_residual", "array_terminal_residual"),
)
CONTEXT_FEATURES = {
    "QPI16999.1": ("replication_annotated", "DNA primase"),
    "QPI17093.1": ("replication_annotated", "DNA-polymerase catalytic subunit"),
    "QPI17172.1": ("replication_annotated", "DnaB-like replicative helicase"),
    "QPI17055.1": ("structural_annotated", "precursor of major head subunit"),
    "QPI17057.1": ("structural_annotated", "terminase large subunit"),
    "QPI17177.1": ("structural_annotated", "tail protein"),
}
RT_FEATURE_ID = "QPI16926.1"
PARTNER_FEATURE_ID = "QPI16927.1"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def serialise(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.12g}"
    if isinstance(value, (list, tuple)):
        return ";".join(serialise(item) for item in value)
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: serialise(row.get(column)) for column in columns})


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def share(values: list[float]) -> list[float | None]:
    denominator = sum(values)
    if denominator <= 0:
        return [None for _ in values]
    return [value / denominator for value in values]


def require_float(value: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite numeric value: {value}")
    return result


def load_complete_libraries() -> tuple[dict[str, dict[str, object]], list[str]]:
    complete: dict[str, dict[str, object]] = {}
    missing: list[str] = []
    with ENA_META.open(encoding="utf-8") as handle:
        ena_by_run = {row["run_accession"]: row for row in json.load(handle)}
    if set(ena_by_run) != set(EXPECTED_BY_RUN):
        raise RuntimeError("Frozen ENA metadata does not match the 12-library mapping")
    required_windows = {
        "array_whole",
        "array_distal_residual",
        "array_terminal_residual",
        "RT_CDS",
        "partner_CDS",
        *(f"{unit}_spacer" for unit in UNITS),
        *(f"{unit}_core_to_core" for unit in UNITS),
    }
    for run, expected_time, expected_culture in EXPECTED_LIBRARIES:
        run_dir = RUNS_DIR / run
        status_path = run_dir / "run_status.json"
        if not status_path.exists():
            missing.append(run)
            continue
        with status_path.open(encoding="utf-8") as handle:
            status = json.load(handle)
        if status.get("status") != "complete":
            missing.append(run)
            continue
        if int(status["time_min"]) != expected_time or int(status["culture"]) != expected_culture:
            raise RuntimeError(f"Frozen time/culture mismatch for {run}")
        if status.get("fastq_validation") != "pass":
            raise RuntimeError(f"FASTQ validation is not pass for {run}")
        if status.get("parser_reached_eof") is not True:
            raise RuntimeError(f"Parser EOF is not recorded for {run}")
        if any(status.get(key) != 0 for key in ("fastp_returncode", "bowtie2_returncode", "findstr_returncode")):
            raise RuntimeError(f"Nonzero pipeline return code in complete status for {run}")
        before_reads = status.get("fastp_summary", {}).get("before_filtering", {}).get("total_reads")
        after_reads = status.get("fastp_summary", {}).get("after_filtering", {}).get("total_reads")
        expected_total_reads = int(ena_by_run[run]["read_count"]) * 2
        if before_reads != expected_total_reads:
            raise RuntimeError(f"fastp before-read count mismatch for {run}: {before_reads} != {expected_total_reads}")
        if not isinstance(after_reads, int) or after_reads <= 0 or after_reads > before_reads:
            raise RuntimeError(f"Invalid fastp after-read count for {run}: {after_reads}")
        windows = read_tsv(run_dir / "window_counts.tsv")
        by_window = {row["window_id"]: row for row in windows}
        if not required_windows.issubset(by_window):
            absent = sorted(required_windows - set(by_window))
            raise RuntimeError(f"Missing frozen windows for {run}: {absent}")
        features = read_tsv(run_dir / "feature_counts.tsv")
        by_feature = {row["feature_id"]: row for row in features}
        absent_features = sorted(set(CONTEXT_FEATURES) - set(by_feature))
        if absent_features:
            raise RuntimeError(f"Missing context features for {run}: {absent_features}")
        whole = by_window["array_whole"]
        partition_rows = [by_window[window_id] for _, window_id in PARTITIONS]
        if sum(int(row["sense_midpoint_fragments"]) for row in partition_rows) != int(
            whole["sense_midpoint_fragments"]
        ):
            raise RuntimeError(f"Whole-array midpoint partition mismatch for {run}")
        if sum(int(row["sense_span_bases"]) for row in partition_rows) != int(whole["sense_span_bases"]):
            raise RuntimeError(f"Whole-array span partition mismatch for {run}")
        complete[run] = {
            "status": status,
            "windows": by_window,
            "features": by_feature,
        }
    return complete, missing


def build_tables(
    complete: dict[str, dict[str, object]],
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    library_rows: list[dict[str, object]] = []
    unit_rows: list[dict[str, object]] = []
    partition_rows: list[dict[str, object]] = []
    context_rows: list[dict[str, object]] = []
    for run, time_min, culture in EXPECTED_LIBRARIES:
        if run not in complete:
            continue
        payload = complete[run]
        status = payload["status"]
        windows = payload["windows"]
        features = payload["features"]
        spacer_windows = [windows[f"{unit}_spacer"] for unit in UNITS]
        core_windows = [windows[f"{unit}_core_to_core"] for unit in UNITS]
        spacer_counts = [int(row["sense_midpoint_fragments"]) for row in spacer_windows]
        spacer_densities = [require_float(row["sense_midpoint_density_per_nt"]) for row in spacer_windows]
        spacer_depths = [require_float(row["sense_mean_fragment_depth"]) for row in spacer_windows]
        core_counts = [int(row["sense_midpoint_fragments"]) for row in core_windows]
        core_densities = [require_float(row["sense_midpoint_density_per_nt"]) for row in core_windows]
        core_depths = [require_float(row["sense_mean_fragment_depth"]) for row in core_windows]
        spacer_density_shares = share(spacer_densities)
        spacer_depth_shares = share(spacer_depths)
        core_density_shares = share(core_densities)
        core_depth_shares = share(core_depths)
        spacer_total = sum(spacer_counts)
        whole = windows["array_whole"]
        distal = windows["array_distal_residual"]
        terminal = windows["array_terminal_residual"]
        whole_midpoints = int(whole["sense_midpoint_fragments"])
        whole_span = int(whole["sense_span_bases"])
        outside_midpoints = int(distal["sense_midpoint_fragments"]) + int(terminal["sense_midpoint_fragments"])
        outside_span = int(distal["sense_span_bases"]) + int(terminal["sense_span_bases"])
        before_reads = status.get("fastp_summary", {}).get("before_filtering", {}).get("total_reads")
        after_reads = status.get("fastp_summary", {}).get("after_filtering", {}).get("total_reads")
        rt_feature = features[RT_FEATURE_ID]
        partner_feature = features[PARTNER_FEATURE_ID]
        library_rows.append(
            {
                "run_accession": run,
                "sample_alias": status["sample_alias"],
                "time_min": time_min,
                "culture": culture,
                "fastp_before_total_reads": before_reads,
                "fastp_after_total_reads": after_reads,
                "bowtie2_overall_alignment_rate_percent": status.get("bowtie2_overall_alignment_rate_percent"),
                "retained_phage_fragments": status["parser_counters"]["retained_phage_fragments"],
                "retained_phage_fragments_plus": status["parser_counters"].get("retained_phage_fragments_plus", 0),
                "retained_phage_fragments_minus": status["parser_counters"].get("retained_phage_fragments_minus", 0),
                "spacer_sense_midpoint_fragments": spacer_total,
                "unit_estimable_ge100": spacer_total >= 100,
                "whole_array_sense_midpoint_fragments": whole_midpoints,
                "whole_array_share_of_retained_phage_fragments": status.get(
                    "array_share_of_retained_phage_fragments"
                ),
                "whole_array_sense_span_bases": whole_span,
                "outside_U1_U4_midpoint_share": outside_midpoints / whole_midpoints if whole_midpoints else None,
                "outside_U1_U4_span_share": outside_span / whole_span if whole_span else None,
                "RT_sense_midpoint_fragments": int(windows["RT_CDS"]["sense_midpoint_fragments"]),
                "partner_sense_midpoint_fragments": int(windows["partner_CDS"]["sense_midpoint_fragments"]),
                "RT_sense_tpm_259_features": require_float(rt_feature["sense_tpm_259_features"]),
                "partner_sense_tpm_259_features": require_float(
                    partner_feature["sense_tpm_259_features"]
                ),
            }
        )
        for index, unit in enumerate(UNITS):
            unit_rows.append(
                {
                    "run_accession": run,
                    "sample_alias": status["sample_alias"],
                    "time_min": time_min,
                    "culture": culture,
                    "unit": unit,
                    "unit_estimable_ge100": spacer_total >= 100,
                    "spacer_total_sense_midpoint_fragments": spacer_total,
                    "spacer_sense_midpoint_fragments": spacer_counts[index],
                    "spacer_midpoint_density_per_nt": spacer_densities[index],
                    "spacer_midpoint_density_share": spacer_density_shares[index],
                    "spacer_mean_fragment_depth": spacer_depths[index],
                    "spacer_depth_share": spacer_depth_shares[index],
                    "core_sense_midpoint_fragments": core_counts[index],
                    "core_midpoint_density_per_nt": core_densities[index],
                    "core_midpoint_density_share": core_density_shares[index],
                    "core_mean_fragment_depth": core_depths[index],
                    "core_depth_share": core_depth_shares[index],
                }
            )
        for partition, window_id in PARTITIONS:
            row = windows[window_id]
            midpoint = int(row["sense_midpoint_fragments"])
            span = int(row["sense_span_bases"])
            partition_rows.append(
                {
                    "run_accession": run,
                    "sample_alias": status["sample_alias"],
                    "time_min": time_min,
                    "culture": culture,
                    "partition": partition,
                    "window_id": window_id,
                    "sense_midpoint_fragments": midpoint,
                    "whole_array_midpoint_share": midpoint / whole_midpoints if whole_midpoints else None,
                    "sense_span_bases": span,
                    "whole_array_span_share": span / whole_span if whole_span else None,
                }
            )
        for feature_id, (group, expected_product) in CONTEXT_FEATURES.items():
            row = features[feature_id]
            if row["product"] != expected_product:
                raise RuntimeError(
                    f"Frozen context annotation mismatch for {feature_id}: {row['product']} != {expected_product}"
                )
            context_rows.append(
                {
                    "run_accession": run,
                    "sample_alias": status["sample_alias"],
                    "time_min": time_min,
                    "culture": culture,
                    "context_group": group,
                    "feature_id": feature_id,
                    "product": row["product"],
                    "sense_midpoint_fragments": int(row["sense_midpoint_fragments"]),
                    "antisense_midpoint_fragments": int(row["antisense_midpoint_fragments"]),
                    "sense_tpm_259_features": require_float(row["sense_tpm_259_features"]),
                }
            )
    return library_rows, unit_rows, partition_rows, context_rows


def build_contrasts(
    unit_rows: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    index = {
        (int(row["time_min"]), int(row["culture"]), str(row["unit"])): row
        for row in unit_rows
    }
    paired_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for comparison, start_time, end_time, role in COMPARISONS:
        for endpoint in ENDPOINTS:
            for unit in UNITS:
                available_changes: list[float] = []
                all_estimable = True
                for culture in (1, 2, 3):
                    start = index.get((start_time, culture, unit))
                    end = index.get((end_time, culture, unit))
                    if start is None or end is None:
                        paired_rows.append(
                            {
                                "comparison": comparison,
                                "comparison_role": role,
                                "endpoint": endpoint,
                                "unit": unit,
                                "culture": culture,
                                "start_time_min": start_time,
                                "end_time_min": end_time,
                                "start_value": None,
                                "end_value": None,
                                "delta": None,
                                "both_libraries_estimable": None,
                                "pair_status": "missing_incomplete_run",
                            }
                        )
                        all_estimable = False
                        continue
                    start_value = start[endpoint]
                    end_value = end[endpoint]
                    estimable = bool(start["unit_estimable_ge100"]) and bool(end["unit_estimable_ge100"])
                    if start_value is None or end_value is None:
                        delta = None
                        pair_status = "undefined_zero_denominator"
                    else:
                        delta = float(end_value) - float(start_value)
                        available_changes.append(delta)
                        pair_status = "complete"
                    all_estimable = all_estimable and estimable
                    paired_rows.append(
                        {
                            "comparison": comparison,
                            "comparison_role": role,
                            "endpoint": endpoint,
                            "unit": unit,
                            "culture": culture,
                            "start_time_min": start_time,
                            "end_time_min": end_time,
                            "start_value": start_value,
                            "end_value": end_value,
                            "delta": delta,
                            "both_libraries_estimable": estimable,
                            "pair_status": pair_status,
                        }
                    )
                evaluable = len(available_changes) == 3
                same_positive = evaluable and all(delta > 0 for delta in available_changes)
                same_negative = evaluable and all(delta < 0 for delta in available_changes)
                same_direction = same_positive or same_negative
                direction = "increase" if same_positive else "decrease" if same_negative else "mixed_or_zero"
                median_absolute_delta = (
                    statistics.median(abs(delta) for delta in available_changes) if evaluable else None
                )
                meets_numeric_rule = bool(
                    evaluable
                    and all_estimable
                    and same_direction
                    and median_absolute_delta is not None
                    and (
                        median_absolute_delta > 0.05
                        or math.isclose(median_absolute_delta, 0.05, abs_tol=1e-12)
                    )
                )
                summary_rows.append(
                    {
                        "comparison": comparison,
                        "comparison_role": role,
                        "endpoint": endpoint,
                        "endpoint_role": (
                            "primary"
                            if endpoint == "spacer_midpoint_density_share"
                            else "confirmatory"
                            if endpoint == "spacer_depth_share"
                            else "repeat_sensitive_sensitivity"
                        ),
                        "unit": unit,
                        "complete_culture_pairs": len(available_changes),
                        "all_six_libraries_estimable": all_estimable if evaluable else False,
                        "all_three_changes_same_direction": same_direction,
                        "direction": direction,
                        "median_absolute_share_change": median_absolute_delta,
                        "meets_predeclared_numeric_rule": meets_numeric_rule,
                        "evaluation_status": "complete" if evaluable else "incomplete",
                    }
                )
    sensitivity_rows: list[dict[str, object]] = []
    by_summary = {
        (row["comparison"], row["unit"], row["endpoint"]): row for row in summary_rows
    }
    for comparison, _, _, role in COMPARISONS:
        for unit in UNITS:
            row: dict[str, object] = {
                "comparison": comparison,
                "comparison_role": role,
                "unit": unit,
            }
            for endpoint in ENDPOINTS:
                source = by_summary[(comparison, unit, endpoint)]
                prefix = endpoint.replace("_share", "")
                row[f"{prefix}_median_absolute_delta"] = source["median_absolute_share_change"]
                row[f"{prefix}_direction"] = source["direction"]
                row[f"{prefix}_same_direction"] = source["all_three_changes_same_direction"]
                row[f"{prefix}_meets_numeric_rule"] = source["meets_predeclared_numeric_rule"]
            sensitivity_rows.append(row)
    return paired_rows, summary_rows, sensitivity_rows


def build_timepoint_summary(unit_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, int], list[float]] = defaultdict(list)
    for row in unit_rows:
        for endpoint in ENDPOINTS:
            value = row[endpoint]
            if value is not None:
                grouped[(str(row["unit"]), endpoint, int(row["time_min"]))].append(float(value))
    rows: list[dict[str, object]] = []
    for unit in UNITS:
        for endpoint in ENDPOINTS:
            for time_min in TIMEPOINTS:
                values = grouped.get((unit, endpoint, time_min), [])
                rows.append(
                    {
                        "unit": unit,
                        "endpoint": endpoint,
                        "time_min": time_min,
                        "available_cultures": len(values),
                        "culture_values": values,
                        "culture_mean": (sum(values) / len(values)) if values else None,
                        "culture_median": statistics.median(values) if values else None,
                    }
                )
    return rows


def build_context_group_timepoint_summary(
    context_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Summarize the six preselected annotations without treating genes as replicates."""
    rows: list[dict[str, object]] = []
    for time_min in TIMEPOINTS:
        for context_group in ("replication_annotated", "structural_annotated"):
            culture_medians: dict[int, float] = {}
            for culture in (1, 2, 3):
                values = [
                    float(row["sense_tpm_259_features"])
                    for row in context_rows
                    if int(row["time_min"]) == time_min
                    and int(row["culture"]) == culture
                    and row["context_group"] == context_group
                ]
                if len(values) != 3:
                    raise RuntimeError(
                        f"Expected three selected {context_group} features at "
                        f"{time_min} min culture {culture}; observed {len(values)}"
                    )
                culture_medians[culture] = statistics.median(values)
            rows.append(
                {
                    "time_min": time_min,
                    "context_group": context_group,
                    "selected_feature_count": 3,
                    "culture_1_median_sense_tpm_259_features": culture_medians[1],
                    "culture_2_median_sense_tpm_259_features": culture_medians[2],
                    "culture_3_median_sense_tpm_259_features": culture_medians[3],
                    "median_across_culture_level_feature_medians": statistics.median(
                        culture_medians.values()
                    ),
                    "scope": "descriptive annotation context only",
                }
            )
    return rows


def build_rt_partner_timepoint_summary(
    library_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for time_min in TIMEPOINTS:
        time_rows = sorted(
            (row for row in library_rows if int(row["time_min"]) == time_min),
            key=lambda row: int(row["culture"]),
        )
        if len(time_rows) != 3:
            raise RuntimeError(f"Expected three libraries at {time_min} min")
        for target, field in (
            ("RT", "RT_sense_tpm_259_features"),
            ("direct_partner", "partner_sense_tpm_259_features"),
        ):
            values = [float(row[field]) for row in time_rows]
            rows.append(
                {
                    "time_min": time_min,
                    "target": target,
                    "culture_1_sense_tpm_259_features": values[0],
                    "culture_2_sense_tpm_259_features": values[1],
                    "culture_3_sense_tpm_259_features": values[2],
                    "median_across_cultures": statistics.median(values),
                    "scope": "descriptive context only",
                }
            )
    return rows


def overall_conclusion(
    complete_all: bool,
    unit_rows: list[dict[str, object]],
    contrast_rows: list[dict[str, object]],
    timepoint_rows: list[dict[str, object]],
) -> dict[str, object]:
    if not complete_all:
        return {
            "status": "not_evaluated_incomplete",
            "reason": "All 12 complete libraries are required before a temporal conclusion.",
        }
    primary = [
        row
        for row in contrast_rows
        if row["comparison"] == "5_to_55" and row["endpoint"] == "spacer_midpoint_density_share"
    ]
    supported_units = [str(row["unit"]) for row in primary if row["meets_predeclared_numeric_rule"]]
    secondary_supported = {
        comparison: [
            str(row["unit"])
            for row in contrast_rows
            if row["comparison"] == comparison
            and row["endpoint"] == "spacer_midpoint_density_share"
            and row["meets_predeclared_numeric_rule"]
        ]
        for comparison in ("5_to_15", "15_to_55")
    }
    any_secondary_supported = any(secondary_supported.values())
    postinfection_estimability = {
        (str(row["run_accession"]), int(row["time_min"]), int(row["culture"])): bool(
            row["unit_estimable_ge100"]
        )
        for row in unit_rows
        if int(row["time_min"]) in POSTINFECTION_TIMEPOINTS
    }
    all_postinfection_libraries_estimable = (
        len(postinfection_estimability) == 9 and all(postinfection_estimability.values())
    )
    ranges: dict[str, float] = {}
    for unit in UNITS:
        values = [
            float(row["culture_mean"])
            for row in timepoint_rows
            if row["unit"] == unit
            and row["endpoint"] == "spacer_midpoint_density_share"
            and int(row["time_min"]) in POSTINFECTION_TIMEPOINTS
            and row["culture_mean"] is not None
        ]
        if len(values) != 3:
            raise RuntimeError(f"Missing complete post-infection means for {unit}")
        ranges[unit] = max(values) - min(values)
    if supported_units:
        label = "supported_material_temporal_change"
    elif any_secondary_supported:
        label = "indeterminate_secondary_material_change"
    elif not all_postinfection_libraries_estimable:
        label = "indeterminate_low_estimability"
    elif all(value < 0.05 for value in ranges.values()):
        label = "stable_by_predeclared_rule"
    else:
        label = "indeterminate_by_predeclared_rule"
    return {
        "status": "evaluated",
        "primary_comparison": "5_to_55",
        "supported_material_change_units": supported_units,
        "secondary_supported_material_change_units": secondary_supported,
        "secondary_findings_do_not_replace_primary_5_to_55": True,
        "all_nine_postinfection_libraries_estimable": all_postinfection_libraries_estimable,
        "postinfection_timepoint_mean_ranges": ranges,
        "conclusion": label,
        "rule": "three culture-paired changes in the same nonzero direction, median absolute share change >= 0.05, and all six libraries >=100 spacer sense midpoint fragments; stable additionally requires no unit to meet this rule in any of the three predeclared comparisons and every unit mean range across 5/15/55 min <0.05",
    }


def render_figures(
    unit_rows: list[dict[str, object]],
    partition_rows: list[dict[str, object]],
    paired_rows: list[dict[str, object]],
    contrast_rows: list[dict[str, object]],
) -> list[str]:
    sys.path.insert(0, str(LOCAL_PACKAGES))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    unit_colours = {"U1": "#0072B2", "U2": "#009E73", "U3": "#CC79A7", "U4": "#D55E00"}
    culture_markers = {1: "o", 2: "s", 3: "^"}
    culture_styles = {1: "-", 2: "--", 3: ":"}
    rc = {
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.labelsize": 8.5,
        "axes.titlesize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.linewidth": 0.7,
        "lines.linewidth": 1.2,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    }
    index = {
        (int(row["time_min"]), int(row["culture"]), str(row["unit"])): row
        for row in unit_rows
    }
    primary_summary = {
        str(row["unit"]): row
        for row in contrast_rows
        if row["comparison"] == "5_to_55" and row["endpoint"] == "spacer_midpoint_density_share"
    }
    outputs: list[str] = []
    with plt.rc_context(rc):
        fig, axes = plt.subplots(2, 2, figsize=(170 / 25.4, 125 / 25.4), sharex=True, sharey=True, layout="constrained")
        for panel_index, (axis, unit) in enumerate(zip(axes.flat, UNITS)):
            for culture in (1, 2, 3):
                post_rows = [index[(time_min, culture, unit)] for time_min in POSTINFECTION_TIMEPOINTS]
                axis.plot(
                    POSTINFECTION_TIMEPOINTS,
                    [float(row["spacer_midpoint_density_share"]) for row in post_rows],
                    color=unit_colours[unit],
                    marker=culture_markers[culture],
                    linestyle=culture_styles[culture],
                    markersize=4,
                    markeredgewidth=0.7,
                )
                background = index[(0, culture, unit)]
                if background["unit_estimable_ge100"] and background["spacer_midpoint_density_share"] is not None:
                    axis.scatter(
                        [0],
                        [float(background["spacer_midpoint_density_share"])],
                        marker=culture_markers[culture],
                        facecolors="white",
                        edgecolors=unit_colours[unit],
                        linewidths=0.9,
                        s=20,
                        zorder=3,
                    )
            summary = primary_summary[unit]
            median_delta = float(summary["median_absolute_share_change"])
            direction = str(summary["direction"])
            direction_text = {"increase": "same ↑", "decrease": "same ↓"}.get(direction, "mixed/zero")
            axis.text(
                0.02,
                0.72 if unit == "U3" else 0.96,
                f"median |Δ5→55|={median_delta:.3f}; {direction_text}",
                transform=axis.transAxes,
                ha="left",
                va="top",
                fontsize=8,
            )
            axis.set_title(f"{unit} repeat-excluded spacer", color=unit_colours[unit], fontweight="bold")
            axis.set_xlim(-2, 58)
            axis.set_ylim(0, 1)
            axis.set_xticks(TIMEPOINTS)
            axis.grid(axis="y", color="#D9D9D9", linewidth=0.5)
            axis.spines[["top", "right"]].set_visible(False)
            axis.text(-0.12, 1.05, chr(ord("A") + panel_index), transform=axis.transAxes, fontweight="bold", fontsize=10)
        axes[1, 0].set_xlabel("Minutes after infection")
        axes[1, 1].set_xlabel("Minutes after infection")
        axes[0, 0].set_ylabel("Within-chain density share")
        axes[1, 0].set_ylabel("Within-chain density share")
        handles = [
            Line2D([0], [0], color="#4D4D4D", marker=culture_markers[culture], linestyle=culture_styles[culture], label=f"Culture {culture}")
            for culture in (1, 2, 3)
        ]
        handles.append(
            Line2D(
                [0],
                [0],
                color="#4D4D4D",
                marker="o",
                markerfacecolor="white",
                linestyle="None",
                label="0 min (none estimable)",
            )
        )
        fig.legend(handles=handles, loc="outside lower center", ncol=4, frameon=False)
        stem = FIGURE_DIR / "p6_rnaseq_unit_density_shares"
        for suffix in ("pdf", "svg", "png"):
            path = stem.with_suffix(f".{suffix}")
            fig.savefig(path, dpi=300 if suffix == "png" else None)
            outputs.append(str(path.relative_to(ROOT)))
        plt.close(fig)

        fig = plt.figure(figsize=(170 / 25.4, 122 / 25.4), layout="constrained")
        grid = fig.add_gridspec(2, 2, height_ratios=(1, 1.05))
        ax_confirm = fig.add_subplot(grid[0, 0])
        ax_sensitivity = fig.add_subplot(grid[0, 1])
        ax_partition = fig.add_subplot(grid[1, :])
        paired_index = {
            (str(row["endpoint"]), str(row["unit"]), int(row["culture"])): row
            for row in paired_rows
            if row["comparison"] == "5_to_55"
        }
        confirm_points: list[tuple[str, int, float, float]] = []
        sensitivity_points: list[tuple[str, int, float, float]] = []
        for unit in UNITS:
            for culture in (1, 2, 3):
                primary = paired_index[("spacer_midpoint_density_share", unit, culture)]
                confirm = paired_index[("spacer_depth_share", unit, culture)]
                core = paired_index[("core_midpoint_density_share", unit, culture)]
                if (
                    primary["pair_status"] == "complete"
                    and confirm["pair_status"] == "complete"
                    and primary["both_libraries_estimable"]
                    and confirm["both_libraries_estimable"]
                ):
                    confirm_points.append((unit, culture, float(primary["delta"]), float(confirm["delta"])))
                if (
                    primary["pair_status"] == "complete"
                    and core["pair_status"] == "complete"
                    and primary["both_libraries_estimable"]
                    and core["both_libraries_estimable"]
                ):
                    sensitivity_points.append((unit, culture, float(primary["delta"]), float(core["delta"])))
        for unit, culture, x_value, y_value in confirm_points:
            ax_confirm.scatter(
                [x_value],
                [y_value],
                color=unit_colours[unit],
                marker=culture_markers[culture],
                s=23,
                linewidths=0.45,
                edgecolors="white",
            )
        for unit, culture, x_value, y_value in sensitivity_points:
            ax_sensitivity.scatter(
                [x_value],
                [y_value],
                color=unit_colours[unit],
                marker=culture_markers[culture],
                s=23,
                linewidths=0.45,
                edgecolors="white",
            )
        all_deltas = [abs(value) for _, _, x_value, y_value in confirm_points + sensitivity_points for value in (x_value, y_value)]
        limit = max(0.05, max(all_deltas, default=0.05) * 1.12)
        for axis, point_count in ((ax_confirm, len(confirm_points)), (ax_sensitivity, len(sensitivity_points))):
            axis.plot([-limit, limit], [-limit, limit], color="#777777", linewidth=0.8, linestyle="--", zorder=0)
            axis.axhline(0, color="#BDBDBD", linewidth=0.6, zorder=0)
            axis.axvline(0, color="#BDBDBD", linewidth=0.6, zorder=0)
            axis.set_xlim(-limit, limit)
            axis.set_ylim(-limit, limit)
            axis.set_aspect("equal", adjustable="box")
            axis.grid(color="#E0E0E0", linewidth=0.45)
            axis.spines[["top", "right"]].set_visible(False)
            axis.text(
                0.02,
                0.98,
                f"evaluable pairs: {point_count}/12",
                transform=axis.transAxes,
                ha="left",
                va="top",
                fontsize=8,
            )
        ax_confirm.set_xlabel("Δ spacer density share")
        ax_confirm.set_ylabel("Δ spacer depth share")
        ax_sensitivity.set_xlabel("Δ spacer density share")
        ax_sensitivity.set_ylabel("Δ core-to-core density share")
        ax_confirm.text(-0.16, 1.06, "A", transform=ax_confirm.transAxes, fontweight="bold", fontsize=10)
        ax_sensitivity.text(-0.16, 1.06, "B", transform=ax_sensitivity.transAxes, fontweight="bold", fontsize=10)
        unit_handles = [
            Line2D([0], [0], marker="o", color="none", markerfacecolor=unit_colours[unit], markeredgecolor="white", label=unit)
            for unit in UNITS
        ]
        culture_handles = [
            Line2D([0], [0], marker=culture_markers[culture], color="#4D4D4D", linestyle="None", label=f"Culture {culture}")
            for culture in (1, 2, 3)
        ]
        ax_confirm.legend(handles=unit_handles, loc="lower right", frameon=False, ncol=2)
        ax_sensitivity.legend(handles=culture_handles, loc="lower right", frameon=False)

        partition_index = {
            (str(row["run_accession"]), str(row["partition"])): row for row in partition_rows
        }
        ordered_runs = [run for run, _, _ in EXPECTED_LIBRARIES]
        x = list(range(len(ordered_runs)))
        bottoms = [0.0] * len(ordered_runs)
        partition_colours = {
            "distal_residual": "#BDBDBD",
            "U1": unit_colours["U1"],
            "U2": unit_colours["U2"],
            "U3": unit_colours["U3"],
            "U4": unit_colours["U4"],
            "terminal_residual": "#666666",
        }
        partition_hatches = {"distal_residual": "//", "terminal_residual": "\\\\"}
        for partition, _ in PARTITIONS:
            values = [float(partition_index[(run, partition)]["whole_array_midpoint_share"]) for run in ordered_runs]
            ax_partition.bar(
                x,
                values,
                bottom=bottoms,
                width=0.72,
                color=partition_colours[partition],
                edgecolor="white",
                linewidth=0.35,
                hatch=partition_hatches.get(partition),
                label=partition.replace("_", " "),
            )
            bottoms = [bottom + value for bottom, value in zip(bottoms, values)]
        labels = [f"{time_min}-{culture}" for _, time_min, culture in EXPECTED_LIBRARIES]
        ax_partition.set_xticks(x, labels, rotation=45, ha="right")
        ax_partition.set_ylim(0, 1)
        ax_partition.set_ylabel("Whole-array midpoint share")
        ax_partition.set_xlabel("Time (min) – culture")
        ax_partition.axvline(2.5, color="#4D4D4D", linewidth=0.8, linestyle="--")
        ax_partition.text(1, 1.025, "0 min background", ha="center", va="bottom", fontsize=8)
        ax_partition.text(7, 1.025, "post-infection", ha="center", va="bottom", fontsize=8)
        ax_partition.spines[["top", "right"]].set_visible(False)
        ax_partition.grid(axis="y", color="#E0E0E0", linewidth=0.45)
        partition_handles, partition_labels = ax_partition.get_legend_handles_labels()
        fig.legend(
            handles=partition_handles,
            labels=partition_labels,
            loc="outside lower center",
            ncol=6,
            frameon=False,
        )
        ax_partition.text(-0.075, 1.07, "C", transform=ax_partition.transAxes, fontweight="bold", fontsize=10)
        stem = FIGURE_DIR / "p6_rnaseq_confirmation_sensitivity_partition"
        for suffix in ("pdf", "svg", "png"):
            path = stem.with_suffix(f".{suffix}")
            fig.savefig(path, dpi=300 if suffix == "png" else None)
            outputs.append(str(path.relative_to(ROOT)))
        plt.close(fig)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="Write explicitly provisional machine tables from complete runs; never render final figures.",
    )
    args = parser.parse_args()
    complete, missing = load_complete_libraries()
    complete_all = not missing and len(complete) == len(EXPECTED_LIBRARIES)
    if not complete_all and not args.allow_incomplete:
        raise SystemExit(f"Incomplete run set ({len(complete)}/12); missing or non-complete: {', '.join(missing)}")
    output_dir = FINAL_DIR if complete_all else PROVISIONAL_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    library_rows, unit_rows, partition_rows, context_rows = build_tables(complete)
    paired_rows, contrast_rows, sensitivity_rows = build_contrasts(unit_rows)
    timepoint_rows = build_timepoint_summary(unit_rows)
    context_group_timepoint_rows = build_context_group_timepoint_summary(context_rows)
    rt_partner_timepoint_rows = build_rt_partner_timepoint_summary(library_rows)
    conclusion = overall_conclusion(complete_all, unit_rows, contrast_rows, timepoint_rows)

    write_tsv(
        output_dir / "library_metrics.tsv",
        library_rows,
        [
            "run_accession",
            "sample_alias",
            "time_min",
            "culture",
            "fastp_before_total_reads",
            "fastp_after_total_reads",
            "bowtie2_overall_alignment_rate_percent",
            "retained_phage_fragments",
            "retained_phage_fragments_plus",
            "retained_phage_fragments_minus",
            "spacer_sense_midpoint_fragments",
            "unit_estimable_ge100",
            "whole_array_sense_midpoint_fragments",
            "whole_array_share_of_retained_phage_fragments",
            "whole_array_sense_span_bases",
            "outside_U1_U4_midpoint_share",
            "outside_U1_U4_span_share",
            "RT_sense_midpoint_fragments",
            "partner_sense_midpoint_fragments",
            "RT_sense_tpm_259_features",
            "partner_sense_tpm_259_features",
        ],
    )
    write_tsv(output_dir / "unit_metrics.tsv", unit_rows, list(unit_rows[0]) if unit_rows else [])
    write_tsv(
        output_dir / "whole_array_partitions.tsv",
        partition_rows,
        list(partition_rows[0]) if partition_rows else [],
    )
    write_tsv(
        output_dir / "context_feature_metrics.tsv",
        context_rows,
        list(context_rows[0]) if context_rows else [],
    )
    write_tsv(output_dir / "paired_contrasts.tsv", paired_rows, list(paired_rows[0]) if paired_rows else [])
    write_tsv(output_dir / "contrast_summary.tsv", contrast_rows, list(contrast_rows[0]) if contrast_rows else [])
    write_tsv(
        output_dir / "sensitivity_summary.tsv",
        sensitivity_rows,
        list(sensitivity_rows[0]) if sensitivity_rows else [],
    )
    write_tsv(
        output_dir / "timepoint_summary.tsv",
        timepoint_rows,
        list(timepoint_rows[0]) if timepoint_rows else [],
    )
    write_tsv(
        output_dir / "context_group_timepoint_summary.tsv",
        context_group_timepoint_rows,
        list(context_group_timepoint_rows[0]) if context_group_timepoint_rows else [],
    )
    write_tsv(
        output_dir / "rt_partner_timepoint_summary.tsv",
        rt_partner_timepoint_rows,
        list(rt_partner_timepoint_rows[0]) if rt_partner_timepoint_rows else [],
    )
    figure_outputs = render_figures(unit_rows, partition_rows, paired_rows, contrast_rows) if complete_all else []
    summary = {
        "status": "complete" if complete_all else "incomplete_provisional",
        "completed_run_count": len(complete),
        "expected_run_count": len(EXPECTED_LIBRARIES),
        "completed_runs": [run for run, _, _ in EXPECTED_LIBRARIES if run in complete],
        "missing_or_noncomplete_runs": missing,
        "final_temporal_conclusion_allowed": complete_all,
        "zero_minute_handling": "descriptive background, retained separately and not subtracted",
        "independent_unit": "culture (n=3 paired across time)",
        "primary_endpoint": "repeat-excluded spacer midpoint density share",
        "confirmatory_endpoint": "repeat-excluded spacer fragment-depth share",
        "sensitivity_endpoint": "core-to-core midpoint-density and depth shares",
        "primary_comparison": "5_to_55",
        "secondary_comparisons": ["5_to_15", "15_to_55"],
        "predeclared_material_change_rule": {
            "all_three_culture_pairs_same_nonzero_direction": True,
            "median_absolute_share_change_at_least": 0.05,
            "minimum_spacer_sense_midpoint_fragments_per_involved_library": 100,
        },
        "fastp_field_boundary": {
            "accepted": ["before_filtering.total_reads", "after_filtering.total_reads"],
            "excluded": "base-dependent, mean-length, GC, and related quality-base summary fields from fastp 1.3.3 Windows",
            "used_by_fragment_counts_density_tpm_or_conclusion": False,
            "raw_reports_modified": False,
        },
        "context_feature_selection": {
            "rule": "Three exact GenBank product annotations for replication machinery (primase, DNA polymerase catalytic subunit, replicative helicase) and three for virion assembly/packaging (major head precursor, terminase large subunit, tail protein); descriptive context only.",
            "group_timepoint_summary_definition": "Within each culture and annotation group, take the median TPM-like value across the three selected features; then report the median of those three culture-level medians. Genes are not treated as biological replicates.",
            "group_timepoint_summary_file": "context_group_timepoint_summary.tsv",
            "features": {
                feature_id: {"group": group, "product": product}
                for feature_id, (group, product) in CONTEXT_FEATURES.items()
            },
        },
        "rt_partner_context": {
            "feature_ids": {"RT": RT_FEATURE_ID, "direct_partner": PARTNER_FEATURE_ID},
            "timepoint_summary_file": "rt_partner_timepoint_summary.tsv",
            "scope": "descriptive context only; no additional temporal decision rule",
        },
        "conclusion": conclusion,
        "figure_outputs": figure_outputs,
        "provisional_outputs_are_not_final_figures": not complete_all,
    }
    write_json(output_dir / "analysis_summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
