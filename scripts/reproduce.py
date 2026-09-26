#!/usr/bin/env python3
"""Recompute key saved-table statistics and render the accepted article figures.

The default route is deliberately lightweight and offline.  It reads the
published processed tables, recomputes compact audit statistics, and renders
the ten accepted figures.  It does not download reads, rerun RNA folding,
repeat the shuffled calibration, or run the full RNA-seq alignment workflow.
"""

from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import math
import platform
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


FIGURE_INPUTS = (
    "data/processed/p2_loci/fixed_10_loci.tsv",
    "data/processed/p2_loci/coding_context.tsv",
    "data/processed/p2_comparison/unified_array_annotations.tsv",
    "data/processed/p2_comparison/unified_repeat_copies.tsv",
    "data/processed/p2_comparison/unified_units.tsv",
    "data/processed/p2_comparison/boundary_alternatives.tsv",
    "data/processed/p2_comparison/locus_pairwise_identity.tsv",
    "data/processed/p2_comparison/locus_pair_unit_correspondence_summary.tsv",
    "data/processed/p2_comparison/supported_unit_correspondences.tsv",
    "data/processed/p2_comparison/edge_scan_audit.tsv",
    "data/processed/p2_comparison/cross_group_medoid_copy_mapping.tsv",
    "data/processed/p2_comparison/cross_group_medoid_array_alignment_blocks.tsv",
    "data/processed/p3_calibration/null_replicates.tsv",
    "data/processed/p3_calibration/p3_calibration_summary.json",
    "data/processed/p6_revision_20260925/s1_two_repeat_pair_audit.tsv",
    "data/processed/p6_revision_20260925/s1_locus_pair_identity_range_diagnostic.tsv",
    "data/processed/p6_rna/dinucleotide_null_replicates.tsv",
    "data/processed/p6_rna/p6_rna_summary.json",
    "data/processed/p6_rna/representative_alignment_bpp.tsv",
    "data/processed/p6_rna/representative_alignment_metadata.tsv",
    "data/processed/p6_rnaseq/aggregate/unit_metrics.tsv",
    "data/processed/p6_rnaseq/aggregate/paired_contrasts.tsv",
    "data/processed/p6_rnaseq/aggregate/whole_array_partitions.tsv",
    "data/processed/p6_protein/protein_feature_summary.tsv",
    "data/processed/p6_protein/protein_pairwise_identity.tsv",
)

STATISTIC_INPUTS = (
    "data/processed/p2_loci/fixed_10_loci.tsv",
    "data/processed/p2_comparison/supported_unit_correspondences.tsv",
    "data/processed/p3_calibration/p3_calibration_summary.json",
    "data/processed/p6_rna/unit_pair_structure_similarity.tsv",
    "data/processed/p6_rna/p6_rna_summary.json",
    "data/processed/p6_rnaseq/aggregate/analysis_summary.json",
    "data/processed/p6_rnaseq/aggregate/paired_contrasts.tsv",
    "data/processed/p6_rnaseq/aggregate/whole_array_partitions.tsv",
    "data/processed/p6_protein/protein_feature_summary.tsv",
    "data/processed/p6_protein/protein_pairwise_identity.tsv",
)

IDENTITY_FIELDS = (
    "context_global_identity",
    "saved_core_to_core_global_identity",
    "saved_trimmed_spacer_global_identity",
)


def parse_args() -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=default_root)
    parser.add_argument("--output", type=Path, default=Path("results/reproduced"))
    parser.add_argument(
        "--figures",
        choices=("all", "key", "none"),
        default="all",
        help="Render all ten figures, Figures 5 and 7 only, or no figures.",
    )
    parser.add_argument(
        "--check-inputs",
        action="store_true",
        help="Validate required files and report the plan without writing outputs.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def rank_average(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        average = (start + 1 + end) / 2.0
        for position in range(start, end):
            ranks[order[position]] = average
        start = end
    return ranks


def pearson(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or len(left) < 2:
        raise ValueError("Correlation requires paired vectors of length at least two")
    left_mean = statistics.fmean(left)
    right_mean = statistics.fmean(right)
    numerator = sum((x - left_mean) * (y - right_mean) for x, y in zip(left, right))
    left_ss = sum((x - left_mean) ** 2 for x in left)
    right_ss = sum((y - right_mean) ** 2 for y in right)
    if left_ss == 0.0 or right_ss == 0.0:
        return math.nan
    return numerator / math.sqrt(left_ss * right_ss)


def spearman(left: list[float], right: list[float]) -> float:
    return pearson(rank_average(left), rank_average(right))


def validate_inputs(root: Path, figures: str) -> list[str]:
    required = set(STATISTIC_INPUTS)
    if figures != "none":
        required.update(FIGURE_INPUTS)
    missing = [relative for relative in sorted(required) if not (root / relative).is_file()]
    if missing:
        raise SystemExit("Missing required inputs:\n  " + "\n  ".join(missing))
    return sorted(required)


def structure_statistics(root: Path) -> tuple[dict[str, object], list[dict[str, object]]]:
    rows = [
        row
        for row in read_tsv(root / "data/processed/p6_rna/unit_pair_structure_similarity.tsv")
        if row["context"] == "two_repeat"
    ]
    if len(rows) != 252:
        raise ValueError(f"Expected 252 two-repeat rows, found {len(rows)}")
    supported = [row for row in rows if row["supported_correspondence"].lower() == "true"]
    other = [row for row in rows if row["supported_correspondence"].lower() == "false"]
    if (len(supported), len(other)) != (35, 217):
        raise ValueError(f"Expected 35 supported and 217 other rows, found {len(supported)} and {len(other)}")

    correlation_rows: list[dict[str, object]] = []
    for identity_field in IDENTITY_FIELDS:
        for population, subset in (("all_252", rows), ("supported_35", supported)):
            identities = [float(row[identity_field]) for row in subset]
            overlaps = [float(row["spacer_spacer_bpp_overlap"]) for row in subset]
            correlation_rows.append(
                {
                    "context": "two_repeat",
                    "population": population,
                    "n": len(subset),
                    "identity_field": identity_field,
                    "bpp_field": "spacer_spacer_bpp_overlap",
                    "spearman_rho": f"{spearman(identities, overlaps):.12g}",
                }
            )

    summary = {
        "context": "two_repeat",
        "all_pairs": len(rows),
        "supported_pairs": len(supported),
        "other_pairs": len(other),
        "locus_pairs": len({row["locus_pair_id"] for row in rows}),
        "supported_locus_pairs": len({row["locus_pair_id"] for row in supported}),
        "supported_spacer_bpp_overlap_median": statistics.median(
            float(row["spacer_spacer_bpp_overlap"]) for row in supported
        ),
        "other_spacer_bpp_overlap_median": statistics.median(
            float(row["spacer_spacer_bpp_overlap"]) for row in other
        ),
    }
    return summary, correlation_rows


def temporal_statistics(root: Path) -> tuple[dict[str, object], list[dict[str, object]]]:
    analysis = json.loads(
        (root / "data/processed/p6_rnaseq/aggregate/analysis_summary.json").read_text(encoding="utf-8")
    )
    threshold = float(
        analysis["predeclared_material_change_rule"]["median_absolute_share_change_at_least"]
    )
    primary = analysis["primary_comparison"]
    endpoint_description = analysis["primary_endpoint"]
    endpoint = "spacer_midpoint_density_share"
    rows = [
        row
        for row in read_tsv(root / "data/processed/p6_rnaseq/aggregate/paired_contrasts.tsv")
        if row["comparison"] == primary
        and row["endpoint"] == endpoint
        and row["both_libraries_estimable"].lower() == "true"
        and row["pair_status"] == "complete"
    ]
    by_unit: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_unit[row["unit"]].append(row)
    if sorted(by_unit) != ["U1", "U2", "U3", "U4"] or any(len(group) != 3 for group in by_unit.values()):
        raise ValueError("Expected three complete paired 5-to-55 contrasts for each of U1-U4")

    output_rows: list[dict[str, object]] = []
    for unit in sorted(by_unit):
        group = sorted(by_unit[unit], key=lambda row: int(row["culture"]))
        deltas = [float(row["delta"]) for row in group]
        starts = [float(row["start_value"]) for row in group]
        relative = [100.0 * delta / start if start else math.nan for delta, start in zip(deltas, starts)]
        directions = {0 if delta == 0 else (1 if delta > 0 else -1) for delta in deltas}
        same_nonzero_direction = len(directions) == 1 and 0 not in directions
        median_abs = statistics.median(abs(delta) for delta in deltas)
        output_rows.append(
            {
                "comparison": primary,
                "endpoint": endpoint,
                "unit": unit,
                "cultures_n": len(group),
                "median_delta_share": f"{statistics.median(deltas):.12g}",
                "median_absolute_delta_share": f"{median_abs:.12g}",
                "minimum_relative_change_percent": f"{min(relative):.12g}",
                "maximum_relative_change_percent": f"{max(relative):.12g}",
                "all_three_same_nonzero_direction": str(same_nonzero_direction).lower(),
                "median_absolute_delta_at_least_threshold": str(median_abs >= threshold).lower(),
                "absolute_share_threshold": f"{threshold:.12g}",
            }
        )

    summary = {
        "primary_comparison": primary,
        "endpoint": endpoint,
        "endpoint_description": endpoint_description,
        "independent_unit": analysis["independent_unit"],
        "absolute_share_threshold": threshold,
        "complete_pairs": len(rows),
        "saved_conclusion": analysis["conclusion"]["conclusion"],
        "saved_supported_material_change_units": analysis["conclusion"]["supported_material_change_units"],
    }
    return summary, output_rows


def general_counts(root: Path) -> dict[str, object]:
    loci = read_tsv(root / "data/processed/p2_loci/fixed_10_loci.tsv")
    correspondences = read_tsv(root / "data/processed/p2_comparison/supported_unit_correspondences.tsv")
    partitions = read_tsv(root / "data/processed/p6_rnaseq/aggregate/whole_array_partitions.tsv")
    protein_features = read_tsv(root / "data/processed/p6_protein/protein_feature_summary.tsv")
    protein_pairs = read_tsv(root / "data/processed/p6_protein/protein_pairwise_identity.tsv")
    relationships = Counter(row["relationship"] for row in correspondences)
    return {
        "fixed_loci": len(loci),
        "loci_by_subgroup": dict(sorted(Counter(row["analysis_subgroup"] for row in loci).items())),
        "supported_correspondences": {
            "within_group": relationships.get("within_group", 0),
            "cross_group": relationships.get("cross_group", 0),
            "total": len(correspondences),
        },
        "rnaseq_libraries": len({row["run_accession"] for row in partitions}),
        "whole_array_partitions_per_library": sorted(
            Counter(row["run_accession"] for row in partitions).values()
        ),
        "protein_loci": len(protein_features),
        "protein_pair_rows": len(protein_pairs),
    }


def dependency_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for distribution in ("numpy", "matplotlib", "Pillow", "biopython", "ViennaRNA"):
        try:
            versions[distribution] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[distribution] = None
    return versions


def render_figures(root: Path, output: Path, selection: str) -> None:
    renderer = root / "scripts/p6_revision_20260925_r6/figures/render_all.py"
    command = [
        sys.executable,
        str(renderer),
        "--project-root",
        str(root),
        "--output-dir",
        str(output / "figures"),
        "--figures",
        selection,
    ]
    subprocess.run(command, cwd=root, check=True)


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    output = args.output if args.output.is_absolute() else root / args.output
    output = output.resolve()
    inputs = validate_inputs(root, args.figures)
    if args.check_inputs:
        print(json.dumps({"status": "inputs_ok", "root": str(root), "inputs": inputs}, indent=2))
        return

    output.mkdir(parents=True, exist_ok=True)
    structure_summary, correlations = structure_statistics(root)
    temporal_summary, temporal_rows = temporal_statistics(root)
    counts = general_counts(root)

    tables = output / "tables"
    write_tsv(
        tables / "structure_identity_correlations.tsv",
        correlations,
        ["context", "population", "n", "identity_field", "bpp_field", "spearman_rho"],
    )
    write_tsv(
        tables / "primary_temporal_changes.tsv",
        temporal_rows,
        [
            "comparison",
            "endpoint",
            "unit",
            "cultures_n",
            "median_delta_share",
            "median_absolute_delta_share",
            "minimum_relative_change_percent",
            "maximum_relative_change_percent",
            "all_three_same_nonzero_direction",
            "median_absolute_delta_at_least_threshold",
            "absolute_share_threshold",
        ],
    )
    key_statistics = {
        "analysis_mode": "offline saved-table recomputation",
        "calculation_note": (
            "Counts, medians, correlations, and paired changes below were recalculated from TSV rows. "
            "Saved JSON summaries are used only for frozen analysis parameters and the saved analysis interpretation."
        ),
        "counts": counts,
        "structure": structure_summary,
        "temporal": temporal_summary,
    }
    (output / "key_statistics.json").write_text(
        json.dumps(key_statistics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    if args.figures != "none":
        render_figures(root, output, args.figures)

    environment = {
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "packages": dependency_versions(),
    }
    (output / "environment_versions.json").write_text(
        json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    p3 = json.loads(
        (root / "data/processed/p3_calibration/p3_calibration_summary.json").read_text(encoding="utf-8")
    )
    rna = json.loads((root / "data/processed/p6_rna/p6_rna_summary.json").read_text(encoding="utf-8"))
    generated_files = sorted(
        [path.relative_to(output).as_posix() for path in output.rglob("*") if path.is_file()]
        + ["reproduction_summary.json"]
    )
    try:
        output_display = output.relative_to(root).as_posix()
    except ValueError:
        output_display = "."
    summary = {
        "status": "complete",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "recommended_command": "python scripts/reproduce.py --output results/reproduced",
        "project_root": ".",
        "output": output_display,
        "figure_selection": args.figures,
        "environment_file": "environment_versions.json",
        "inputs": inputs,
        "consistency_checks": {
            "fixed_loci_expected_10": counts["fixed_loci"] == 10,
            "two_repeat_pairs_expected_252": structure_summary["all_pairs"] == 252,
            "two_repeat_partition_expected_35_plus_217": (
                structure_summary["supported_pairs"] == 35 and structure_summary["other_pairs"] == 217
            ),
            "rnaseq_libraries_expected_12": counts["rnaseq_libraries"] == 12,
            "protein_pair_rows_expected_45": counts["protein_pair_rows"] == 45,
        },
        "published_parameters_reused": {
            "p3_seed": p3["seed"],
            "p3_replicates": p3["replicates"],
            "p3_block_size_nt": p3["block_size_nt"],
            "rna_null_seed": rna["dinucleotide_null"]["seed"],
            "rna_null_replicates": rna["dinucleotide_null"]["replicates"],
            "rna_bpp_cutoff": rna["bpp_cutoff"],
            "viennarna_saved_version": rna["viennarna"]["version"],
            "viennarna_temperature_celsius": rna["viennarna"]["temperature_celsius"],
        },
        "randomness_in_this_run": "none; saved seeds are reported but their expensive analyses are not rerun",
        "not_run_by_this_entry_point": [
            "raw FASTQ download and validation",
            "fastp and Bowtie2 alignment",
            "full per-library RNA-seq aggregation",
            "ViennaRNA ensemble folding",
            "blockwise shuffled calibration",
            "remote protein retrieval",
        ],
        "generated_files": generated_files,
    }
    (output / "reproduction_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": "complete", "output": str(output), "files": len(generated_files)}, indent=2))


if __name__ == "__main__":
    main()
