#!/usr/bin/env python3
"""Consolidated structural validation for completed P6 RNA outputs."""

from __future__ import division

import csv
import json
import os
import re
import statistics
from collections import Counter, defaultdict


def read_tsv(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def dinucleotide_counts(sequence):
    return Counter(sequence[index:index + 2] for index in range(len(sequence) - 1))


def encoded_counts(counter):
    return ";".join("{}:{}".format(key, counter[key]) for key in sorted(counter))


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    data = os.path.join(root, "data", "processed", "p6_rna")
    figures = os.path.join(root, "figures", "p6_rna")

    phase = read_tsv(os.path.join(data, "repeat_phase_audit.tsv"))
    objects = read_tsv(os.path.join(data, "unit_sequence_objects.tsv"))
    ensembles = read_tsv(os.path.join(data, "unit_structure_ensembles.tsv"))
    pairs = read_tsv(os.path.join(data, "unit_pair_structure_similarity.tsv"))
    locus = read_tsv(os.path.join(data, "locus_pair_structure_summary.tsv"))
    manifest = read_tsv(os.path.join(data, "dinucleotide_shuffle_manifest.tsv"))
    null_pairs = read_tsv(os.path.join(data, "dinucleotide_null_pair_scores.tsv"))
    null_replicates = read_tsv(os.path.join(data, "dinucleotide_null_replicates.tsv"))
    validation = read_tsv(os.path.join(data, "dinucleotide_shuffle_validation.tsv"))
    with open(os.path.join(data, "p6_rna_summary.json"), "r", encoding="utf-8") as handle:
        summary = json.load(handle)

    require(len(phase) == 47, "expected 47 phase rows")
    group_offsets = defaultdict(set)
    for row in phase:
        group_offsets[row["analysis_subgroup"]].add(int(row["MarsHill14_best_offset_from_copy_start"]))
    require(group_offsets["P0_seed7"] == set([0]), "reference phase offset mismatch")
    require(group_offsets["PB50_near3"] == set([-2]), "PB50-near phase offset mismatch")
    require(sum(row["analysis_subgroup"] == "P0_seed7" for row in phase) == 35, "reference phase count mismatch")
    require(sum(row["analysis_subgroup"] == "PB50_near3" for row in phase) == 12, "PB phase count mismatch")

    require(len(objects) == 74, "expected 74 sequence objects")
    unit_contexts = defaultdict(set)
    original_spacer = {}
    for row in objects:
        unit_contexts[row["unit_id"]].add(row["context"])
        if row["context"] == "two_repeat":
            start = int(row["spacer_start_0based"])
            end = int(row["spacer_end_exclusive_0based"])
            original_spacer[row["unit_id"]] = row["sequence_dna"][start:end]
    require(len(unit_contexts) == 37, "expected 37 units")
    require(all(value == set(["one_repeat", "two_repeat"]) for value in unit_contexts.values()),
            "each unit must have both contexts")
    require(len(ensembles) == 74, "expected 74 observed folds")
    require(set((row["unit_id"], row["context"]) for row in ensembles) ==
            set((row["unit_id"], row["context"]) for row in objects),
            "fold/object key mismatch")

    require(len(pairs) == 504, "expected 504 observed pair-context rows")
    for context in ("two_repeat", "one_repeat"):
        subset = [row for row in pairs if row["context"] == context]
        require(len(subset) == 252, "expected 252 rows per context")
        require(sum(row["supported_correspondence"] == "true" for row in subset) == 35,
                "expected 35 supported pairs per context")
        require(sum(row["supported_correspondence"] == "false" for row in subset) == 217,
                "expected 217 wrong-unit pairs per context")
        expected = statistics.median(float(row["spacer_spacer_bpp_overlap"]) for row in subset
                                     if row["supported_correspondence"] == "true")
        require(abs(expected - summary["observed"][context]["supported_spacer_spacer_bpp_overlap_median"]) < 1e-7,
                "observed summary mismatch for {}".format(context))
    require(len(locus) == 42, "expected 21 locus pairs x 2 contexts")
    require(sum(int(row["supported_pairs"]) > 0 for row in locus if row["context"] == "two_repeat") == 18,
            "expected 18 locus pairs with supported correspondence")

    require(len(manifest) == 2368, "expected 64 x 37 shuffle rows")
    manifest_by_replicate = defaultdict(list)
    for row in manifest:
        manifest_by_replicate[int(row["replicate_id"])].append(row)
        original = original_spacer[row["unit_id"]]
        shuffled = row["shuffled_spacer_sequence"]
        require(len(shuffled) == len(original) == int(row["spacer_length_nt"]), "shuffle length mismatch")
        require(shuffled[0] == original[0] and shuffled[-1] == original[-1], "Eulerian endpoints changed")
        require(encoded_counts(dinucleotide_counts(original)) == row["dinucleotide_counts"],
                "saved original dinucleotide counts mismatch")
        require(dinucleotide_counts(shuffled) == dinucleotide_counts(original), "dinucleotide counts changed")
    require(set(manifest_by_replicate) == set(range(1, 65)), "replicate IDs mismatch")
    require(all(len(rows) == 37 for rows in manifest_by_replicate.values()), "replicate unit count mismatch")

    require(len(null_pairs) == 4480, "expected 64 x 35 x 2 null scores")
    null_pair_counts = Counter(int(row["replicate_id"]) for row in null_pairs)
    require(all(null_pair_counts[index] == 70 for index in range(1, 65)), "null pair count mismatch")
    require(len(null_replicates) == 128, "expected 64 x 2 panel rows")
    require(set((int(row["replicate_id"]), row["context"]) for row in null_replicates) ==
            set((replicate, context) for replicate in range(1, 65)
                for context in ("two_repeat", "one_repeat")),
            "null panel/context keys mismatch")
    require(len(validation) == 64, "expected 64 validation rows")
    for row in validation:
        require(int(row["units_shuffled"]) == 37, "validation unit count mismatch")
        for field in (
            "all_lengths_preserved", "all_dinucleotide_counts_preserved",
            "all_repeat_sequences_unchanged", "one_shared_unit_panel_used_for_all_supported_pairs",
        ):
            require(row[field] == "true", "{} failed in replicate {}".format(field, row["replicate_id"]))
        require(int(row["folded_objects"]) == 74, "validation fold count mismatch")
        require(int(row["pair_scores"]) == 70, "validation pair score mismatch")

    figure_stems = [
        "figure_p6_1_supported_vs_wrong_units",
        "figure_p6_2_dinucleotide_panel_calibration",
        "figure_p6_3_representative_pair_maps",
    ]
    figure_rows = []
    for stem in figure_stems:
        for extension in ("pdf", "svg", "png"):
            path = os.path.join(figures, stem + "." + extension)
            require(os.path.isfile(path) and os.path.getsize(path) > 0, "missing figure {}".format(path))
        svg = os.path.join(figures, stem + ".svg")
        with open(svg, "r", encoding="utf-8") as handle:
            text = handle.read()
        width = float(re.search(r'width="([0-9.]+)pt"', text).group(1))
        font_sizes = [float(value) for value in re.findall(r"font-size:([0-9.]+)", text)]
        require(abs(width - 481.889764) < 0.01, "figure width is not 170 mm")
        require(font_sizes and min(font_sizes) >= 8.0, "figure contains text below 8 pt")
        figure_rows.append({
            "stem": stem,
            "svg_width_pt": width,
            "minimum_svg_font_size_pt": min(font_sizes),
        })

    require(os.path.isfile(os.path.join(root, "reports", "P6_RNA_REPORT.md")), "P6 report missing")

    receipt = {
        "status": "PASS",
        "date": "2026-09-24",
        "phase_rows": len(phase),
        "sequence_objects": len(objects),
        "observed_folds": len(ensembles),
        "observed_pair_context_rows": len(pairs),
        "supported_pairs_per_context": 35,
        "wrong_unit_pairs_per_context": 217,
        "shuffle_rows": len(manifest),
        "null_pair_context_rows": len(null_pairs),
        "null_panel_context_rows": len(null_replicates),
        "validated_replicates": len(validation),
        "all_shuffle_dinucleotides_recomputed_equal": True,
        "figures": figure_rows,
        "scientific_analysis_rerun_during_validation": False,
    }
    output = os.path.join(data, "validation_receipt.json")
    with open(output, "w", encoding="utf-8") as handle:
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
