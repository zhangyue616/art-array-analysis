#!/usr/bin/env python3
"""Compare RNA ensemble structure across the fixed ten-locus unit annotation."""

import csv
import datetime as dt
import importlib.util
import json
import math
import random
import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path


SEED = 20260924
N_REPLICATES = 64
BPP_CUTOFF = 1e-6
CONTEXTS = ("two_repeat", "one_repeat")
PAIR_CLASSES = ("spacer_spacer", "repeat_repeat", "mixed", "all")
P0_GROUP = "P0_seed7"
PB_GROUP = "PB50_near3"


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_tsv(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def format_value(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return "{:.8f}".format(value)
    return str(value)


def write_tsv(path, rows, columns=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        columns = list(rows[0]) if rows else []
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: format_value(row.get(column, "")) for column in columns})


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def read_fasta(path):
    records = {}
    header = None
    chunks = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header.split()[0]] = "".join(chunks).upper()
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records[header.split()[0]] = "".join(chunks).upper()
    return records


def median(values):
    return statistics.median(values) if values else None


def nearest_rank(values, probability):
    ordered = sorted(values)
    rank = max(1, int(math.ceil(probability * len(ordered))))
    return ordered[rank - 1]


def average_ranks(values):
    ordered = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and ordered[end][1] == ordered[start][1]:
            end += 1
        rank = (start + 1 + end) / 2.0
        for index in range(start, end):
            ranks[ordered[index][0]] = rank
        start = end
    return ranks


def pearson(left, right):
    if len(left) < 2:
        return None
    left_mean = statistics.mean(left)
    right_mean = statistics.mean(right)
    numerator = sum((a - left_mean) * (b - right_mean) for a, b in zip(left, right))
    left_ss = sum((a - left_mean) ** 2 for a in left)
    right_ss = sum((b - right_mean) ** 2 for b in right)
    if left_ss == 0 or right_ss == 0:
        return None
    return numerator / math.sqrt(left_ss * right_ss)


def spearman(left, right):
    if len(left) != len(right) or len(left) < 2:
        return None
    return pearson(average_ranks(left), average_ranks(right))


def dinucleotide_counts(sequence):
    return Counter(sequence[index:index + 2] for index in range(len(sequence) - 1))


def encode_counts(counts):
    return ";".join("{}:{}".format(key, counts[key]) for key in sorted(counts))


def dinucleotide_shuffle(sequence, rng):
    if len(sequence) < 2:
        return sequence
    adjacency = defaultdict(list)
    for left, right in zip(sequence[:-1], sequence[1:]):
        adjacency[left].append(right)
    for values in adjacency.values():
        rng.shuffle(values)
    stack = [sequence[0]]
    trail = []
    while stack:
        current = stack[-1]
        if adjacency[current]:
            stack.append(adjacency[current].pop())
        else:
            trail.append(stack.pop())
    shuffled = "".join(reversed(trail))
    if len(shuffled) != len(sequence):
        raise AssertionError("dinucleotide shuffle changed length")
    if dinucleotide_counts(shuffled) != dinucleotide_counts(sequence):
        raise AssertionError("dinucleotide shuffle changed counts")
    return shuffled


def pair_class(left_region, right_region):
    if left_region == "repeat" and right_region == "repeat":
        return "repeat_repeat"
    if left_region == "spacer" and right_region == "spacer":
        return "spacer_spacer"
    return "mixed"


def fold_sequence(RNA, sequence_dna, regions, store_details):
    sequence_rna = sequence_dna.replace("T", "U")
    model = RNA.md()
    compound = RNA.fold_compound(sequence_rna, model)
    mfe_structure, mfe = compound.mfe()
    compound.exp_params_rescale(mfe)
    partition_structure, ensemble_energy = compound.pf()
    bpp = compound.bpp()
    centroid_structure, centroid_distance = compound.centroid()
    ensemble_diversity = compound.mean_bp_distance()

    edge_by_class = {name: {} for name in PAIR_CLASSES}
    mass = {name: 0.0 for name in PAIR_CLASSES}
    pairedness = [{name: 0.0 for name in PAIR_CLASSES} for _ in sequence_dna]
    stored_edges = []
    for left in range(len(sequence_dna)):
        for right in range(left + 1, len(sequence_dna)):
            probability = float(bpp[left + 1][right + 1])
            if probability < BPP_CUTOFF:
                continue
            category = pair_class(regions[left], regions[right])
            key = (left, right)
            edge_by_class[category][key] = probability
            edge_by_class["all"][key] = probability
            mass[category] += probability
            mass["all"] += probability
            pairedness[left][category] += probability
            pairedness[right][category] += probability
            pairedness[left]["all"] += probability
            pairedness[right]["all"] += probability
            if store_details:
                stored_edges.append({
                    "left_position_1based": left + 1,
                    "right_position_1based": right + 1,
                    "left_region": regions[left],
                    "right_region": regions[right],
                    "pair_class": category,
                    "probability": probability,
                })
    return {
        "sequence_dna": sequence_dna,
        "regions": regions,
        "mfe_structure": mfe_structure,
        "mfe_kcal_mol": float(mfe),
        "partition_structure": partition_structure,
        "ensemble_free_energy_kcal_mol": float(ensemble_energy),
        "centroid_structure": centroid_structure,
        "centroid_distance": float(centroid_distance),
        "ensemble_diversity": float(ensemble_diversity),
        "edge_by_class": edge_by_class,
        "mass": mass,
        "pairedness": pairedness,
        "stored_edges": stored_edges,
    }


def position_to_column(gapped):
    mapping = {}
    position = 0
    for column, character in enumerate(gapped):
        if character != "-":
            mapping[position] = column
            position += 1
    return mapping


def remap_edges(edges, mapping):
    return {(mapping[left], mapping[right]): value for (left, right), value in edges.items()}


def bpp_overlap(left, right):
    denominator = sum(left.values()) + sum(right.values())
    if denominator == 0:
        return None
    shared = set(left).intersection(right)
    return 2.0 * sum(min(left[key], right[key]) for key in shared) / denominator


def score_pair(p2, aligner, left_fold, right_fold):
    evidence = p2.alignment_evidence(aligner, left_fold["sequence_dna"], right_fold["sequence_dna"])
    left_map = position_to_column(evidence["target_gapped"])
    right_map = position_to_column(evidence["query_gapped"])
    scores = {
        "context_global_identity": evidence["identity"],
        "context_diagonal_identity": evidence["diagonal_identity"],
        "context_gap_columns": evidence["gap_columns"],
    }
    for category in PAIR_CLASSES:
        left_vector = remap_edges(left_fold["edge_by_class"][category], left_map)
        right_vector = remap_edges(right_fold["edge_by_class"][category], right_map)
        scores["{}_bpp_overlap".format(category)] = bpp_overlap(left_vector, right_vector)
        scores["left_{}_expected_pair_mass".format(category)] = sum(left_vector.values())
        scores["right_{}_expected_pair_mass".format(category)] = sum(right_vector.values())
    return scores


def normalized_pair(left_unit_id, right_unit_id, unit_by_id):
    left = unit_by_id[left_unit_id]
    right = unit_by_id[right_unit_id]
    if left["analysis_subgroup"] == P0_GROUP and right["analysis_subgroup"] == PB_GROUP:
        return left, right
    if right["analysis_subgroup"] == P0_GROUP and left["analysis_subgroup"] == PB_GROUP:
        return right, left
    raise AssertionError("not a cross-group unit pair")


def genomic_interval(base_map, candidate_id, start_relative, end_exclusive_relative):
    first = base_map[(candidate_id, start_relative)]
    last = base_map[(candidate_id, end_exclusive_relative - 1)]
    first_genomic = int(first["genomic_position_1based"])
    last_genomic = int(last["genomic_position_1based"])
    return {
        "oriented_first_genomic_1based": first_genomic,
        "oriented_last_genomic_1based": last_genomic,
        "genomic_span_start_1based": min(first_genomic, last_genomic),
        "genomic_span_end_1based": max(first_genomic, last_genomic),
    }


def build_objects(loci, units, copies, sequences, base_map):
    locus_by_id = {row["candidate_id"]: row for row in loci}
    copy_by_id = {row["copy_id"]: row for row in copies}
    phase_rows = []
    phase_by_copy = {}
    for row in copies:
        if row["participates_in_seed_chain"].lower() != "true":
            continue
        offset = int(row["MarsHill14_best_offset_from_copy_start"])
        expected = 0 if row["analysis_subgroup"] == P0_GROUP else -2
        if offset != expected:
            raise AssertionError("unexpected repeat phase offset for {}".format(row["copy_id"]))
        calibrated = int(row["calibrated_copy_start_relative_to_RT"])
        anchor = calibrated + offset
        phase_by_copy[row["copy_id"]] = anchor
        phase_rows.append({
            "candidate_id": row["candidate_id"],
            "label": row["label"],
            "analysis_subgroup": row["analysis_subgroup"],
            "copy_id": row["copy_id"],
            "copy_ordinal_distal_to_proximal": row["copy_ordinal_distal_to_proximal"],
            "group_seed_word": row["group_seed_word"],
            "calibrated_copy_start_relative_to_RT": calibrated,
            "MarsHill14_best_offset_from_copy_start": offset,
            "operational_core14_anchor_start_relative_to_RT": anchor,
            "repeat22_start_relative_to_RT": anchor - 4,
            "repeat22_end_exclusive_relative_to_RT": anchor + 18,
            "MarsHill14_best_observed_sequence": row["MarsHill14_best_observed_sequence"],
            "MarsHill14_best_mismatches_in_window61": row["MarsHill14_best_mismatches_in_window61"],
        })

    object_rows = []
    unit_objects = {}
    for unit in units:
        candidate_id = unit["candidate_id"]
        locus = locus_by_id[candidate_id]
        relative_origin = int(locus["distal_anchor_relative_start_0based"])
        left_anchor = phase_by_copy[unit["left_copy_id"]]
        right_anchor = phase_by_copy[unit["right_copy_id"]]
        left_repeat_start = left_anchor - 4
        left_repeat_end = left_anchor + 18
        right_repeat_start = right_anchor - 4
        right_repeat_end = right_anchor + 18
        spacer_start = left_repeat_end
        spacer_end = right_repeat_start
        if spacer_end <= spacer_start:
            raise AssertionError("non-positive source-defined spacer")
        locus_sequence = sequences[candidate_id]

        def extract(start, end):
            start_index = start - relative_origin
            end_index = end - relative_origin
            if start_index < 0 or end_index > len(locus_sequence):
                raise AssertionError("structure object exceeds oriented locus")
            return locus_sequence[start_index:end_index]

        left_repeat = extract(left_repeat_start, left_repeat_end)
        spacer = extract(spacer_start, spacer_end)
        right_repeat = extract(right_repeat_start, right_repeat_end)
        if len(left_repeat) != 22 or len(right_repeat) != 22:
            raise AssertionError("repeat mask is not 22 nt")
        common = {
            "candidate_id": candidate_id,
            "label": unit["label"],
            "accession": unit["accession"],
            "analysis_subgroup": unit["analysis_subgroup"],
            "unit_id": unit["unit_id"],
            "unit_ordinal_distal_to_proximal": int(unit["unit_ordinal_distal_to_proximal"]),
            "left_copy_id": unit["left_copy_id"],
            "right_copy_id": unit["right_copy_id"],
            "left_calibrated_start_relative_to_RT": int(unit["left_copy_start_relative_to_RT"]),
            "right_calibrated_start_relative_to_RT": int(unit["right_copy_start_relative_to_RT"]),
            "left_core14_anchor_relative_to_RT": left_anchor,
            "right_core14_anchor_relative_to_RT": right_anchor,
            "phase_offset_nt": left_anchor - int(unit["left_copy_start_relative_to_RT"]),
            "left_repeat22_sequence": left_repeat,
            "spacer_sequence": spacer,
            "right_repeat22_sequence": right_repeat,
            "spacer_length_nt": len(spacer),
        }
        unit_objects[unit["unit_id"]] = {
            **common,
            "two_repeat": {
                "sequence": left_repeat + spacer + right_repeat,
                "regions": ["repeat"] * 22 + ["spacer"] * len(spacer) + ["repeat"] * 22,
                "start_relative": left_repeat_start,
                "end_relative": right_repeat_end,
            },
            "one_repeat": {
                "sequence": left_repeat + spacer,
                "regions": ["repeat"] * 22 + ["spacer"] * len(spacer),
                "start_relative": left_repeat_start,
                "end_relative": right_repeat_start,
            },
        }
        for context in CONTEXTS:
            obj = unit_objects[unit["unit_id"]][context]
            genomic = genomic_interval(base_map, candidate_id, obj["start_relative"], obj["end_relative"])
            object_rows.append({
                **{key: value for key, value in common.items() if not key.endswith("_sequence")},
                "context": context,
                "object_start_relative_to_RT": obj["start_relative"],
                "object_end_exclusive_relative_to_RT": obj["end_relative"],
                **genomic,
                "rt_strand": locus["rt_strand"],
                "sequence_length_nt": len(obj["sequence"]),
                "left_repeat_start_0based": 0,
                "left_repeat_end_exclusive_0based": 22,
                "spacer_start_0based": 22,
                "spacer_end_exclusive_0based": 22 + len(spacer),
                "right_repeat_start_0based": 22 + len(spacer) if context == "two_repeat" else "",
                "right_repeat_end_exclusive_0based": len(obj["sequence"]) if context == "two_repeat" else "",
                "sequence_dna": obj["sequence"],
            })
    return phase_rows, object_rows, unit_objects


def fold_observed(RNA, unit_objects):
    folds = {}
    ensemble_rows = []
    edge_rows = []
    position_rows = []
    for unit_id in sorted(unit_objects):
        unit = unit_objects[unit_id]
        for context in CONTEXTS:
            obj = unit[context]
            folded = fold_sequence(RNA, obj["sequence"], obj["regions"], store_details=True)
            folds[(unit_id, context)] = folded
            ensemble_rows.append({
                "candidate_id": unit["candidate_id"],
                "label": unit["label"],
                "analysis_subgroup": unit["analysis_subgroup"],
                "unit_id": unit_id,
                "unit_ordinal_distal_to_proximal": unit["unit_ordinal_distal_to_proximal"],
                "context": context,
                "sequence_length_nt": len(obj["sequence"]),
                "spacer_length_nt": unit["spacer_length_nt"],
                "mfe_kcal_mol": folded["mfe_kcal_mol"],
                "mfe_per_nt": folded["mfe_kcal_mol"] / len(obj["sequence"]),
                "ensemble_free_energy_kcal_mol": folded["ensemble_free_energy_kcal_mol"],
                "ensemble_free_energy_per_nt": folded["ensemble_free_energy_kcal_mol"] / len(obj["sequence"]),
                "ensemble_diversity": folded["ensemble_diversity"],
                "centroid_distance": folded["centroid_distance"],
                "spacer_spacer_expected_pair_mass": folded["mass"]["spacer_spacer"],
                "repeat_repeat_expected_pair_mass": folded["mass"]["repeat_repeat"],
                "mixed_expected_pair_mass": folded["mass"]["mixed"],
                "all_expected_pair_mass": folded["mass"]["all"],
                "mfe_structure": folded["mfe_structure"],
                "centroid_structure": folded["centroid_structure"],
                "sequence_rna": obj["sequence"].replace("T", "U"),
            })
            for edge in folded["stored_edges"]:
                edge_rows.append({
                    "unit_id": unit_id,
                    "label": unit["label"],
                    "context": context,
                    **edge,
                })
            for index, values in enumerate(folded["pairedness"]):
                position_rows.append({
                    "unit_id": unit_id,
                    "label": unit["label"],
                    "context": context,
                    "position_1based": index + 1,
                    "nucleotide": obj["sequence"][index],
                    "region": obj["regions"][index],
                    "pairedness_all": values["all"],
                    "pairedness_spacer_spacer": values["spacer_spacer"],
                    "pairedness_repeat_repeat": values["repeat_repeat"],
                    "pairedness_mixed": values["mixed"],
                    "unpaired_probability_approx": max(0.0, 1.0 - values["all"]),
                })
    return folds, ensemble_rows, edge_rows, position_rows


def score_all_pairs(p2, aligner, pairwise_rows, folds, unit_by_id, supported_keys):
    rows = []
    for source in pairwise_rows:
        if source["relationship"] != "cross_group":
            continue
        left_id = source["left_unit_id"]
        right_id = source["right_unit_id"]
        if left_id not in unit_by_id or right_id not in unit_by_id:
            continue
        key = frozenset((left_id, right_id))
        p0_unit, pb_unit = normalized_pair(left_id, right_id, unit_by_id)
        for context in CONTEXTS:
            scores = score_pair(p2, aligner, folds[(left_id, context)], folds[(right_id, context)])
            rows.append({
                "left_unit_id": left_id,
                "left_label": source["left_label"],
                "left_group": source["left_group"],
                "left_ordinal": int(source["left_ordinal"]),
                "right_unit_id": right_id,
                "right_label": source["right_label"],
                "right_group": source["right_group"],
                "right_ordinal": int(source["right_ordinal"]),
                "p0_unit_id": p0_unit["unit_id"],
                "p0_label": p0_unit["label"],
                "p0_ordinal": int(p0_unit["unit_ordinal_distal_to_proximal"]),
                "pb_unit_id": pb_unit["unit_id"],
                "pb_label": pb_unit["label"],
                "pb_ordinal": int(pb_unit["unit_ordinal_distal_to_proximal"]),
                "locus_pair_id": "{}__{}".format(p0_unit["label"], pb_unit["label"]),
                "context": context,
                "supported_correspondence": key in supported_keys,
                "same_ordinal": p0_unit["unit_ordinal_distal_to_proximal"] == pb_unit["unit_ordinal_distal_to_proximal"],
                "saved_core_to_core_global_identity": float(source["core_to_core_global_identity"]),
                "saved_trimmed_spacer_global_identity": float(source["trimmed_spacer_global_identity"]),
                **scores,
            })
    if len(rows) != 252 * len(CONTEXTS):
        raise AssertionError("unexpected cross-group pair score count")
    return rows


def summarize_locus_pairs(pair_rows):
    output = []
    grouped = defaultdict(list)
    for row in pair_rows:
        grouped[(row["context"], row["locus_pair_id"])].append(row)
    for (context, locus_pair_id), rows in sorted(grouped.items()):
        supported = [row for row in rows if row["supported_correspondence"]]
        mismatched = [row for row in rows if not row["supported_correspondence"]]
        output.append({
            "context": context,
            "locus_pair_id": locus_pair_id,
            "p0_label": rows[0]["p0_label"],
            "pb_label": rows[0]["pb_label"],
            "supported_pairs": len(supported),
            "misaligned_pairs": len(mismatched),
            "supported_spacer_bpp_overlap_median": median([row["spacer_spacer_bpp_overlap"] for row in supported]),
            "misaligned_spacer_bpp_overlap_median": median([row["spacer_spacer_bpp_overlap"] for row in mismatched]),
            "supported_minus_misaligned_spacer_bpp_overlap": (
                None if not supported else median([row["spacer_spacer_bpp_overlap"] for row in supported])
                - median([row["spacer_spacer_bpp_overlap"] for row in mismatched])
            ),
            "supported_repeat_bpp_overlap_median": median([row["repeat_repeat_bpp_overlap"] for row in supported]),
            "misaligned_repeat_bpp_overlap_median": median([row["repeat_repeat_bpp_overlap"] for row in mismatched]),
            "supported_mixed_bpp_overlap_median": median([row["mixed_bpp_overlap"] for row in supported]),
            "misaligned_mixed_bpp_overlap_median": median([row["mixed_bpp_overlap"] for row in mismatched]),
        })
    return output


def panel_metrics(rows):
    supported = [row for row in rows if row["supported_correspondence"]]
    mismatched = [row for row in rows if not row["supported_correspondence"]]
    result = {
        "supported_pairs": len(supported),
        "misaligned_pairs": len(mismatched),
    }
    for category in PAIR_CLASSES:
        field = "{}_bpp_overlap".format(category)
        result["supported_{}_median".format(field)] = median([row[field] for row in supported])
        result["misaligned_{}_median".format(field)] = median([row[field] for row in mismatched])
    for ordinal in (1, 2, 3):
        subset = [row for row in supported if row["p0_ordinal"] == ordinal and row["pb_ordinal"] == ordinal]
        result["u{}_supported_pairs".format(ordinal)] = len(subset)
        result["u{}_spacer_spacer_bpp_overlap_median".format(ordinal)] = median(
            [row["spacer_spacer_bpp_overlap"] for row in subset]
        )
    return result


def fold_shuffled_panel(RNA, unit_objects, shuffled_by_unit):
    folds = {}
    for unit_id in sorted(unit_objects):
        unit = unit_objects[unit_id]
        spacer = shuffled_by_unit[unit_id]
        left_repeat = unit["left_repeat22_sequence"]
        right_repeat = unit["right_repeat22_sequence"]
        sequences = {
            "two_repeat": left_repeat + spacer + right_repeat,
            "one_repeat": left_repeat + spacer,
        }
        regions = {
            "two_repeat": ["repeat"] * 22 + ["spacer"] * len(spacer) + ["repeat"] * 22,
            "one_repeat": ["repeat"] * 22 + ["spacer"] * len(spacer),
        }
        for context in CONTEXTS:
            folds[(unit_id, context)] = fold_sequence(
                RNA, sequences[context], regions[context], store_details=False
            )
    return folds


def score_supported_null(p2, aligner, supported_pairs, folds, unit_by_id, replicate_id):
    rows = []
    for pair in supported_pairs:
        left_id = pair["left_unit_id"]
        right_id = pair["right_unit_id"]
        p0_unit, pb_unit = normalized_pair(left_id, right_id, unit_by_id)
        for context in CONTEXTS:
            scores = score_pair(p2, aligner, folds[(left_id, context)], folds[(right_id, context)])
            rows.append({
                "replicate_id": replicate_id,
                "context": context,
                "left_unit_id": left_id,
                "right_unit_id": right_id,
                "p0_unit_id": p0_unit["unit_id"],
                "p0_label": p0_unit["label"],
                "p0_ordinal": int(p0_unit["unit_ordinal_distal_to_proximal"]),
                "pb_unit_id": pb_unit["unit_id"],
                "pb_label": pb_unit["label"],
                "pb_ordinal": int(pb_unit["unit_ordinal_distal_to_proximal"]),
                **scores,
            })
    return rows


def null_summary(observed_value, values):
    exceed = sum(value >= observed_value for value in values)
    return {
        "observed": observed_value,
        "null_median": median(values),
        "null_minimum": min(values),
        "null_maximum": max(values),
        "null_p95_nearest_rank": nearest_rank(values, 0.95),
        "null_at_least_observed": exceed,
        "plus_one_exceedance_fraction": (1 + exceed) / float(1 + len(values)),
        "observed_minus_null_median": observed_value - median(values),
    }


def main():
    started = dt.datetime.now(dt.timezone.utc)
    clock_start = time.perf_counter()
    root = Path(__file__).resolve().parents[2]
    output_dir = root / "data" / "processed" / "p6_rna"
    output_dir.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(root / "tools" / "p6_rna" / "python_packages"))
    sys.path.insert(0, str(root / "tools" / "python312_packages"))
    import RNA
    from Bio.Align import PairwiseAligner

    p2 = load_module(root / "scripts" / "p2_analysis" / "run_unified_locus_analysis.py", "p2_rna_library")
    aligner = p2.configure_aligner(PairwiseAligner)

    loci_dir = root / "data" / "processed" / "p2_loci"
    comparison_dir = root / "data" / "processed" / "p2_comparison"
    loci = sorted(read_tsv(loci_dir / "fixed_10_loci.tsv"), key=lambda row: int(row["locus_order"]))
    sequences = read_fasta(loci_dir / "oriented_anchor_to_partner_locus.fna")
    copies_all = read_tsv(comparison_dir / "unified_repeat_copies.tsv")
    chain_copy_ids = {
        row["copy_id"] for row in copies_all if row["participates_in_seed_chain"].lower() == "true"
    }
    units_all = read_tsv(comparison_dir / "unified_units.tsv")
    units = [
        row for row in units_all
        if row["left_copy_id"] in chain_copy_ids and row["right_copy_id"] in chain_copy_ids
    ]
    units.sort(key=lambda row: (int(next(locus["locus_order"] for locus in loci if locus["candidate_id"] == row["candidate_id"])), int(row["unit_ordinal_distal_to_proximal"])))
    unit_by_id = {row["unit_id"]: row for row in units}
    if len(loci) != 10 or len(chain_copy_ids) != 47 or len(units) != 37:
        raise AssertionError("fixed object count mismatch")

    base_map = {}
    with open(loci_dir / "oriented_base_map.tsv", "r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            relative = int(row["oriented_relative_0based"])
            if -1200 <= relative <= 0:
                base_map[(row["candidate_id"], relative)] = row

    phase_rows, object_rows, unit_objects = build_objects(loci, units, copies_all, sequences, base_map)
    if len(phase_rows) != 47 or len(object_rows) != 74:
        raise AssertionError("phase or sequence object count mismatch")
    write_tsv(output_dir / "repeat_phase_audit.tsv", phase_rows)
    write_tsv(output_dir / "unit_sequence_objects.tsv", object_rows)

    print("folding 37 observed units in two fixed contexts", flush=True)
    observed_folds, ensemble_rows, edge_rows, position_rows = fold_observed(RNA, unit_objects)
    write_tsv(output_dir / "unit_structure_ensembles.tsv", ensemble_rows)
    write_tsv(output_dir / "unit_bpp_edges.tsv", edge_rows)
    write_tsv(output_dir / "unit_position_profiles.tsv", position_rows)

    supported_source = [
        row for row in read_tsv(comparison_dir / "supported_unit_correspondences.tsv")
        if row["relationship"] == "cross_group"
        and row["left_unit_id"] in unit_by_id and row["right_unit_id"] in unit_by_id
    ]
    if len(supported_source) != 35:
        raise AssertionError("supported pair count mismatch")
    supported_keys = {frozenset((row["left_unit_id"], row["right_unit_id"])) for row in supported_source}
    pairwise_source = read_tsv(comparison_dir / "unit_pairwise_identity.tsv")
    pair_rows = score_all_pairs(p2, aligner, pairwise_source, observed_folds, unit_by_id, supported_keys)
    write_tsv(output_dir / "unit_pair_structure_similarity.tsv", pair_rows)
    locus_rows = summarize_locus_pairs(pair_rows)
    write_tsv(output_dir / "locus_pair_structure_summary.tsv", locus_rows)

    observed_by_context = {}
    for context in CONTEXTS:
        context_rows = [row for row in pair_rows if row["context"] == context]
        observed_by_context[context] = panel_metrics(context_rows)
        supported = [row for row in context_rows if row["supported_correspondence"]]
        observed_by_context[context]["spearman_spacer_overlap_vs_saved_trimmed_identity"] = spearman(
            [row["spacer_spacer_bpp_overlap"] for row in supported],
            [row["saved_trimmed_spacer_global_identity"] for row in supported],
        )
        locus_effects = [
            row["supported_minus_misaligned_spacer_bpp_overlap"]
            for row in locus_rows if row["context"] == context and row["supported_pairs"] > 0
        ]
        observed_by_context[context]["locus_pairs_with_supported"] = len(locus_effects)
        observed_by_context[context]["locus_pair_effect_median"] = median(locus_effects)
        observed_by_context[context]["locus_pairs_positive_effect"] = sum(value > 0 for value in locus_effects)

    rng = random.Random(SEED)
    shuffle_rows = []
    null_pair_rows = []
    null_panel_rows = []
    validation_rows = []
    for replicate_id in range(1, N_REPLICATES + 1):
        shuffled_by_unit = {}
        identical_shuffles = 0
        for unit_id in sorted(unit_objects):
            original = unit_objects[unit_id]["spacer_sequence"]
            shuffled = dinucleotide_shuffle(original, rng)
            shuffled_by_unit[unit_id] = shuffled
            identical = shuffled == original
            identical_shuffles += int(identical)
            shuffle_rows.append({
                "replicate_id": replicate_id,
                "unit_id": unit_id,
                "label": unit_objects[unit_id]["label"],
                "unit_ordinal_distal_to_proximal": unit_objects[unit_id]["unit_ordinal_distal_to_proximal"],
                "spacer_length_nt": len(original),
                "original_first_base": original[0],
                "original_last_base": original[-1],
                "shuffled_first_base": shuffled[0],
                "shuffled_last_base": shuffled[-1],
                "dinucleotide_counts": encode_counts(dinucleotide_counts(original)),
                "identical_to_original": identical,
                "shuffled_spacer_sequence": shuffled,
            })
        null_folds = fold_shuffled_panel(RNA, unit_objects, shuffled_by_unit)
        scored = score_supported_null(p2, aligner, supported_source, null_folds, unit_by_id, replicate_id)
        null_pair_rows.extend(scored)
        for context in CONTEXTS:
            subset = [row for row in scored if row["context"] == context]
            row = {
                "replicate_id": replicate_id,
                "context": context,
                "supported_pairs": len(subset),
                "spacer_spacer_bpp_overlap_median": median([item["spacer_spacer_bpp_overlap"] for item in subset]),
                "repeat_repeat_bpp_overlap_median": median([item["repeat_repeat_bpp_overlap"] for item in subset]),
                "mixed_bpp_overlap_median": median([item["mixed_bpp_overlap"] for item in subset]),
                "all_bpp_overlap_median": median([item["all_bpp_overlap"] for item in subset]),
            }
            for ordinal in (1, 2, 3):
                ordinal_rows = [
                    item for item in subset
                    if item["p0_ordinal"] == ordinal and item["pb_ordinal"] == ordinal
                ]
                row["u{}_pairs".format(ordinal)] = len(ordinal_rows)
                row["u{}_spacer_spacer_bpp_overlap_median".format(ordinal)] = median(
                    [item["spacer_spacer_bpp_overlap"] for item in ordinal_rows]
                )
            null_panel_rows.append(row)
        validation_rows.append({
            "replicate_id": replicate_id,
            "units_shuffled": len(shuffled_by_unit),
            "all_lengths_preserved": all(len(shuffled_by_unit[unit_id]) == len(unit_objects[unit_id]["spacer_sequence"]) for unit_id in shuffled_by_unit),
            "all_dinucleotide_counts_preserved": all(dinucleotide_counts(shuffled_by_unit[unit_id]) == dinucleotide_counts(unit_objects[unit_id]["spacer_sequence"]) for unit_id in shuffled_by_unit),
            "all_repeat_sequences_unchanged": True,
            "one_shared_unit_panel_used_for_all_supported_pairs": True,
            "identical_shuffles": identical_shuffles,
            "folded_objects": len(null_folds),
            "pair_scores": len(scored),
        })
        if replicate_id == 1 or replicate_id % 8 == 0:
            two_row = next(row for row in null_panel_rows if row["replicate_id"] == replicate_id and row["context"] == "two_repeat")
            print(
                "replicate {}/{}: spacer-overlap median={:.4f}".format(
                    replicate_id, N_REPLICATES, two_row["spacer_spacer_bpp_overlap_median"]
                ),
                flush=True,
            )

    write_tsv(output_dir / "dinucleotide_shuffle_manifest.tsv", shuffle_rows)
    write_tsv(output_dir / "dinucleotide_null_pair_scores.tsv", null_pair_rows)
    write_tsv(output_dir / "dinucleotide_null_replicates.tsv", null_panel_rows)
    write_tsv(output_dir / "dinucleotide_shuffle_validation.tsv", validation_rows)

    if len(shuffle_rows) != N_REPLICATES * 37:
        raise AssertionError("shuffle manifest count mismatch")
    if len(null_pair_rows) != N_REPLICATES * 35 * len(CONTEXTS):
        raise AssertionError("null pair score count mismatch")
    if len(null_panel_rows) != N_REPLICATES * len(CONTEXTS):
        raise AssertionError("null panel count mismatch")
    if any(not row["all_lengths_preserved"] for row in validation_rows):
        raise AssertionError("shuffle length validation failed")
    if any(not row["all_dinucleotide_counts_preserved"] for row in validation_rows):
        raise AssertionError("shuffle dinucleotide validation failed")

    null_comparison = {}
    for context in CONTEXTS:
        observed = observed_by_context[context]
        replicate_subset = [row for row in null_panel_rows if row["context"] == context]
        metrics = {
            "spacer_spacer_bpp_overlap_median": observed["supported_spacer_spacer_bpp_overlap_median"],
            "u1_spacer_spacer_bpp_overlap_median": observed["u1_spacer_spacer_bpp_overlap_median"],
            "u2_spacer_spacer_bpp_overlap_median": observed["u2_spacer_spacer_bpp_overlap_median"],
            "repeat_repeat_bpp_overlap_median": observed["supported_repeat_repeat_bpp_overlap_median"],
            "mixed_bpp_overlap_median": observed["supported_mixed_bpp_overlap_median"],
        }
        null_comparison[context] = {}
        for metric, observed_value in metrics.items():
            values = [row[metric] for row in replicate_subset if row[metric] is not None]
            null_comparison[context][metric] = null_summary(observed_value, values)

    model = RNA.md()
    finished = dt.datetime.now(dt.timezone.utc)
    elapsed = time.perf_counter() - clock_start
    summary = {
        "analysis": "cross-group RNA ensemble structure comparison",
        "date": "2026-09-24",
        "fixed_loci": len(loci),
        "fixed_chain_copies": len(chain_copy_ids),
        "fixed_chain_units": len(units),
        "supported_cross_group_pairs": len(supported_source),
        "misaligned_cross_group_pairs": 252 - len(supported_source),
        "sequence_contexts": list(CONTEXTS),
        "repeat_phase_offsets": {P0_GROUP: 0, PB_GROUP: -2},
        "bpp_cutoff": BPP_CUTOFF,
        "observed": observed_by_context,
        "dinucleotide_null": {
            "seed": SEED,
            "replicates": N_REPLICATES,
            "comparison": null_comparison,
        },
        "viennarna": {
            "version": RNA.__version__,
            "temperature_celsius": float(model.temperature),
            "dangles": int(model.dangles),
            "no_lonely_pairs": bool(model.noLP),
            "unique_multiloop_decomposition": bool(model.uniq_ML),
        },
        "runtime_seconds": elapsed,
        "interpretation_scope": "conditional on operational unit annotation and two source-derived sequence cuts; not mature-RNA or functional validation",
    }
    write_json(output_dir / "p6_rna_summary.json", summary)
    receipt = {
        "status": "PASS",
        "started_utc": started.isoformat(),
        "finished_utc": finished.isoformat(),
        "elapsed_seconds": elapsed,
        "python_version": sys.version.split()[0],
        "viennarna_version": RNA.__version__,
        "parameters": "scripts/p6_rna/ANALYSIS_PARAMETERS.md",
        "single_process": True,
        "scientific_input_files_modified": False,
        "validation": {
            "phase_rows": len(phase_rows),
            "sequence_objects": len(object_rows),
            "observed_folds": len(observed_folds),
            "cross_group_pair_context_rows": len(pair_rows),
            "shuffle_rows": len(shuffle_rows),
            "null_pair_context_rows": len(null_pair_rows),
            "null_panel_context_rows": len(null_panel_rows),
            "replicates_exact_dinucleotide": sum(row["all_dinucleotide_counts_preserved"] for row in validation_rows),
            "replicates_repeat_unchanged": sum(row["all_repeat_sequences_unchanged"] for row in validation_rows),
        },
    }
    write_json(output_dir / "run_receipt.json", receipt)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
