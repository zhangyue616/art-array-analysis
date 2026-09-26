#!/usr/bin/env python3
"""P3 panel-level spacer-order calibration for cross-group unit mapping."""

import collections
import datetime as dt
import importlib.util
import math
import os
from pathlib import Path
import random
import statistics


SEED = 20260924
N_REPLICATES = 64
BLOCK_SIZE = 50
BOUNDARY_TOLERANCE_NT = 10
TIE_TOLERANCE = 1e-12
P0_GROUP = "P0_seed7"
PB_GROUP = "PB50_near3"


def load_p2_module(root):
    path = root / "scripts" / "p2_analysis" / "run_unified_locus_analysis.py"
    spec = importlib.util.spec_from_file_location("p2_analysis_library", str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def nearest_rank(values, probability):
    ordered = sorted(values)
    rank = max(1, int(math.ceil(probability * len(ordered))))
    return ordered[rank - 1]


def distribution_summary(values, observed):
    defined = [value for value in values if value is not None]
    if not defined:
        return {
            "defined_replicates": 0,
            "missing_replicates": len(values),
            "median": None,
            "minimum": None,
            "maximum": None,
            "p95_nearest_rank": None,
            "null_at_least_observed": None,
            "panel_level_exceedance_fraction_plus_one": None,
        }
    exceed = sum(value >= observed for value in defined)
    return {
        "defined_replicates": len(defined),
        "missing_replicates": len(values) - len(defined),
        "median": statistics.median(defined),
        "minimum": min(defined),
        "maximum": max(defined),
        "p95_nearest_rank": nearest_rank(defined, 0.95),
        "null_at_least_observed": exceed,
        "panel_level_exceedance_fraction_plus_one": (1 + exceed) / float(1 + len(defined)),
    }


def normalized_cross_pair(left, right):
    if left["analysis_subgroup"] == P0_GROUP and right["analysis_subgroup"] == PB_GROUP:
        return left, right
    if right["analysis_subgroup"] == P0_GROUP and left["analysis_subgroup"] == PB_GROUP:
        return right, left
    raise AssertionError("pair is not a P0-by-PB50-group comparison")


def observed_metrics(p2, p2_dir, unit_by_id, retained_unit_ids):
    supported = []
    for row in p2.read_tsv(str(p2_dir / "supported_unit_correspondences.tsv")):
        if row["relationship"] != "cross_group":
            continue
        if row["left_unit_id"] not in retained_unit_ids or row["right_unit_id"] not in retained_unit_ids:
            continue
        supported.append(row)

    directed = []
    for row in p2.read_tsv(str(p2_dir / "unit_free_best_hits.tsv")):
        if row["query_group"] == row["target_group"]:
            continue
        if row["query_unit_id"] not in retained_unit_ids or row["best_target_unit_id"] not in retained_unit_ids:
            continue
        if row["sequence_and_coordinate_supported_correspondence"].lower() == "true":
            directed.append(row)

    locus_pairs_any = set()
    u1_pairs = set()
    u2_pairs = set()
    spacer_identities = []
    for row in supported:
        left = unit_by_id[row["left_unit_id"]]
        right = unit_by_id[row["right_unit_id"]]
        p0_unit, pb_unit = normalized_cross_pair(left, right)
        locus_pair = (p0_unit["candidate_id"], pb_unit["candidate_id"])
        locus_pairs_any.add(locus_pair)
        p0_ordinal = int(p0_unit["unit_ordinal_distal_to_proximal"])
        pb_ordinal = int(pb_unit["unit_ordinal_distal_to_proximal"])
        if p0_ordinal == 1 and pb_ordinal == 1:
            u1_pairs.add(locus_pair)
        if p0_ordinal == 2 and pb_ordinal == 2:
            u2_pairs.add(locus_pair)
        spacer_identities.append(float(row["trimmed_spacer_global_identity"]))

    result = {
        "supported_pair_count": len(supported),
        "supported_directed_count": len(directed),
        "cross_locus_pairs_with_any_support": len(locus_pairs_any),
        "u1_same_position_locus_pair_coverage": len(u1_pairs),
        "u2_same_position_locus_pair_coverage": len(u2_pairs),
        "supported_trimmed_spacer_identity_median": statistics.median(spacer_identities),
        "supported_trimmed_spacer_identity_minimum": min(spacer_identities),
        "supported_trimmed_spacer_identity_maximum": max(spacer_identities),
        "source": "P2 saved supported_unit_correspondences.tsv and unit_free_best_hits.tsv; restricted to 37 chain units",
    }
    if result["supported_pair_count"] != 35 or result["supported_directed_count"] != 70:
        raise AssertionError("saved P2 chain-only observed correspondence count changed")
    if result["u1_same_position_locus_pair_coverage"] != 15 or result["u2_same_position_locus_pair_coverage"] != 18:
        raise AssertionError("saved P2 U1/U2 observed coverage changed")
    return result


def build_shuffle_blocks(loci, units, locus_sequences, locus_relative_start):
    locus_order = {row["candidate_id"]: int(row["locus_order"]) for row in loci}
    ordered_units = sorted(
        units,
        key=lambda row: (
            locus_order[row["candidate_id"]],
            int(row["unit_ordinal_distal_to_proximal"]),
        ),
    )
    blocks = []
    occupied = collections.defaultdict(set)
    for unit in ordered_units:
        candidate_id = unit["candidate_id"]
        sequence = locus_sequences[candidate_id]
        relative_start = locus_relative_start[candidate_id]
        left = int(unit["left_copy_start_relative_to_RT"])
        right = int(unit["right_copy_start_relative_to_RT"])
        trimmed_start_rel = left + 40
        trimmed_end_rel = right - 26
        trimmed_start = trimmed_start_rel - relative_start
        trimmed_end = trimmed_end_rel - relative_start
        saved = unit["P0_trimmed_spacer_sequence"]
        if trimmed_end - trimmed_start != int(unit["P0_trimmed_spacer_length_nt"]):
            raise AssertionError("trimmed spacer coordinate/length mismatch for {}".format(unit["unit_id"]))
        if sequence[trimmed_start:trimmed_end] != saved:
            raise AssertionError("trimmed spacer sequence mismatch for {}".format(unit["unit_id"]))
        for position in range(trimmed_start, trimmed_end):
            if position in occupied[candidate_id]:
                raise AssertionError("trimmed spacer intervals overlap in {}".format(unit["label"]))
            occupied[candidate_id].add(position)
        for block_offset in range(0, len(saved), BLOCK_SIZE):
            block_sequence = saved[block_offset:block_offset + BLOCK_SIZE]
            counts = collections.Counter(block_sequence)
            blocks.append({
                "candidate_id": candidate_id,
                "label": unit["label"],
                "analysis_subgroup": unit["analysis_subgroup"],
                "unit_id": unit["unit_id"],
                "unit_ordinal_distal_to_proximal": int(unit["unit_ordinal_distal_to_proximal"]),
                "block_index_within_unit": block_offset // BLOCK_SIZE + 1,
                "locus_index_start_0based": trimmed_start + block_offset,
                "locus_index_end_exclusive_0based": trimmed_start + block_offset + len(block_sequence),
                "relative_start_to_RT": trimmed_start_rel + block_offset,
                "relative_end_exclusive_to_RT": trimmed_start_rel + block_offset + len(block_sequence),
                "length_nt": len(block_sequence),
                "A_count": counts.get("A", 0),
                "C_count": counts.get("C", 0),
                "G_count": counts.get("G", 0),
                "T_count": counts.get("T", 0),
                "original_sequence": block_sequence,
            })
    return ordered_units, blocks, occupied


def shuffle_panel(rng, loci, locus_sequences, blocks, occupied):
    shuffled_chars = {candidate_id: list(sequence) for candidate_id, sequence in locus_sequences.items()}
    all_block_composition_preserved = True
    all_block_lengths_preserved = True
    changed_positions = 0
    for block in blocks:
        candidate_id = block["candidate_id"]
        start = block["locus_index_start_0based"]
        end = block["locus_index_end_exclusive_0based"]
        before = list(locus_sequences[candidate_id][start:end])
        after = list(before)
        rng.shuffle(after)
        shuffled_chars[candidate_id][start:end] = after
        all_block_lengths_preserved = all_block_lengths_preserved and len(before) == len(after)
        all_block_composition_preserved = all_block_composition_preserved and collections.Counter(before) == collections.Counter(after)
        changed_positions += sum(left != right for left, right in zip(before, after))

    panel = {candidate_id: "".join(chars) for candidate_id, chars in shuffled_chars.items()}
    nontrimmed_unchanged = True
    for locus in loci:
        candidate_id = locus["candidate_id"]
        original = locus_sequences[candidate_id]
        replacement = panel[candidate_id]
        if len(original) != len(replacement):
            all_block_lengths_preserved = False
        mask = occupied[candidate_id]
        if any(original[index] != replacement[index] for index in range(len(original)) if index not in mask):
            nontrimmed_unchanged = False
    return panel, {
        "all_block_lengths_preserved": all_block_lengths_preserved,
        "all_block_mononucleotide_composition_preserved": all_block_composition_preserved,
        "all_nontrimmed_positions_unchanged": nontrimmed_unchanged,
        "changed_positions": changed_positions,
    }


def analyze_replicate(p2, aligner, loci, units, panel, locus_relative_start, locus_by_id, replicate_id):
    p0_loci = [row for row in loci if row["analysis_subgroup"] == P0_GROUP]
    pb_loci = [row for row in loci if row["analysis_subgroup"] == PB_GROUP]
    units_by_locus = collections.defaultdict(list)
    for unit in units:
        units_by_locus[unit["candidate_id"]].append(unit)
    for values in units_by_locus.values():
        values.sort(key=lambda row: int(row["unit_ordinal_distal_to_proximal"]))

    locus_pair_evidence = {}
    for p0_locus in p0_loci:
        for pb_locus in pb_loci:
            p0_id = p0_locus["candidate_id"]
            pb_id = pb_locus["candidate_id"]
            locus_pair_evidence[(p0_id, pb_id)] = p2.alignment_evidence(
                aligner, panel[p0_id], panel[pb_id]
            )

    replicate_units = []
    unit_by_id = {}
    for unit in units:
        row = dict(unit)
        candidate_id = row["candidate_id"]
        relative_start = locus_relative_start[candidate_id]
        left = int(row["left_copy_start_relative_to_RT"])
        right = int(row["right_copy_start_relative_to_RT"])
        left_index = left - relative_start
        right_index = right - relative_start
        row["core_to_core_sequence"] = panel[candidate_id][left_index:right_index]
        row["P0_trimmed_spacer_sequence"] = panel[candidate_id][left_index + 40:right_index - 26]
        if len(row["core_to_core_sequence"]) != int(row["core_to_core_length_nt"]):
            raise AssertionError("replicate unit length changed")
        if len(row["P0_trimmed_spacer_sequence"]) != int(row["P0_trimmed_spacer_length_nt"]):
            raise AssertionError("replicate trimmed spacer length changed")
        replicate_units.append(row)
        unit_by_id[row["unit_id"]] = row

    pair_lookup = {}
    unit_pair_alignment_count = 0
    for p0_unit in [row for row in replicate_units if row["analysis_subgroup"] == P0_GROUP]:
        for pb_unit in [row for row in replicate_units if row["analysis_subgroup"] == PB_GROUP]:
            core = p2.alignment_evidence(
                aligner, p0_unit["core_to_core_sequence"], pb_unit["core_to_core_sequence"]
            )
            spacer = p2.alignment_evidence(
                aligner,
                p0_unit["P0_trimmed_spacer_sequence"],
                pb_unit["P0_trimmed_spacer_sequence"],
            )
            pair = {
                "core_to_core_global_identity": core["identity"],
                "trimmed_spacer_global_identity": spacer["identity"],
            }
            pair_lookup[(p0_unit["unit_id"], pb_unit["unit_id"])] = pair
            pair_lookup[(pb_unit["unit_id"], p0_unit["unit_id"])] = pair
            unit_pair_alignment_count += 1

    directed_best = {}
    directed_rows = []
    for query in replicate_units:
        target_loci = pb_loci if query["analysis_subgroup"] == P0_GROUP else p0_loci
        for target_locus in target_loci:
            target_units = units_by_locus[target_locus["candidate_id"]]
            ranked = []
            for target in target_units:
                pair = pair_lookup[(query["unit_id"], target["unit_id"])]
                ranked.append((
                    float(pair["core_to_core_global_identity"]),
                    target["unit_id"],
                    target,
                    pair,
                ))
            ranked.sort(key=lambda value: (-value[0], value[1]))
            best_value = ranked[0][0]
            tied = [value for value in ranked if abs(value[0] - best_value) < TIE_TOLERANCE]
            runner_up = ranked[1][0] if len(ranked) > 1 else None
            unique = len(tied) == 1
            best = ranked[0]
            directed_best[(query["unit_id"], target_locus["candidate_id"])] = best[2]["unit_id"] if unique else None
            directed_rows.append({
                "replicate_id": replicate_id,
                "query_unit_id": query["unit_id"],
                "query_label": query["label"],
                "query_group": query["analysis_subgroup"],
                "query_ordinal": int(query["unit_ordinal_distal_to_proximal"]),
                "target_candidate_id": target_locus["candidate_id"],
                "target_label": target_locus["label"],
                "target_group": target_locus["analysis_subgroup"],
                "best_target_unit_id": best[2]["unit_id"],
                "best_target_ordinal": int(best[2]["unit_ordinal_distal_to_proximal"]),
                "best_core_to_core_identity": best_value,
                "runner_up_core_to_core_identity": runner_up,
                "best_runner_up_margin": best_value - runner_up if runner_up is not None else None,
                "best_is_unique": unique,
                "best_tie_count": len(tied),
            })

    for row in directed_rows:
        query = unit_by_id[row["query_unit_id"]]
        target = unit_by_id[row["best_target_unit_id"]]
        reverse = directed_best.get((target["unit_id"], query["candidate_id"]))
        reciprocal = bool(row["best_is_unique"] and reverse == query["unit_id"])
        if query["analysis_subgroup"] == P0_GROUP:
            evidence = locus_pair_evidence[(query["candidate_id"], target["candidate_id"])]
            boundary_map = evidence["target_to_query"]
        else:
            evidence = locus_pair_evidence[(target["candidate_id"], query["candidate_id"])]
            boundary_map = evidence["query_to_target"]
        query_start_rel = locus_relative_start[query["candidate_id"]]
        target_start_rel = locus_relative_start[target["candidate_id"]]
        query_left_index = int(query["left_copy_start_relative_to_RT"]) - query_start_rel
        query_right_index = int(query["right_copy_start_relative_to_RT"]) - query_start_rel
        mapped_left_index = boundary_map.get(query_left_index)
        mapped_right_index = boundary_map.get(query_right_index)
        mapped_left_rel = target_start_rel + mapped_left_index if mapped_left_index is not None else None
        mapped_right_rel = target_start_rel + mapped_right_index if mapped_right_index is not None else None
        left_offset = (
            mapped_left_rel - int(target["left_copy_start_relative_to_RT"])
            if mapped_left_rel is not None else None
        )
        right_offset = (
            mapped_right_rel - int(target["right_copy_start_relative_to_RT"])
            if mapped_right_rel is not None else None
        )
        coordinate_supported = (
            left_offset is not None
            and right_offset is not None
            and abs(left_offset) <= BOUNDARY_TOLERANCE_NT
            and abs(right_offset) <= BOUNDARY_TOLERANCE_NT
        )
        row.update({
            "reciprocal_unique_best": reciprocal,
            "mapped_query_left_boundary_in_target_relative": mapped_left_rel,
            "mapped_query_right_boundary_in_target_relative": mapped_right_rel,
            "left_boundary_offset_to_best_target_nt": left_offset,
            "right_boundary_offset_to_best_target_nt": right_offset,
            "whole_locus_boundary_support_within_10nt": coordinate_supported,
            "sequence_and_coordinate_supported_correspondence": bool(reciprocal and coordinate_supported),
        })

    supported_rows = []
    seen = set()
    for row in directed_rows:
        if not row["sequence_and_coordinate_supported_correspondence"]:
            continue
        pair_key = tuple(sorted((row["query_unit_id"], row["best_target_unit_id"])))
        if pair_key in seen:
            continue
        seen.add(pair_key)
        left = unit_by_id[pair_key[0]]
        right = unit_by_id[pair_key[1]]
        pair = pair_lookup[(left["unit_id"], right["unit_id"])]
        p0_unit, pb_unit = normalized_cross_pair(left, right)
        supported_rows.append({
            "replicate_id": replicate_id,
            "left_unit_id": left["unit_id"],
            "left_label": left["label"],
            "left_group": left["analysis_subgroup"],
            "left_ordinal": int(left["unit_ordinal_distal_to_proximal"]),
            "right_unit_id": right["unit_id"],
            "right_label": right["label"],
            "right_group": right["analysis_subgroup"],
            "right_ordinal": int(right["unit_ordinal_distal_to_proximal"]),
            "p0_candidate_id": p0_unit["candidate_id"],
            "p0_label": p0_unit["label"],
            "p0_ordinal": int(p0_unit["unit_ordinal_distal_to_proximal"]),
            "pb_candidate_id": pb_unit["candidate_id"],
            "pb_label": pb_unit["label"],
            "pb_ordinal": int(pb_unit["unit_ordinal_distal_to_proximal"]),
            "core_to_core_global_identity": pair["core_to_core_global_identity"],
            "trimmed_spacer_global_identity": pair["trimmed_spacer_global_identity"],
            "supporting_query_unit_id": row["query_unit_id"],
            "left_boundary_offset_to_best_target_nt": row["left_boundary_offset_to_best_target_nt"],
            "right_boundary_offset_to_best_target_nt": row["right_boundary_offset_to_best_target_nt"],
        })

    locus_pair_rows = []
    u1_coverage = 0
    u2_coverage = 0
    cross_locus_pairs_any = 0
    for p0_locus in p0_loci:
        for pb_locus in pb_loci:
            rows = [
                row for row in supported_rows
                if row["p0_candidate_id"] == p0_locus["candidate_id"]
                and row["pb_candidate_id"] == pb_locus["candidate_id"]
            ]
            ordinal_pairs = sorted({(row["p0_ordinal"], row["pb_ordinal"]) for row in rows})
            has_u1 = (1, 1) in ordinal_pairs
            has_u2 = (2, 2) in ordinal_pairs
            u1_coverage += int(has_u1)
            u2_coverage += int(has_u2)
            cross_locus_pairs_any += int(bool(rows))
            evidence = locus_pair_evidence[(p0_locus["candidate_id"], pb_locus["candidate_id"])]
            locus_pair_rows.append({
                "replicate_id": replicate_id,
                "p0_candidate_id": p0_locus["candidate_id"],
                "p0_label": p0_locus["label"],
                "pb_candidate_id": pb_locus["candidate_id"],
                "pb_label": pb_locus["label"],
                "anchor_through_RT_global_identity": evidence["identity"],
                "anchor_through_RT_gap_columns": evidence["gap_columns"],
                "anchor_through_RT_alignment_score": evidence["score"],
                "supported_pair_count": len(rows),
                "supported_ordinal_pairs_P0_to_PB": ";".join("{}:{}".format(*pair) for pair in ordinal_pairs),
                "u1_same_position_supported": has_u1,
                "u2_same_position_supported": has_u2,
                "supported_trimmed_spacer_identity_median": (
                    statistics.median(row["trimmed_spacer_global_identity"] for row in rows) if rows else None
                ),
            })

    spacer_values = [row["trimmed_spacer_global_identity"] for row in supported_rows]
    replicate_summary = {
        "replicate_id": replicate_id,
        "supported_pair_count": len(supported_rows),
        "supported_directed_count": sum(row["sequence_and_coordinate_supported_correspondence"] for row in directed_rows),
        "cross_locus_pairs_with_any_support": cross_locus_pairs_any,
        "u1_same_position_locus_pair_coverage": u1_coverage,
        "u2_same_position_locus_pair_coverage": u2_coverage,
        "supported_trimmed_spacer_identity_median": statistics.median(spacer_values) if spacer_values else None,
        "supported_trimmed_spacer_identity_minimum": min(spacer_values) if spacer_values else None,
        "supported_trimmed_spacer_identity_maximum": max(spacer_values) if spacer_values else None,
        "whole_locus_alignments_recomputed": len(locus_pair_evidence),
        "unit_pair_alignments_recomputed": unit_pair_alignment_count,
        "directed_free_best_selections_recomputed": len(directed_rows),
    }
    return replicate_summary, directed_rows, supported_rows, locus_pair_rows


def main():
    started = dt.datetime.now(dt.timezone.utc)
    root = Path(__file__).resolve().parents[2]
    p2_dir = root / "data" / "processed" / "p2_comparison"
    loci_dir = root / "data" / "processed" / "p2_loci"
    output_dir = root / "data" / "processed" / "p3_calibration"
    output_dir.mkdir(parents=True, exist_ok=True)

    p2 = load_p2_module(root)
    PairwiseAligner = p2.load_biopython(str(root))
    aligner = p2.configure_aligner(PairwiseAligner)

    loci = sorted(p2.read_tsv(str(loci_dir / "fixed_10_loci.tsv")), key=lambda row: int(row["locus_order"]))
    coordinate_rows = p2.read_tsv(str(loci_dir / "oriented_coordinate_map.tsv"))
    coordinate_by_id = {row["candidate_id"]: row for row in coordinate_rows}
    locus_by_id = {row["candidate_id"]: row for row in loci}
    full_sequences = p2.read_fasta(str(loci_dir / "oriented_anchor_to_partner_locus.fna"))

    locus_sequences = {}
    locus_relative_start = {}
    for locus in loci:
        candidate_id = locus["candidate_id"]
        relative_start = int(locus["distal_anchor_relative_start_0based"])
        rt_end = int(coordinate_by_id[candidate_id]["rt_relative_end_0based"])
        locus_sequences[candidate_id] = full_sequences[candidate_id][:rt_end - relative_start + 1]
        locus_relative_start[candidate_id] = relative_start

    all_copies = p2.read_tsv(str(p2_dir / "unified_repeat_copies.tsv"))
    chain_copy_ids = {
        row["copy_id"] for row in all_copies
        if row["participates_in_seed_chain"].lower() == "true"
    }
    all_units = p2.read_tsv(str(p2_dir / "unified_units.tsv"))
    units = [
        row for row in all_units
        if row["left_copy_id"] in chain_copy_ids and row["right_copy_id"] in chain_copy_ids
    ]
    retained_unit_ids = {row["unit_id"] for row in units}
    unit_by_id = {row["unit_id"]: row for row in units}
    if len(loci) != 10 or len(chain_copy_ids) != 47 or len(units) != 37:
        raise AssertionError("P3 fixed-object count mismatch")
    if sum(row["analysis_subgroup"] == P0_GROUP for row in units) != 28:
        raise AssertionError("P0 chain-unit count mismatch")
    if sum(row["analysis_subgroup"] == PB_GROUP for row in units) != 9:
        raise AssertionError("PB50-group chain-unit count mismatch")

    observed = observed_metrics(p2, p2_dir, unit_by_id, retained_unit_ids)
    ordered_units, blocks, occupied = build_shuffle_blocks(
        loci, units, locus_sequences, locus_relative_start
    )
    if len(blocks) == 0:
        raise AssertionError("no shuffle blocks")

    rng = random.Random(SEED)
    replicate_rows = []
    all_directed_rows = []
    all_supported_rows = []
    all_locus_pair_rows = []
    validation_rows = []
    for replicate_id in range(1, N_REPLICATES + 1):
        panel, shuffle_checks = shuffle_panel(rng, loci, locus_sequences, blocks, occupied)
        replicate, directed, supported, locus_pairs = analyze_replicate(
            p2, aligner, loci, ordered_units, panel, locus_relative_start,
            locus_by_id, replicate_id,
        )
        replicate_rows.append(replicate)
        all_directed_rows.extend(directed)
        all_supported_rows.extend(supported)
        all_locus_pair_rows.extend(locus_pairs)
        validation = {
            "replicate_id": replicate_id,
            "fixed_loci": len(loci),
            "fixed_chain_copies": len(chain_copy_ids),
            "fixed_chain_units": len(units),
            "shuffle_blocks": len(blocks),
            "shuffled_positions": sum(len(occupied[row["candidate_id"]]) for row in loci),
            "changed_positions": shuffle_checks["changed_positions"],
            **shuffle_checks,
            "one_shared_panel_state_used_for_all_21_locus_pairs": True,
            "whole_locus_alignments_recomputed": replicate["whole_locus_alignments_recomputed"],
            "unit_pair_alignments_recomputed": replicate["unit_pair_alignments_recomputed"],
            "directed_free_best_selections_recomputed": replicate["directed_free_best_selections_recomputed"],
        }
        validation_rows.append(validation)
        if replicate_id % 8 == 0 or replicate_id == 1:
            print(
                "replicate {}/{}: supported_pairs={} U1={} U2={} spacer_median={}".format(
                    replicate_id,
                    N_REPLICATES,
                    replicate["supported_pair_count"],
                    replicate["u1_same_position_locus_pair_coverage"],
                    replicate["u2_same_position_locus_pair_coverage"],
                    "NA" if replicate["supported_trimmed_spacer_identity_median"] is None else "{:.6f}".format(replicate["supported_trimmed_spacer_identity_median"]),
                ),
                flush=True,
            )

    if len(replicate_rows) != N_REPLICATES:
        raise AssertionError("replicate count mismatch")
    if any(not row["all_block_lengths_preserved"] for row in validation_rows):
        raise AssertionError("block length preservation failed")
    if any(not row["all_block_mononucleotide_composition_preserved"] for row in validation_rows):
        raise AssertionError("block composition preservation failed")
    if any(not row["all_nontrimmed_positions_unchanged"] for row in validation_rows):
        raise AssertionError("nontrimmed sequence changed")
    if any(row["whole_locus_alignments_recomputed"] != 21 for row in validation_rows):
        raise AssertionError("whole-locus alignment count mismatch")
    if any(row["unit_pair_alignments_recomputed"] != 252 for row in validation_rows):
        raise AssertionError("unit-pair alignment count mismatch")
    if any(row["directed_free_best_selections_recomputed"] != 147 for row in validation_rows):
        raise AssertionError("directed selection count mismatch")

    metric_specs = {
        "supported_pair_count": observed["supported_pair_count"],
        "u1_same_position_locus_pair_coverage": observed["u1_same_position_locus_pair_coverage"],
        "u2_same_position_locus_pair_coverage": observed["u2_same_position_locus_pair_coverage"],
        "supported_trimmed_spacer_identity_median": observed["supported_trimmed_spacer_identity_median"],
    }
    null_distributions = {}
    for metric, observed_value in metric_specs.items():
        values = [row[metric] for row in replicate_rows]
        null_distributions[metric] = {
            "observed": observed_value,
            **distribution_summary(values, observed_value),
        }

    summary = {
        "analysis": "P3 cross-group unit correspondence calibration under blockwise spacer-order shuffling",
        "parameters": "scripts/p3_calibration/ANALYSIS_PARAMETERS.md",
        "seed": SEED,
        "replicates": N_REPLICATES,
        "block_size_nt": BLOCK_SIZE,
        "fixed_loci": len(loci),
        "fixed_chain_copies": len(chain_copy_ids),
        "fixed_chain_units": len(units),
        "excluded_weak_edge_units": len(all_units) - len(units),
        "cross_group_locus_pairs_per_replicate": 21,
        "cross_group_unit_pairs_aligned_per_replicate": 252,
        "directed_free_best_selections_per_replicate": 147,
        "observed_reused_from_P2": observed,
        "null_distributions": null_distributions,
        "scope": "conditional on P2 chain annotation and geometry; not a copy-discovery validation",
    }

    p2.write_tsv(str(output_dir / "observed_metrics.tsv"), [observed], list(observed))
    p2.write_tsv(str(output_dir / "shuffle_block_manifest.tsv"), blocks, list(blocks[0]))
    p2.write_tsv(str(output_dir / "null_replicates.tsv"), replicate_rows, list(replicate_rows[0]))
    p2.write_tsv(str(output_dir / "null_directed_best_hits.tsv"), all_directed_rows, list(all_directed_rows[0]))
    supported_columns = [
        "replicate_id", "left_unit_id", "left_label", "left_group", "left_ordinal",
        "right_unit_id", "right_label", "right_group", "right_ordinal",
        "p0_candidate_id", "p0_label", "p0_ordinal", "pb_candidate_id", "pb_label",
        "pb_ordinal", "core_to_core_global_identity", "trimmed_spacer_global_identity",
        "supporting_query_unit_id", "left_boundary_offset_to_best_target_nt",
        "right_boundary_offset_to_best_target_nt",
    ]
    p2.write_tsv(str(output_dir / "null_supported_pairs.tsv"), all_supported_rows, supported_columns)
    p2.write_tsv(str(output_dir / "null_locus_pair_summary.tsv"), all_locus_pair_rows, list(all_locus_pair_rows[0]))
    p2.write_tsv(str(output_dir / "shuffle_validation_by_replicate.tsv"), validation_rows, list(validation_rows[0]))
    p2.write_json(str(output_dir / "p3_calibration_summary.json"), summary)

    validation_summary = {
        "status": "PASS",
        "replicates": len(replicate_rows),
        "replicates_with_exact_block_composition": sum(row["all_block_mononucleotide_composition_preserved"] for row in validation_rows),
        "replicates_with_nontrimmed_sequence_unchanged": sum(row["all_nontrimmed_positions_unchanged"] for row in validation_rows),
        "replicates_with_one_shared_panel_state": sum(row["one_shared_panel_state_used_for_all_21_locus_pairs"] for row in validation_rows),
        "whole_locus_alignments_recomputed_total": sum(row["whole_locus_alignments_recomputed"] for row in validation_rows),
        "unit_pair_alignments_recomputed_total": sum(row["unit_pair_alignments_recomputed"] for row in validation_rows),
        "directed_free_best_selections_recomputed_total": sum(row["directed_free_best_selections_recomputed"] for row in validation_rows),
        "observed_chain_only_supported_pairs_reused": observed["supported_pair_count"],
        "observed_chain_only_directed_supported_reused": observed["supported_directed_count"],
    }
    p2.write_json(str(output_dir / "validation_summary.json"), validation_summary)

    finished = dt.datetime.now(dt.timezone.utc)
    receipt = {
        "status": "PASS",
        "analysis": summary["analysis"],
        "started_utc": started.isoformat(),
        "finished_utc": finished.isoformat(),
        "elapsed_seconds": (finished - started).total_seconds(),
        "parameters_predeclared_before_run": True,
        "parameters_path": "scripts/p3_calibration/ANALYSIS_PARAMETERS.md",
        "random_seed": SEED,
        "replicates_completed": N_REPLICATES,
        "p2_files_modified": False,
        "observed_alignment_rerun": False,
        "observed_source": observed["source"],
        "weak_edge_sensitivity_rerun": False,
        "notes": [
            "Each replicate used one ten-locus shuffled panel for all 21 cross-group locus pairs.",
            "Full-locus boundary mapping, free unit selection, uniqueness, and reciprocity were recomputed within every replicate.",
            "The null conditions on fixed P2 chain-copy geometry and does not validate copy discovery.",
        ],
    }
    p2.write_json(str(output_dir / "run_receipt.json"), receipt)
    print("P3 calibration complete: {}".format(output_dir), flush=True)


if __name__ == "__main__":
    main()
