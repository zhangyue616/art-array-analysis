#!/usr/bin/env python3
"""P2 unified annotation and locus/unit comparison for ten Staph ART loci.

The analysis follows ANALYSIS_PARAMETERS.md.  It deliberately uses the
source/feature assignments supplied by the P2 data lane and writes a separate
P2 result layer; no P0/P1 table is modified.
"""

import argparse
import csv
import itertools
import json
import math
import os
import statistics
import sys
from collections import Counter, defaultdict


SEED_LEN = 10
MARS_HILL_14 = "ATATGAATACGTAT"
WINDOW26_LEFT = 6
WINDOW26_RIGHT = 20
WINDOW61_LEFT = 20
WINDOW61_RIGHT = 41
MIN_GROUP_SEED_LOCI = 2


def project_root_from_script():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def load_biopython(project_root):
    local_packages = os.path.join(project_root, "tools", "python312_packages")
    if local_packages not in sys.path:
        sys.path.insert(0, local_packages)
    from Bio.Align import PairwiseAligner
    return PairwiseAligner


def read_tsv(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, columns):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: format_value(row.get(column, "")) for column in columns})


def format_value(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return "{:.6f}".format(value)
    if isinstance(value, (list, tuple, set)):
        return ";".join(str(item) for item in value)
    return str(value)


def write_json(path, value):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def read_fasta(path):
    records = {}
    header = None
    sequence = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header.split()[0]] = "".join(sequence).upper()
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
    if header is not None:
        records[header.split()[0]] = "".join(sequence).upper()
    return records


def write_fasta(path, records, width=80):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(">{}\n".format(header))
            for start in range(0, len(sequence), width):
                handle.write(sequence[start:start + width] + "\n")


def gc_fraction(sequence):
    if not sequence:
        return float("nan")
    return sum(base in "GC" for base in sequence) / float(len(sequence))


def shannon_entropy(word):
    counts = Counter(word)
    total = float(len(word))
    return -sum((count / total) * math.log(count / total, 2) for count in counts.values())


def hamming(left, right):
    return sum(a != b for a, b in zip(left, right))


def consensus(sequences):
    if not sequences:
        return ""
    width = min(len(sequence) for sequence in sequences)
    output = []
    for column in range(width):
        counts = Counter(sequence[column] for sequence in sequences)
        output.append(sorted(counts, key=lambda base: (-counts[base], base))[0])
    return "".join(output)


def profile_match(sequence, profile):
    return hamming(sequence, profile) if len(sequence) == len(profile) else None


def profile_match_count(sequence, profile):
    if len(sequence) != len(profile):
        return -1
    return sum(left == right for left, right in zip(sequence, profile))


def configure_aligner(PairwiseAligner):
    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    return aligner


def alignment_evidence(aligner, target, query):
    alignment = aligner.align(target, query)[0]
    coordinates = alignment.coordinates
    identical = 0
    aligned_columns = 0
    gap_columns = 0
    diagonal_columns = 0
    query_to_target = {}
    target_to_query = {}
    block_rows = []
    target_chunks = []
    query_chunks = []
    for block_index in range(coordinates.shape[1] - 1):
        target_start = int(coordinates[0, block_index])
        target_end = int(coordinates[0, block_index + 1])
        query_start = int(coordinates[1, block_index])
        query_end = int(coordinates[1, block_index + 1])
        target_step = target_end - target_start
        query_step = query_end - query_start
        if target_step and query_step:
            if target_step != query_step:
                raise AssertionError("unexpected unequal diagonal alignment block")
            kind = "aligned"
            target_piece = target[target_start:target_end]
            query_piece = query[query_start:query_end]
            identical += sum(a == b for a, b in zip(target_piece, query_piece))
            diagonal_columns += target_step
            aligned_columns += target_step
            for offset in range(target_step):
                target_to_query[target_start + offset] = query_start + offset
                query_to_target[query_start + offset] = target_start + offset
        elif target_step:
            kind = "query_gap"
            target_piece = target[target_start:target_end]
            query_piece = "-" * target_step
            gap_columns += target_step
            aligned_columns += target_step
        elif query_step:
            kind = "target_gap"
            target_piece = "-" * query_step
            query_piece = query[query_start:query_end]
            gap_columns += query_step
            aligned_columns += query_step
        else:
            continue
        target_chunks.append(target_piece)
        query_chunks.append(query_piece)
        block_rows.append({
            "block_index": block_index + 1,
            "block_type": kind,
            "target_start_0based": target_start,
            "target_end_exclusive_0based": target_end,
            "query_start_0based": query_start,
            "query_end_exclusive_0based": query_end,
            "alignment_columns": max(target_step, query_step),
        })
    identity = identical / float(aligned_columns) if aligned_columns else float("nan")
    diagonal_identity = identical / float(diagonal_columns) if diagonal_columns else float("nan")
    return {
        "score": float(alignment.score),
        "identical": identical,
        "aligned_columns": aligned_columns,
        "diagonal_columns": diagonal_columns,
        "gap_columns": gap_columns,
        "identity": identity,
        "diagonal_identity": diagonal_identity,
        "target_gapped": "".join(target_chunks),
        "query_gapped": "".join(query_chunks),
        "query_to_target": query_to_target,
        "target_to_query": target_to_query,
        "blocks": block_rows,
    }


def valid_seed_word(word, exact_count):
    if exact_count < 3 or len(set(word)) < 3:
        return False
    if any(base * 6 in word for base in "ACGT"):
        return False
    return True


def seed_words(sequence):
    counts = Counter(sequence[index:index + SEED_LEN] for index in range(len(sequence) - SEED_LEN + 1))
    return sorted(word for word, count in counts.items() if valid_seed_word(word, count))


def seed_occurrences(sequence, word, max_mismatch=1):
    raw = []
    for start in range(len(sequence) - len(word) + 1):
        observed = sequence[start:start + len(word)]
        mismatches = hamming(word, observed)
        if mismatches <= max_mismatch:
            raw.append({"start": start, "mismatches": mismatches, "observed": observed})
    if not raw:
        return []
    components = []
    current = []
    for item in raw:
        if current and item["start"] - current[-1]["start"] >= len(word):
            components.append(current)
            current = []
        current.append(item)
    if current:
        components.append(current)
    return [min(component, key=lambda item: (item["mismatches"], item["start"])) for component in components]


def chain_information_content(sequence, positions):
    windows = []
    for position in positions:
        left = position - WINDOW26_LEFT
        right = position + WINDOW26_RIGHT
        if left < 0 or right > len(sequence):
            return 0.0
        windows.append(sequence[left:right])
    background = Counter(sequence)
    total = float(len(sequence))
    information = 0.0
    for column in range(WINDOW26_LEFT + WINDOW26_RIGHT):
        counts = Counter(window[column] for window in windows)
        for base, count in counts.items():
            observed = count / float(len(windows))
            expected = background[base] / total
            if expected > 0:
                information += observed * math.log(observed / expected, 2)
    return information


def enumerate_chains(sequence, occurrences):
    positions = [item["start"] for item in occurrences]
    chains = set()

    def extend(chain, period, used_skip):
        if len(chain) >= 3:
            chains.add(tuple(chain))
        last = chain[-1]
        for candidate in positions:
            if candidate <= last:
                continue
            gap = candidate - last
            if gap < 60 or gap > 600:
                continue
            for multiplier in (1, 2):
                if multiplier == 2 and used_skip:
                    continue
                adjusted = gap / float(multiplier)
                if abs(adjusted - period) <= 0.30 * period:
                    new_period = (period * (len(chain) - 1) + adjusted) / float(len(chain))
                    extend(chain + [candidate], new_period, used_skip or multiplier == 2)

    for left_index in range(len(positions)):
        for right_index in range(left_index + 1, len(positions)):
            gap = positions[right_index] - positions[left_index]
            for multiplier in (1, 2):
                period = gap / float(multiplier)
                if 120 <= period <= 350:
                    extend([positions[left_index], positions[right_index]], period, multiplier == 2)

    rows = []
    for chain in chains:
        gaps = [chain[index + 1] - chain[index] for index in range(len(chain) - 1)]
        best_multipliers = []
        preliminary = statistics.median(gaps)
        for gap in gaps:
            options = [(abs(gap - preliminary), 1), (abs(gap / 2.0 - preliminary), 2)]
            best_multipliers.append(min(options)[1])
        if best_multipliers.count(2) > 1:
            continue
        adjusted_gaps = [gap / float(multiplier) for gap, multiplier in zip(gaps, best_multipliers)]
        period = statistics.median(adjusted_gaps)
        if not (120 <= period <= 350):
            continue
        if any(abs(adjusted - period) > 0.30 * period for adjusted in adjusted_gaps):
            continue
        normalized_deviation = sum(abs(adjusted - period) / period for adjusted in adjusted_gaps) / len(adjusted_gaps)
        information = chain_information_content(sequence, chain)
        score = (len(chain) - 1) * information
        rows.append({
            "positions": chain,
            "period": period,
            "gaps": gaps,
            "multipliers": best_multipliers,
            "normalized_spacing_deviation": normalized_deviation,
            "information_content": information,
            "score": score,
        })
    rows.sort(key=lambda row: (
        -row["score"], -len(row["positions"]), row["normalized_spacing_deviation"],
        -row["positions"][-1], row["positions"],
    ))
    return rows


def choose_group_seed(group_loci, interval_sequences):
    words = sorted(set(itertools.chain.from_iterable(seed_words(interval_sequences[locus["candidate_id"]]) for locus in group_loci)))
    ranking = []
    candidates = {}
    for word in words:
        by_locus = {}
        for locus in group_loci:
            sequence = interval_sequences[locus["candidate_id"]]
            exact_count = sum(sequence[index:index + SEED_LEN] == word for index in range(len(sequence) - SEED_LEN + 1))
            if not valid_seed_word(word, exact_count):
                continue
            occurrences = seed_occurrences(sequence, word, max_mismatch=1)
            chains = enumerate_chains(sequence, occurrences)
            if chains:
                by_locus[locus["candidate_id"]] = {
                    "occurrences": occurrences,
                    "chains": chains,
                    "best_chain": chains[0],
                    "exact_count": exact_count,
                }
        if len(by_locus) < MIN_GROUP_SEED_LOCI:
            continue
        mean_score = statistics.mean(item["best_chain"]["score"] for item in by_locus.values())
        ranking.append({
            "seed_word": word,
            "locus_support": len(by_locus),
            "mean_best_chain_score": mean_score,
            "seed_entropy": shannon_entropy(word),
        })
        candidates[word] = by_locus
    ranking.sort(key=lambda row: (-row["locus_support"], -row["mean_best_chain_score"], -row["seed_entropy"], row["seed_word"]))
    if not ranking:
        raise RuntimeError("no group seed passed the fixed chain rules")
    return ranking[0]["seed_word"], ranking, candidates


def cluster_mapped_occurrences(items, tolerance=10):
    clusters = []
    for item in sorted(items, key=lambda value: (value["mapped_medoid_index"], value["label"], value["interval_index"])):
        choices = []
        for cluster_index, cluster in enumerate(clusters):
            if item["candidate_id"] in {member["candidate_id"] for member in cluster["members"]}:
                continue
            distance = abs(item["mapped_medoid_index"] - cluster["center"])
            if distance <= tolerance:
                choices.append((distance, cluster_index))
        if choices:
            _, cluster_index = min(choices)
            clusters[cluster_index]["members"].append(item)
            clusters[cluster_index]["center"] = statistics.median(
                member["mapped_medoid_index"] for member in clusters[cluster_index]["members"]
            )
        else:
            clusters.append({"center": item["mapped_medoid_index"], "members": [item]})
    clusters.sort(key=lambda cluster: cluster["center"])
    return clusters


def infer_repeat_block(copy_windows):
    if not copy_windows:
        return {"left_offset": -6, "right_offset_exclusive": 20, "consensus": "", "column_consensus": []}
    fractions = []
    consensuses = []
    for column in range(61):
        counts = Counter(window[column] for window in copy_windows)
        base, count = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0]
        fractions.append(count / float(len(copy_windows)))
        consensuses.append(base)
    best = None
    for left in range(0, 21):
        for right in range(30, 62):
            low = sum(fractions[column] < 0.80 for column in range(left, right))
            if low > 1:
                continue
            candidate = (right - left, -abs((left + right - 1) / 2.0 - 24.5), -left, left, right)
            if best is None or candidate > best:
                best = candidate
    if best is None:
        left, right = 20, 30
    else:
        left, right = best[-2], best[-1]
    return {
        "left_offset": left - WINDOW61_LEFT,
        "right_offset_exclusive": right - WINDOW61_LEFT,
        "consensus": "".join(consensuses[left:right]),
        "column_consensus": fractions,
    }


def interval_overlap(left_start, left_end, right_start, right_end):
    return max(left_start, right_start) <= min(left_end, right_end)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=project_root_from_script())
    args = parser.parse_args()
    root = os.path.abspath(args.project_root)
    input_dir = os.path.join(root, "data", "processed", "p2_loci")
    output_dir = os.path.join(root, "data", "processed", "p2_comparison")
    os.makedirs(output_dir, exist_ok=True)

    PairwiseAligner = load_biopython(root)
    aligner = configure_aligner(PairwiseAligner)

    loci = read_tsv(os.path.join(input_dir, "fixed_10_loci.tsv"))
    coordinate_rows = read_tsv(os.path.join(input_dir, "oriented_coordinate_map.tsv"))
    coding_rows = read_tsv(os.path.join(input_dir, "coding_context.tsv"))
    interval_sequences = read_fasta(os.path.join(input_dir, "oriented_anchor_to_RT_interval.fna"))
    full_sequences = read_fasta(os.path.join(input_dir, "oriented_anchor_to_partner_locus.fna"))
    coordinate_by_id = {row["candidate_id"]: row for row in coordinate_rows}
    locus_by_id = {row["candidate_id"]: row for row in loci}
    coding_by_id = defaultdict(list)
    for row in coding_rows:
        coding_by_id[row["candidate_id"]].append(row)

    if len(loci) != 10 or len(interval_sequences) != 10 or len(full_sequences) != 10:
        raise AssertionError("P2 fixed set must contain exactly ten loci and both sequence objects")

    # Slice the full anchor-to-partner FASTA to anchor-start through RT-end.
    locus_sequences = {}
    locus_relative_start = {}
    for locus in loci:
        candidate_id = locus["candidate_id"]
        coordinate = coordinate_by_id[candidate_id]
        # oriented_anchor_to_partner_locus.fna begins at the selected distal
        # anchor, whereas full_locus_relative_start refers to the separate
        # 6-kb FASTA object.
        full_start = int(locus["distal_anchor_relative_start_0based"])
        rt_end = int(coordinate["rt_relative_end_0based"])
        full_sequence = full_sequences[candidate_id]
        length = rt_end - full_start + 1
        locus_sequences[candidate_id] = full_sequence[:length]
        locus_relative_start[candidate_id] = full_start
        expected_interval_length = int(locus["distal_anchor_to_RT_intervening_nt"])
        if len(interval_sequences[candidate_id]) != expected_interval_length:
            raise AssertionError("intergenic interval length mismatch for {}".format(locus["label"]))

    # All pairwise locus identities and group medoids.
    locus_pair_rows = []
    pair_cache = {}
    for left_index, left in enumerate(loci):
        for right in loci[left_index + 1:]:
            left_id = left["candidate_id"]
            right_id = right["candidate_id"]
            evidence = alignment_evidence(aligner, locus_sequences[left_id], locus_sequences[right_id])
            intergenic = alignment_evidence(aligner, interval_sequences[left_id], interval_sequences[right_id])
            pair_cache[(left_id, right_id)] = evidence
            pair_cache[(right_id, left_id)] = evidence
            locus_pair_rows.append({
                "left_candidate_id": left_id,
                "left_label": left["label"],
                "right_candidate_id": right_id,
                "right_label": right["label"],
                "relationship": "within_group" if left["analysis_subgroup"] == right["analysis_subgroup"] else "cross_group",
                "left_group": left["analysis_subgroup"],
                "right_group": right["analysis_subgroup"],
                "anchor_through_RT_left_length_nt": len(locus_sequences[left_id]),
                "anchor_through_RT_right_length_nt": len(locus_sequences[right_id]),
                "anchor_through_RT_global_identity": evidence["identity"],
                "anchor_through_RT_diagonal_identity": evidence["diagonal_identity"],
                "anchor_through_RT_gap_columns": evidence["gap_columns"],
                "intergenic_left_length_nt": len(interval_sequences[left_id]),
                "intergenic_right_length_nt": len(interval_sequences[right_id]),
                "intergenic_global_identity": intergenic["identity"],
                "intergenic_diagonal_identity": intergenic["diagonal_identity"],
                "intergenic_gap_columns": intergenic["gap_columns"],
                "method": "Biopython PairwiseAligner global affine; match=2,mismatch=-1,open=-5,extend=-1; identity includes gaps",
            })

    groups = defaultdict(list)
    for locus in loci:
        groups[locus["analysis_subgroup"]].append(locus)
    medoids = {}
    for group_name, group_loci in groups.items():
        mean_identities = {}
        for locus in group_loci:
            values = []
            for other in group_loci:
                if other["candidate_id"] == locus["candidate_id"]:
                    continue
                values.append(pair_cache[(locus["candidate_id"], other["candidate_id"])]["identity"])
            mean_identities[locus["candidate_id"]] = statistics.mean(values)
        medoids[group_name] = max(group_loci, key=lambda locus: (mean_identities[locus["candidate_id"]], -int(locus["locus_order"])))

    # Medoid alignment maps and coordinate-block evidence.
    medoid_alignments = {}
    alignment_summary_rows = []
    alignment_block_rows = []
    alignment_text = []
    for group_name, group_loci in groups.items():
        medoid = medoids[group_name]
        medoid_id = medoid["candidate_id"]
        for locus in group_loci:
            candidate_id = locus["candidate_id"]
            evidence = alignment_evidence(aligner, locus_sequences[medoid_id], locus_sequences[candidate_id])
            medoid_alignments[(group_name, candidate_id)] = evidence
            alignment_summary_rows.append({
                "analysis_subgroup": group_name,
                "medoid_candidate_id": medoid_id,
                "medoid_label": medoid["label"],
                "query_candidate_id": candidate_id,
                "query_label": locus["label"],
                "medoid_length_nt": len(locus_sequences[medoid_id]),
                "query_length_nt": len(locus_sequences[candidate_id]),
                "global_identity": evidence["identity"],
                "diagonal_identity": evidence["diagonal_identity"],
                "gap_columns": evidence["gap_columns"],
                "alignment_score": evidence["score"],
            })
            for block in evidence["blocks"]:
                row = dict(block)
                row.update({
                    "analysis_subgroup": group_name,
                    "medoid_candidate_id": medoid_id,
                    "medoid_label": medoid["label"],
                    "query_candidate_id": candidate_id,
                    "query_label": locus["label"],
                })
                alignment_block_rows.append(row)
            alignment_text.extend([
                "# group={} medoid={} query={} identity={:.6f}".format(group_name, medoid["label"], locus["label"], evidence["identity"]),
                ">{}|medoid".format(medoid["label"]),
                evidence["target_gapped"],
                ">{}|query".format(locus["label"]),
                evidence["query_gapped"],
                "",
            ])

    seed_ranking_rows = []
    group_seed = {}
    group_seed_candidates = {}
    for group_name, group_loci in groups.items():
        selected_word, ranking, candidates = choose_group_seed(group_loci, interval_sequences)
        group_seed[group_name] = selected_word
        group_seed_candidates[group_name] = candidates[selected_word]
        for rank, row in enumerate(ranking, 1):
            seed_ranking_rows.append({
                "analysis_subgroup": group_name,
                "rank": rank,
                "selected_group_seed": rank == 1,
                **row,
            })

    copy_rows = []
    annotation_rows = []
    seed_view_rows = []
    boundary_alternative_rows = []
    edge_scan_rows = []
    canonical_occurrence_audit_rows = []
    group_slot_details = {}

    for group_name, group_loci in groups.items():
        word = group_seed[group_name]
        medoid = medoids[group_name]
        medoid_id = medoid["candidate_id"]
        items = []
        per_locus_occurrences = {}
        for locus in group_loci:
            candidate_id = locus["candidate_id"]
            sequence = interval_sequences[candidate_id]
            occurrences = seed_occurrences(sequence, word, max_mismatch=1)
            chains = enumerate_chains(sequence, occurrences)
            participating = set(itertools.chain.from_iterable(chain["positions"] for chain in chains))
            per_locus_occurrences[candidate_id] = {item["start"]: item for item in occurrences}
            query_start_rel = -len(sequence)
            query_locus_rel_start = locus_relative_start[candidate_id]
            medoid_locus_rel_start = locus_relative_start[medoid_id]
            mapping = medoid_alignments[(group_name, candidate_id)]["query_to_target"]
            for occurrence in occurrences:
                interval_index = occurrence["start"]
                relative_start = query_start_rel + interval_index
                query_locus_index = relative_start - query_locus_rel_start
                mapped_locus_index = mapping.get(query_locus_index)
                if mapped_locus_index is None:
                    continue
                mapped_relative = medoid_locus_rel_start + mapped_locus_index
                item = {
                    "analysis_subgroup": group_name,
                    "candidate_id": candidate_id,
                    "label": locus["label"],
                    "interval_index": interval_index,
                    "relative_start": relative_start,
                    "mismatches": occurrence["mismatches"],
                    "observed": occurrence["observed"],
                    "participates_in_valid_chain": interval_index in participating,
                    "mapped_medoid_index": mapped_relative,
                }
                if item["participates_in_valid_chain"]:
                    items.append(item)
            # Preserve all seed views and their relation to the selected phase.
            for candidate_word in seed_words(sequence):
                candidate_occurrences = seed_occurrences(sequence, candidate_word, max_mismatch=1)
                candidate_chains = enumerate_chains(sequence, candidate_occurrences)
                if not candidate_chains:
                    continue
                best = candidate_chains[0]
                selected_positions = [entry["start"] for entry in occurrences]
                shifts = []
                for position in best["positions"]:
                    nearest = min(selected_positions, key=lambda value: abs(value - position))
                    if abs(nearest - position) <= 10:
                        shifts.append(nearest - position)
                constant_shift = int(round(statistics.median(shifts))) if shifts else None
                matched = sum(abs((position + (constant_shift or 0)) - nearest) <= 10 for position, nearest in zip(best["positions"], selected_positions[:len(best["positions"])])) if shifts else 0
                seed_view_rows.append({
                    "analysis_subgroup": group_name,
                    "candidate_id": candidate_id,
                    "label": locus["label"],
                    "candidate_seed_word": candidate_word,
                    "selected_group_seed_word": word,
                    "candidate_chain_positions_relative_to_RT": [query_start_rel + value for value in best["positions"]],
                    "candidate_chain_score": best["score"],
                    "candidate_chain_copy_count": len(best["positions"]),
                    "inferred_shift_to_selected_phase_nt": constant_shift,
                    "is_selected_seed_phase": candidate_word == word,
                    "view_audit_note": "raw seed phase retained; repeat/unit coordinates use group-selected seed phase",
                })

        clusters = cluster_mapped_occurrences(items, tolerance=10)
        support_threshold = int(math.ceil(len(group_loci) / 2.0))
        supported = [cluster for cluster in clusters if len({item["candidate_id"] for item in cluster["members"]}) >= support_threshold]
        if len(supported) < 3:
            raise AssertionError("fewer than three group-supported repeat slots for {}".format(group_name))
        # Keep the dominant periodic run if any isolated majority-supported word exists.
        best_run = None
        for left in range(len(supported)):
            for right in range(left + 3, len(supported) + 1):
                subset = supported[left:right]
                centers = [cluster["center"] for cluster in subset]
                gaps = [centers[index + 1] - centers[index] for index in range(len(centers) - 1)]
                period = statistics.median(gaps)
                if not (120 <= period <= 350):
                    continue
                if any(gap < 60 or gap > 600 or min(abs(gap - period), abs(gap / 2.0 - period)) > 0.30 * period for gap in gaps):
                    continue
                support_sum = sum(len({item["candidate_id"] for item in cluster["members"]}) for cluster in subset)
                candidate = (support_sum, len(subset), -sum(abs(gap - period) for gap in gaps), left, right)
                if best_run is None or candidate > best_run:
                    best_run = candidate
        if best_run is None:
            raise AssertionError("no periodic group-supported repeat run for {}".format(group_name))
        supported = supported[best_run[-2]:best_run[-1]]
        group_slot_details[group_name] = supported

        main_item_keys = {
            (member["candidate_id"], member["interval_index"])
            for cluster in supported for member in cluster["members"]
        }
        main_positions_by_locus = defaultdict(list)
        for cluster in supported:
            for member in cluster["members"]:
                main_positions_by_locus[member["candidate_id"]].append(member["relative_start"])
        for cluster_index, cluster in enumerate(clusters, 1):
            cluster_support = len({member["candidate_id"] for member in cluster["members"]})
            cluster_is_main = cluster in supported
            for member in cluster["members"]:
                native_main_positions = main_positions_by_locus[member["candidate_id"]]
                nearest_main = min(native_main_positions, key=lambda value: abs(value - member["relative_start"])) if native_main_positions else None
                nearest_distance = abs(nearest_main - member["relative_start"]) if nearest_main is not None else None
                selected = (member["candidate_id"], member["interval_index"]) in main_item_keys
                canonical_occurrence_audit_rows.append({
                    "analysis_subgroup": group_name,
                    "candidate_id": member["candidate_id"],
                    "label": member["label"],
                    "group_seed_word": word,
                    "observed_seed_word": member["observed"],
                    "seed_start_relative_to_RT": member["relative_start"],
                    "seed_mismatches": member["mismatches"],
                    "mapped_medoid_relative_start": member["mapped_medoid_index"],
                    "mapped_cluster_index": cluster_index,
                    "mapped_cluster_support_loci": cluster_support,
                    "group_majority_threshold": support_threshold,
                    "selected_in_main_annotation": selected,
                    "nearest_selected_copy_relative_to_RT": nearest_main,
                    "distance_to_nearest_selected_copy_nt": nearest_distance,
                    "can_be_additional_chain_copy_under_60nt_min_gap": nearest_distance is None or nearest_distance >= 60,
                    "interpretation": "main group-supported slot" if selected else (
                        "competing local seed register; not an additional copy under the 60-nt minimum gap"
                        if nearest_distance is not None and nearest_distance < 60
                        else "unselected low-support periodic seed occurrence"
                    ),
                })
                if not selected:
                    side = "distal_register" if nearest_main == min(native_main_positions) else "local_register"
                    boundary_alternative_rows.append({
                        "analysis_subgroup": group_name,
                        "candidate_id": member["candidate_id"],
                        "label": member["label"],
                        "side": side,
                        "alternative_type": "competing_seed_register" if nearest_distance is not None and nearest_distance < 60 else "unselected_seed_occurrence",
                        "expected_start_relative_to_RT": nearest_main,
                        "candidate_start_relative_to_RT": member["relative_start"],
                        "offset_from_expected_nt": member["relative_start"] - nearest_main if nearest_main is not None else None,
                        "profile_match_26": None,
                        "seed_mismatches": member["mismatches"],
                        "candidate_fixed26_sequence": interval_sequences[member["candidate_id"]][member["interval_index"] - 6:member["interval_index"] + 20],
                        "group_edge_candidate_loci": cluster_support,
                        "group_majority_threshold": support_threshold,
                        "status": "specific_unresolved_register_alternative_not_additional_copy" if nearest_distance is not None and nearest_distance < 60 else "low_support_unselected_occurrence",
                        "reason": "separated {} nt from selected copy; mapped support {}/{} loci".format(nearest_distance, cluster_support, len(group_loci)),
                    })

        # First pass: one occurrence per supported slot and locus.
        selected_by_locus = defaultdict(list)
        for slot_index, cluster in enumerate(supported, 1):
            by_candidate = {item["candidate_id"]: item for item in cluster["members"]}
            for locus in group_loci:
                item = by_candidate.get(locus["candidate_id"])
                if item is not None:
                    chosen = dict(item)
                    chosen["slot_index"] = slot_index
                    chosen["copy_class"] = "seed_exact" if item["mismatches"] == 0 else "seed_1_mismatch"
                    chosen["group_slot_support_loci"] = len(by_candidate)
                    selected_by_locus[locus["candidate_id"]].append(chosen)

        # A selected group slot missing in one locus receives the one-pass local profile search.
        for locus in group_loci:
            candidate_id = locus["candidate_id"]
            sequence = interval_sequences[candidate_id]
            current = selected_by_locus[candidate_id]
            windows = [sequence[item["interval_index"] - 6:item["interval_index"] + 20] for item in current if 6 <= item["interval_index"] <= len(sequence) - 20]
            locus_profile = consensus(windows)
            evidence = medoid_alignments[(group_name, candidate_id)]
            inverse = evidence["target_to_query"]
            query_locus_rel_start = locus_relative_start[candidate_id]
            medoid_rel_start = locus_relative_start[medoid_id]
            present_slots = {item["slot_index"] for item in current}
            for slot_index, cluster in enumerate(supported, 1):
                if slot_index in present_slots:
                    continue
                target_locus_index = int(round(cluster["center"] - medoid_rel_start))
                expected_query_locus_index = inverse.get(target_locus_index)
                if expected_query_locus_index is None:
                    continue
                expected_relative = query_locus_rel_start + expected_query_locus_index
                expected_interval = expected_relative + len(sequence)
                best = None
                for start in range(max(6, expected_interval - 40), min(len(sequence) - 20, expected_interval + 40) + 1):
                    window = sequence[start - 6:start + 20]
                    score = profile_match_count(window, locus_profile)
                    seed_mm = hamming(sequence[start:start + 10], word)
                    if score < 18 and seed_mm > 2:
                        continue
                    candidate = (score, -abs(start - expected_interval), -start, start, seed_mm, window)
                    if best is None or candidate > best:
                        best = candidate
                if best is not None:
                    start = best[3]
                    selected_by_locus[candidate_id].append({
                        "analysis_subgroup": group_name,
                        "candidate_id": candidate_id,
                        "label": locus["label"],
                        "interval_index": start,
                        "relative_start": start - len(sequence),
                        "mismatches": best[4],
                        "observed": sequence[start:start + 10],
                        "participates_in_valid_chain": False,
                        "mapped_medoid_index": cluster["center"],
                        "slot_index": slot_index,
                        "copy_class": "profile_rescued_expected_slot",
                        "group_slot_support_loci": len({item["candidate_id"] for item in cluster["members"]}),
                    })

        # Calculate group repeat block from phase-calibrated 61-nt windows.
        all_windows61 = []
        for locus in group_loci:
            sequence = interval_sequences[locus["candidate_id"]]
            for item in selected_by_locus[locus["candidate_id"]]:
                start = item["interval_index"]
                if start >= 20 and start + 41 <= len(sequence):
                    all_windows61.append(sequence[start - 20:start + 41])
        repeat_block = infer_repeat_block(all_windows61)

        for locus in group_loci:
            candidate_id = locus["candidate_id"]
            sequence = interval_sequences[candidate_id]
            selected = sorted(selected_by_locus[candidate_id], key=lambda item: item["relative_start"])
            profile26 = consensus([sequence[item["interval_index"] - 6:item["interval_index"] + 20] for item in selected])
            positions = [item["relative_start"] for item in selected]
            gaps = [positions[index + 1] - positions[index] for index in range(len(positions) - 1)]
            annotation_rows.append({
                "locus_order": locus["locus_order"],
                "candidate_id": candidate_id,
                "label": locus["label"],
                "accession": locus["accession"],
                "analysis_subgroup": group_name,
                "group_seed_word": word,
                "seed_chain_copy_count": len(selected),
                "promoted_degenerate_edge_copy_count": 0,
                "copy_count": len(selected),
                "copy_starts_relative_to_RT": positions,
                "unit_lengths_nt": gaps,
                "array_span_first_to_last_start_nt": positions[-1] - positions[0],
                "distal_copy_relative_start": positions[0],
                "proximal_copy_relative_start": positions[-1],
                "repeat_block_left_offset_from_seed": repeat_block["left_offset"],
                "repeat_block_right_offset_exclusive_from_seed": repeat_block["right_offset_exclusive"],
                "repeat_block_length_nt": repeat_block["right_offset_exclusive"] - repeat_block["left_offset"],
                "repeat_block_consensus": repeat_block["consensus"],
                "fixed26_locus_consensus": profile26,
                "main_annotation_status": "unique_group_supported_periodic_chain",
                "biological_endpoint_status": "operational_P2_boundary; biological_transcript_endpoint_unverified",
                "annotation_conflict_protein_ids": locus["annotation_conflict_protein_ids"],
                "assignment_status": locus["assignment_status"],
                "reported_art_support": locus["reported_art_support"],
            })
            for item in selected:
                start = item["interval_index"]
                relative_start = item["relative_start"]
                window26 = sequence[start - 6:start + 20]
                window61 = sequence[start - 20:start + 41]
                block_start = start + repeat_block["left_offset"]
                block_end = start + repeat_block["right_offset_exclusive"]
                repeat_sequence = sequence[block_start:block_end]
                left_flank = sequence[max(0, block_start - 20):block_start]
                right_flank = sequence[block_end:min(len(sequence), block_end + 20)]
                overlap_features = []
                block_relative_start = relative_start + repeat_block["left_offset"]
                block_relative_end = relative_start + repeat_block["right_offset_exclusive"] - 1
                for feature in coding_by_id[candidate_id]:
                    feature_start = int(feature["oriented_relative_start_0based"])
                    feature_end = int(feature["oriented_relative_end_0based"])
                    if interval_overlap(block_relative_start, block_relative_end, feature_start, feature_end):
                        overlap_features.append(feature["protein_id"] or feature["locus_tag"] or feature["feature_role"])
                old14_candidates = []
                for local_start in range(len(window61) - len(MARS_HILL_14) + 1):
                    observed14 = window61[local_start:local_start + len(MARS_HILL_14)]
                    old14_candidates.append((hamming(observed14, MARS_HILL_14), local_start, observed14))
                old14_best = min(old14_candidates)
                copy_rows.append({
                    "candidate_id": candidate_id,
                    "label": locus["label"],
                    "accession": locus["accession"],
                    "analysis_subgroup": group_name,
                    "copy_ordinal_distal_to_proximal": item["slot_index"],
                    "copy_id": "{}|copy{:02d}".format(locus["label"], item["slot_index"]),
                    "group_seed_word": word,
                    "raw_seed_word_observed": item["observed"],
                    "raw_seed_start_relative_to_RT": relative_start,
                    "calibrated_copy_start_relative_to_RT": relative_start,
                    "seed_phase_shift_nt": 0,
                    "seed_mismatches": item["mismatches"],
                    "copy_class": item["copy_class"],
                    "annotation_confidence_tier": "group_supported_periodic_seed_chain",
                    "biological_copy_status": "operational_sequence_annotation; transcript_endpoint_and_function_unverified",
                    "participates_in_seed_chain": item["participates_in_valid_chain"],
                    "group_slot_support_loci": item["group_slot_support_loci"],
                    "profile_match_26": profile_match_count(window26, profile26),
                    "fixed26_sequence": window26,
                    "window61_sequence": window61,
                    "repeat_block_relative_start": block_relative_start,
                    "repeat_block_relative_end": block_relative_end,
                    "repeat_block_sequence": repeat_sequence,
                    "left_flank20_sequence": left_flank,
                    "right_flank20_sequence": right_flank,
                    "MarsHill14_best_mismatches_in_window61": old14_best[0],
                    "MarsHill14_best_offset_from_copy_start": old14_best[1] - WINDOW61_LEFT,
                    "MarsHill14_best_observed_sequence": old14_best[2],
                    "coding_overlap_features": overlap_features,
                })

            # One-pass edge alternatives.  They are not promoted without group-majority support.
            usual_spacing = statistics.median(gaps)
            for side, expected_relative in (("distal", positions[0] - usual_spacing), ("proximal", positions[-1] + usual_spacing)):
                expected_index = int(round(expected_relative + len(sequence)))
                best = None
                for start in range(max(6, expected_index - 40), min(len(sequence) - 20, expected_index + 40) + 1):
                    if any(abs(start - item["interval_index"]) < 10 for item in selected):
                        continue
                    window = sequence[start - 6:start + 20]
                    profile_score = profile_match_count(window, profile26)
                    seed_mm = hamming(sequence[start:start + 10], word)
                    candidate = (profile_score, -seed_mm, -abs(start - expected_index), -start, start)
                    if best is None or candidate > best:
                        best = candidate
                if best is not None:
                    start = best[-1]
                    meets_floor = best[0] >= 18 or -best[1] <= 2
                    scan_row = {
                        "analysis_subgroup": group_name,
                        "candidate_id": candidate_id,
                        "label": locus["label"],
                        "side": side,
                        "expected_start_relative_to_RT": int(round(expected_relative)),
                        "candidate_start_relative_to_RT": start - len(sequence),
                        "offset_from_expected_nt": int(round(start - expected_index)),
                        "profile_match_26": best[0],
                        "seed_mismatches": -best[1],
                        "candidate_fixed26_sequence": sequence[start - 6:start + 20],
                        "meets_single_locus_candidate_floor": meets_floor,
                    }
                    edge_scan_rows.append(scan_row)
                    if not meets_floor:
                        continue
                    boundary_row = dict(scan_row)
                    boundary_row.update({
                        "alternative_type": "one_pass_expected_edge_candidate",
                        "status": "unresolved_pending_group_edge_support",
                        "reason": "one-pass expected edge search; candidate floor met but main boundary remains unresolved",
                    })
                    boundary_alternative_rows.append(boundary_row)

    # Promote no edge candidate unless the same mapped edge is present in a group majority.
    for group_name, group_loci in groups.items():
        threshold = int(math.ceil(len(group_loci) / 2.0))
        for side in ("distal", "proximal"):
            rows = [row for row in boundary_alternative_rows if row["analysis_subgroup"] == group_name and row["side"] == side]
            if not rows:
                continue
            support = len(rows)
            for row in rows:
                row["group_edge_candidate_loci"] = support
                row["group_majority_threshold"] = threshold
                row["status"] = "promoted_to_main_group_supported_degenerate_edge" if support >= threshold else "locus_specific_boundary_alternative"

    # Apply the rule predeclared in the P2 parameter record.  These copies are weaker than
    # the preliminary <=1-mismatch chain, but they meet the stated alternative
    # floor and occur at the same expected edge slot in a group majority.
    annotation_by_id = {row["candidate_id"]: row for row in annotation_rows}
    promoted_edge_rows = [
        row for row in boundary_alternative_rows
        if row.get("status") == "promoted_to_main_group_supported_degenerate_edge"
    ]
    for decision in promoted_edge_rows:
        candidate_id = decision["candidate_id"]
        locus = locus_by_id[candidate_id]
        sequence = interval_sequences[candidate_id]
        interval_rel_start = -len(sequence)
        relative_start = int(decision["candidate_start_relative_to_RT"])
        start = relative_start - interval_rel_start
        annotation = annotation_by_id[candidate_id]
        repeat_left = int(annotation["repeat_block_left_offset_from_seed"])
        repeat_right = int(annotation["repeat_block_right_offset_exclusive_from_seed"])
        block_start = start + repeat_left
        block_end = start + repeat_right
        block_relative_start = relative_start + repeat_left
        block_relative_end = relative_start + repeat_right - 1
        window26 = sequence[start - 6:start + 20]
        window61 = sequence[start - 20:start + 41]
        old14_candidates = []
        for local_start in range(len(window61) - len(MARS_HILL_14) + 1):
            observed14 = window61[local_start:local_start + len(MARS_HILL_14)]
            old14_candidates.append((hamming(observed14, MARS_HILL_14), local_start, observed14))
        old14_best = min(old14_candidates)
        overlap_features = []
        for feature in coding_by_id[candidate_id]:
            feature_start = int(feature["oriented_relative_start_0based"])
            feature_end = int(feature["oriented_relative_end_0based"])
            if interval_overlap(block_relative_start, block_relative_end, feature_start, feature_end):
                overlap_features.append(feature["protein_id"] or feature["locus_tag"] or feature["feature_role"])
        existing = [row for row in copy_rows if row["candidate_id"] == candidate_id]
        ordinal = len(existing) + 1 if decision["side"] == "proximal" else 1
        if decision["side"] == "distal":
            for row in existing:
                row["copy_ordinal_distal_to_proximal"] = int(row["copy_ordinal_distal_to_proximal"]) + 1
                row["copy_id"] = "{}|copy{:02d}".format(locus["label"], int(row["copy_ordinal_distal_to_proximal"]))
        copy_id = "{}|copy{:02d}".format(locus["label"], ordinal)
        copy_rows.append({
            "candidate_id": candidate_id,
            "label": locus["label"],
            "accession": locus["accession"],
            "analysis_subgroup": locus["analysis_subgroup"],
            "copy_ordinal_distal_to_proximal": ordinal,
            "copy_id": copy_id,
            "group_seed_word": annotation["group_seed_word"],
            "raw_seed_word_observed": sequence[start:start + 10],
            "raw_seed_start_relative_to_RT": relative_start,
            "calibrated_copy_start_relative_to_RT": relative_start,
            "seed_phase_shift_nt": 0,
            "seed_mismatches": int(decision["seed_mismatches"]),
            "copy_class": "group_supported_degenerate_edge",
            "annotation_confidence_tier": "low_confidence_group_supported_degenerate_edge_candidate",
            "biological_copy_status": "operational_edge_candidate; not an experimentally established biological copy",
            "participates_in_seed_chain": False,
            "group_slot_support_loci": int(decision["group_edge_candidate_loci"]),
            "profile_match_26": int(decision["profile_match_26"]),
            "fixed26_sequence": window26,
            "window61_sequence": window61,
            "repeat_block_relative_start": block_relative_start,
            "repeat_block_relative_end": block_relative_end,
            "repeat_block_sequence": sequence[block_start:block_end],
            "left_flank20_sequence": sequence[max(0, block_start - 20):block_start],
            "right_flank20_sequence": sequence[block_end:min(len(sequence), block_end + 20)],
            "MarsHill14_best_mismatches_in_window61": old14_best[0],
            "MarsHill14_best_offset_from_copy_start": old14_best[1] - WINDOW61_LEFT,
            "MarsHill14_best_observed_sequence": old14_best[2],
            "coding_overlap_features": overlap_features,
        })
        decision["promoted_copy_id"] = copy_id

    # Recompute locus-level main annotations after applying accepted edge copies.
    for annotation in annotation_rows:
        locus_copies = sorted(
            [row for row in copy_rows if row["candidate_id"] == annotation["candidate_id"]],
            key=lambda row: int(row["calibrated_copy_start_relative_to_RT"]),
        )
        positions = [int(row["calibrated_copy_start_relative_to_RT"]) for row in locus_copies]
        gaps = [positions[index + 1] - positions[index] for index in range(len(positions) - 1)]
        annotation["copy_count"] = len(locus_copies)
        annotation["copy_starts_relative_to_RT"] = positions
        annotation["unit_lengths_nt"] = gaps
        annotation["array_span_first_to_last_start_nt"] = positions[-1] - positions[0]
        annotation["distal_copy_relative_start"] = positions[0]
        annotation["proximal_copy_relative_start"] = positions[-1]
        if any(row["copy_class"] == "group_supported_degenerate_edge" for row in locus_copies):
            annotation["main_annotation_status"] = "group_supported_periodic_chain_plus_one_pass_degenerate_edge"
        annotation["promoted_degenerate_edge_copy_count"] = sum(
            row["copy_class"] == "group_supported_degenerate_edge" for row in locus_copies
        )

    alternatives_by_locus = defaultdict(list)
    for row in boundary_alternative_rows:
        alternatives_by_locus[row["candidate_id"]].append(row)
    for row in annotation_rows:
        alternatives = [
            item for item in alternatives_by_locus[row["candidate_id"]]
            if item.get("status") != "promoted_to_main_group_supported_degenerate_edge"
        ]
        row["boundary_alternative_count"] = len(alternatives)
        row["boundary_alternative_summary"] = [
            "{}:{}:{}".format(item.get("side", ""), item.get("candidate_start_relative_to_RT", ""), item.get("status", ""))
            for item in alternatives
        ]

    copy_rows.sort(key=lambda row: (int(locus_by_id[row["candidate_id"]]["locus_order"]), int(row["copy_ordinal_distal_to_proximal"])))
    annotation_rows.sort(key=lambda row: int(row["locus_order"]))

    # Units and sequences from calibrated copy starts.
    units = []
    unit_fastas = []
    copies_by_locus = defaultdict(list)
    for row in copy_rows:
        copies_by_locus[row["candidate_id"]].append(row)
    for locus in loci:
        candidate_id = locus["candidate_id"]
        sequence = interval_sequences[candidate_id]
        interval_rel_start = -len(sequence)
        copies = sorted(copies_by_locus[candidate_id], key=lambda row: int(row["calibrated_copy_start_relative_to_RT"]))
        for unit_index in range(len(copies) - 1):
            left = copies[unit_index]
            right = copies[unit_index + 1]
            left_rel = int(left["calibrated_copy_start_relative_to_RT"])
            right_rel = int(right["calibrated_copy_start_relative_to_RT"])
            left_index = left_rel - interval_rel_start
            right_index = right_rel - interval_rel_start
            unit_sequence = sequence[left_index:right_index]
            unflanked = sequence[left_index + 20:right_index - 6]
            trimmed = sequence[left_index + 40:right_index - 26]
            unit_id = "{}|unit{:02d}".format(locus["label"], unit_index + 1)
            row = {
                "candidate_id": candidate_id,
                "label": locus["label"],
                "accession": locus["accession"],
                "analysis_subgroup": locus["analysis_subgroup"],
                "unit_id": unit_id,
                "unit_ordinal_distal_to_proximal": unit_index + 1,
                "left_copy_id": left["copy_id"],
                "right_copy_id": right["copy_id"],
                "left_copy_start_relative_to_RT": left_rel,
                "right_copy_start_relative_to_RT": right_rel,
                "core_to_core_length_nt": len(unit_sequence),
                "core_to_core_GC": gc_fraction(unit_sequence),
                "core_to_core_sequence": unit_sequence,
                "between_fixed26_length_nt": len(unflanked),
                "between_fixed26_GC": gc_fraction(unflanked),
                "between_fixed26_sequence": unflanked,
                "P0_trimmed_spacer_length_nt": len(trimmed),
                "P0_trimmed_spacer_GC": gc_fraction(trimmed),
                "P0_trimmed_spacer_sequence": trimmed,
            }
            units.append(row)
            unit_fastas.extend([
                (unit_id + "|core_to_core", unit_sequence),
                (unit_id + "|P0_trimmed_spacer", trimmed),
            ])

    # Detailed cross-group medoid boundary map.  This records where a copy
    # boundary maps to an aligned base versus a gap, so proximal local collapse
    # is not inferred from ordinal labels.
    p0_medoid = medoids["P0_seed7"]
    near3_medoid = medoids["PB50_near3"]
    p0_id = p0_medoid["candidate_id"]
    near3_id = near3_medoid["candidate_id"]
    if int(p0_medoid["locus_order"]) < int(near3_medoid["locus_order"]):
        cross_medoid_evidence = pair_cache[(p0_id, near3_id)]
        near3_to_p0 = cross_medoid_evidence["query_to_target"]
        p0_to_near3 = cross_medoid_evidence["target_to_query"]
    else:
        cross_medoid_evidence = pair_cache[(near3_id, p0_id)]
        near3_to_p0 = cross_medoid_evidence["target_to_query"]
        p0_to_near3 = cross_medoid_evidence["query_to_target"]
    p0_start_rel = locus_relative_start[p0_id]
    near3_start_rel = locus_relative_start[near3_id]
    cross_copy_mapping_rows = []
    for source_locus, target_locus, mapping, source_start, target_start in (
        (p0_medoid, near3_medoid, p0_to_near3, p0_start_rel, near3_start_rel),
        (near3_medoid, p0_medoid, near3_to_p0, near3_start_rel, p0_start_rel),
    ):
        for copy in sorted(copies_by_locus[source_locus["candidate_id"]], key=lambda row: int(row["copy_ordinal_distal_to_proximal"])):
            source_relative = int(copy["calibrated_copy_start_relative_to_RT"])
            source_index = source_relative - source_start
            target_index = mapping.get(source_index)
            cross_copy_mapping_rows.append({
                "source_group": source_locus["analysis_subgroup"],
                "source_label": source_locus["label"],
                "source_copy_id": copy["copy_id"],
                "source_copy_ordinal": copy["copy_ordinal_distal_to_proximal"],
                "source_copy_relative_start": source_relative,
                "target_group": target_locus["analysis_subgroup"],
                "target_label": target_locus["label"],
                "mapped_target_relative_start": target_start + target_index if target_index is not None else None,
                "mapping_status": "aligned_base" if target_index is not None else "maps_to_gap",
                "whole_anchor_through_RT_global_identity": cross_medoid_evidence["identity"],
            })
    cross_unit_boundary_rows = []
    for source_locus, target_locus, mapping, source_start, target_start in (
        (p0_medoid, near3_medoid, p0_to_near3, p0_start_rel, near3_start_rel),
        (near3_medoid, p0_medoid, near3_to_p0, near3_start_rel, p0_start_rel),
    ):
        for unit in units_by_locus[source_locus["candidate_id"]] if 'units_by_locus' in locals() else [row for row in units if row["candidate_id"] == source_locus["candidate_id"]]:
            left_relative = int(unit["left_copy_start_relative_to_RT"])
            right_relative = int(unit["right_copy_start_relative_to_RT"])
            left_index = mapping.get(left_relative - source_start)
            right_index = mapping.get(right_relative - source_start)
            cross_unit_boundary_rows.append({
                "source_group": source_locus["analysis_subgroup"],
                "source_label": source_locus["label"],
                "source_unit_id": unit["unit_id"],
                "source_unit_ordinal": unit["unit_ordinal_distal_to_proximal"],
                "source_left_relative": left_relative,
                "source_right_relative": right_relative,
                "target_group": target_locus["analysis_subgroup"],
                "target_label": target_locus["label"],
                "mapped_left_target_relative": target_start + left_index if left_index is not None else None,
                "mapped_right_target_relative": target_start + right_index if right_index is not None else None,
                "boundary_mapping_status": "both_boundaries_aligned" if left_index is not None and right_index is not None else (
                    "left_only" if left_index is not None else ("right_only" if right_index is not None else "neither")
                ),
            })
    cross_medoid_block_rows = []
    for block in cross_medoid_evidence["blocks"]:
        row = dict(block)
        row.update({
            "target_group": p0_medoid["analysis_subgroup"],
            "target_label": p0_medoid["label"],
            "target_relative_start": p0_start_rel + int(block["target_start_0based"]),
            "target_relative_end_exclusive": p0_start_rel + int(block["target_end_exclusive_0based"]),
            "query_group": near3_medoid["analysis_subgroup"],
            "query_label": near3_medoid["label"],
            "query_relative_start": near3_start_rel + int(block["query_start_0based"]),
            "query_relative_end_exclusive": near3_start_rel + int(block["query_end_exclusive_0based"]),
        })
        if max(row["target_relative_end_exclusive"], row["query_relative_end_exclusive"]) >= -1250 and min(row["target_relative_start"], row["query_relative_start"]) <= 100:
            cross_medoid_block_rows.append(row)

    # Pairwise unit alignments, plus free best-hit mapping per locus pair.
    unit_pair_rows = []
    for left_index, left in enumerate(units):
        for right in units[left_index + 1:]:
            if left["candidate_id"] == right["candidate_id"]:
                continue
            full = alignment_evidence(aligner, left["core_to_core_sequence"], right["core_to_core_sequence"])
            trimmed = alignment_evidence(aligner, left["P0_trimmed_spacer_sequence"], right["P0_trimmed_spacer_sequence"])
            unit_pair_rows.append({
                "left_unit_id": left["unit_id"],
                "left_label": left["label"],
                "left_group": left["analysis_subgroup"],
                "left_ordinal": left["unit_ordinal_distal_to_proximal"],
                "right_unit_id": right["unit_id"],
                "right_label": right["label"],
                "right_group": right["analysis_subgroup"],
                "right_ordinal": right["unit_ordinal_distal_to_proximal"],
                "relationship": "within_group" if left["analysis_subgroup"] == right["analysis_subgroup"] else "cross_group",
                "core_to_core_global_identity": full["identity"],
                "core_to_core_diagonal_identity": full["diagonal_identity"],
                "core_to_core_gap_columns": full["gap_columns"],
                "trimmed_spacer_global_identity": trimmed["identity"],
                "trimmed_spacer_diagonal_identity": trimmed["diagonal_identity"],
                "trimmed_spacer_gap_columns": trimmed["gap_columns"],
            })

    pair_lookup = {}
    for row in unit_pair_rows:
        pair_lookup[(row["left_unit_id"], row["right_unit_id"])] = row
        pair_lookup[(row["right_unit_id"], row["left_unit_id"])] = row
    units_by_locus = defaultdict(list)
    for unit in units:
        units_by_locus[unit["candidate_id"]].append(unit)
    best_hit_rows = []
    directed_best = {}
    for query in units:
        for target_locus_id, target_units in units_by_locus.items():
            if target_locus_id == query["candidate_id"]:
                continue
            ranked = []
            for target in target_units:
                pair = pair_lookup[(query["unit_id"], target["unit_id"])]
                ranked.append((float(pair["core_to_core_global_identity"]), target["unit_id"], target, pair))
            ranked.sort(key=lambda value: (-value[0], value[1]))
            best_value = ranked[0][0]
            tied = [value for value in ranked if abs(value[0] - best_value) < 1e-12]
            runner_up = ranked[1][0] if len(ranked) > 1 else float("nan")
            unique = len(tied) == 1
            best = ranked[0]
            key = (query["unit_id"], target_locus_id)
            directed_best[key] = best[2]["unit_id"] if unique else None
            best_hit_rows.append({
                "query_unit_id": query["unit_id"],
                "query_label": query["label"],
                "query_group": query["analysis_subgroup"],
                "query_ordinal": query["unit_ordinal_distal_to_proximal"],
                "target_candidate_id": target_locus_id,
                "target_label": locus_by_id[target_locus_id]["label"],
                "target_group": locus_by_id[target_locus_id]["analysis_subgroup"],
                "best_target_unit_id": best[2]["unit_id"],
                "best_target_ordinal": best[2]["unit_ordinal_distal_to_proximal"],
                "best_core_to_core_identity": best_value,
                "runner_up_core_to_core_identity": runner_up,
                "best_runner_up_margin": best_value - runner_up if len(ranked) > 1 else float("nan"),
                "best_is_unique": unique,
                "best_tie_count": len(tied),
                "same_ordinal": str(query["unit_ordinal_distal_to_proximal"]) == str(best[2]["unit_ordinal_distal_to_proximal"]),
            })
    for row in best_hit_rows:
        query_unit = next(unit for unit in units if unit["unit_id"] == row["query_unit_id"])
        target_unit = next(unit for unit in units if unit["unit_id"] == row["best_target_unit_id"])
        reverse = directed_best.get((target_unit["unit_id"], query_unit["candidate_id"]))
        row["reciprocal_unique_best"] = bool(row["best_is_unique"] and reverse == query_unit["unit_id"])
        query_id = query_unit["candidate_id"]
        target_id = target_unit["candidate_id"]
        query_order = int(locus_by_id[query_id]["locus_order"])
        target_order = int(locus_by_id[target_id]["locus_order"])
        if query_order < target_order:
            evidence = pair_cache[(query_id, target_id)]
            boundary_map = evidence["target_to_query"]
        else:
            evidence = pair_cache[(target_id, query_id)]
            boundary_map = evidence["query_to_target"]
        query_start_rel = locus_relative_start[query_id]
        target_start_rel = locus_relative_start[target_id]
        query_left_index = int(query_unit["left_copy_start_relative_to_RT"]) - query_start_rel
        query_right_index = int(query_unit["right_copy_start_relative_to_RT"]) - query_start_rel
        mapped_left_index = boundary_map.get(query_left_index)
        mapped_right_index = boundary_map.get(query_right_index)
        mapped_left_rel = target_start_rel + mapped_left_index if mapped_left_index is not None else None
        mapped_right_rel = target_start_rel + mapped_right_index if mapped_right_index is not None else None
        left_offset = mapped_left_rel - int(target_unit["left_copy_start_relative_to_RT"]) if mapped_left_rel is not None else None
        right_offset = mapped_right_rel - int(target_unit["right_copy_start_relative_to_RT"]) if mapped_right_rel is not None else None
        coordinate_supported = left_offset is not None and right_offset is not None and abs(left_offset) <= 10 and abs(right_offset) <= 10
        row["mapped_query_left_boundary_in_target_relative"] = mapped_left_rel
        row["mapped_query_right_boundary_in_target_relative"] = mapped_right_rel
        row["left_boundary_offset_to_best_target_nt"] = left_offset
        row["right_boundary_offset_to_best_target_nt"] = right_offset
        row["whole_locus_boundary_support_within_10nt"] = coordinate_supported
        row["sequence_and_coordinate_supported_correspondence"] = bool(
            row["best_is_unique"] and row["reciprocal_unique_best"] and coordinate_supported
        )
        if row["sequence_and_coordinate_supported_correspondence"]:
            row["correspondence_status"] = "unique_reciprocal_sequence_best_plus_whole_locus_boundaries"
        elif row["reciprocal_unique_best"]:
            row["correspondence_status"] = "unique_reciprocal_sequence_best_without_whole_locus_boundary_support"
        else:
            row["correspondence_status"] = "free_sequence_best_only"

    unit_by_id = {unit["unit_id"]: unit for unit in units}
    supported_pair_rows = []
    seen_supported_pairs = set()
    for row in best_hit_rows:
        if not row["sequence_and_coordinate_supported_correspondence"]:
            continue
        pair_key = tuple(sorted((row["query_unit_id"], row["best_target_unit_id"])))
        if pair_key in seen_supported_pairs:
            continue
        seen_supported_pairs.add(pair_key)
        left = unit_by_id[pair_key[0]]
        right = unit_by_id[pair_key[1]]
        pair = pair_lookup[pair_key]
        supported_pair_rows.append({
            "relationship": "within_group" if left["analysis_subgroup"] == right["analysis_subgroup"] else "cross_group",
            "left_unit_id": left["unit_id"],
            "left_label": left["label"],
            "left_group": left["analysis_subgroup"],
            "left_ordinal": left["unit_ordinal_distal_to_proximal"],
            "right_unit_id": right["unit_id"],
            "right_label": right["label"],
            "right_group": right["analysis_subgroup"],
            "right_ordinal": right["unit_ordinal_distal_to_proximal"],
            "core_to_core_global_identity": pair["core_to_core_global_identity"],
            "trimmed_spacer_global_identity": pair["trimmed_spacer_global_identity"],
            "evidence_rule": "unique reciprocal free sequence best plus both whole-locus boundaries within 10 nt",
        })

    locus_correspondence_rows = []
    locus_pairs = {}
    for left_index, left in enumerate(loci):
        for right in loci[left_index + 1:]:
            locus_pairs[(left["candidate_id"], right["candidate_id"])] = (left, right)
    for (left_id, right_id), (left, right) in locus_pairs.items():
        pairs = [
            row for row in supported_pair_rows
            if {unit_by_id[row["left_unit_id"]]["candidate_id"], unit_by_id[row["right_unit_id"]]["candidate_id"]} == {left_id, right_id}
        ]
        ordinal_pairs = []
        for pair in pairs:
            left_unit = unit_by_id[pair["left_unit_id"]]
            right_unit = unit_by_id[pair["right_unit_id"]]
            if left_unit["candidate_id"] == left_id:
                ordinal_pairs.append("{}:{}".format(left_unit["unit_ordinal_distal_to_proximal"], right_unit["unit_ordinal_distal_to_proximal"]))
            else:
                ordinal_pairs.append("{}:{}".format(right_unit["unit_ordinal_distal_to_proximal"], left_unit["unit_ordinal_distal_to_proximal"]))
        locus_correspondence_rows.append({
            "left_candidate_id": left_id,
            "left_label": left["label"],
            "left_group": left["analysis_subgroup"],
            "right_candidate_id": right_id,
            "right_label": right["label"],
            "right_group": right["analysis_subgroup"],
            "relationship": "within_group" if left["analysis_subgroup"] == right["analysis_subgroup"] else "cross_group",
            "supported_unit_pair_count": len(pairs),
            "supported_ordinal_pairs_left_to_right": sorted(ordinal_pairs),
            "core_identity_median": statistics.median(float(pair["core_to_core_global_identity"]) for pair in pairs) if pairs else None,
            "trimmed_spacer_identity_median": statistics.median(float(pair["trimmed_spacer_global_identity"]) for pair in pairs) if pairs else None,
        })

    # LOO group check: seed is reselected in the remaining loci; held-out starts
    # are evaluated against mapped remaining-locus slots.  This is a consistency
    # check, not an independent detector validation.
    loo_rows = []
    for group_name, group_loci in groups.items():
        for heldout in group_loci:
            remaining = [locus for locus in group_loci if locus["candidate_id"] != heldout["candidate_id"]]
            loo_seed, _, _ = choose_group_seed(remaining, interval_sequences)
            loo_medoid = medoids[group_name] if medoids[group_name]["candidate_id"] != heldout["candidate_id"] else remaining[0]
            medoid_id = loo_medoid["candidate_id"]
            held_evidence = alignment_evidence(aligner, locus_sequences[medoid_id], locus_sequences[heldout["candidate_id"]])
            held_mapping = held_evidence["query_to_target"]
            held_occurrences = seed_occurrences(interval_sequences[heldout["candidate_id"]], loo_seed, max_mismatch=2)
            held_rel_start = -len(interval_sequences[heldout["candidate_id"]])
            held_locus_start = locus_relative_start[heldout["candidate_id"]]
            medoid_locus_start = locus_relative_start[medoid_id]
            held_mapped = []
            for occurrence in held_occurrences:
                relative = held_rel_start + occurrence["start"]
                query_index = relative - held_locus_start
                mapped_index = held_mapping.get(query_index)
                if mapped_index is not None:
                    held_mapped.append((medoid_locus_start + mapped_index, relative, occurrence))
            expected_slots = []
            heldout_ordinals = sorted(
                int(row["copy_ordinal_distal_to_proximal"])
                for row in copy_rows if row["candidate_id"] == heldout["candidate_id"]
            )
            for ordinal in heldout_ordinals:
                mapped_values = []
                for locus in remaining:
                    copy = next((row for row in copy_rows if row["candidate_id"] == locus["candidate_id"] and int(row["copy_ordinal_distal_to_proximal"]) == ordinal), None)
                    if copy is None:
                        continue
                    relative = int(copy["calibrated_copy_start_relative_to_RT"])
                    if locus["candidate_id"] == medoid_id:
                        mapped_values.append(relative)
                    else:
                        evidence = alignment_evidence(aligner, locus_sequences[medoid_id], locus_sequences[locus["candidate_id"]])
                        query_index = relative - locus_relative_start[locus["candidate_id"]]
                        mapped_index = evidence["query_to_target"].get(query_index)
                        if mapped_index is not None:
                            mapped_values.append(medoid_locus_start + mapped_index)
                if mapped_values:
                    expected_slots.append((ordinal, statistics.median(mapped_values)))
            for ordinal, expected in expected_slots:
                nearest = min(held_mapped, key=lambda value: abs(value[0] - expected)) if held_mapped else None
                recovered = nearest is not None and abs(nearest[0] - expected) <= 10
                loo_rows.append({
                    "analysis_subgroup": group_name,
                    "heldout_candidate_id": heldout["candidate_id"],
                    "heldout_label": heldout["label"],
                    "loo_group_seed_word": loo_seed,
                    "loo_medoid_label": loo_medoid["label"],
                    "copy_ordinal_distal_to_proximal": ordinal,
                    "expected_medoid_relative_start": expected,
                    "nearest_heldout_medoid_relative_start": nearest[0] if nearest else None,
                    "nearest_heldout_native_relative_start": nearest[1] if nearest else None,
                    "mapped_offset_nt": abs(nearest[0] - expected) if nearest else None,
                    "recovered_within_10nt": recovered,
                    "interpretation": "leave-one-locus-out within-group consistency check; not independent validation",
                })

    # Output tables.
    write_tsv(os.path.join(output_dir, "locus_pairwise_identity.tsv"), locus_pair_rows, list(locus_pair_rows[0]))
    write_tsv(os.path.join(output_dir, "medoid_alignment_summary.tsv"), alignment_summary_rows, list(alignment_summary_rows[0]))
    write_tsv(os.path.join(output_dir, "whole_locus_alignment_blocks.tsv"), alignment_block_rows, list(alignment_block_rows[0]))
    with open(os.path.join(output_dir, "whole_locus_medoid_alignments.txt"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(alignment_text))
    write_tsv(os.path.join(output_dir, "group_seed_ranking.tsv"), seed_ranking_rows, list(seed_ranking_rows[0]))
    write_tsv(os.path.join(output_dir, "seed_view_phase_audit.tsv"), seed_view_rows, list(seed_view_rows[0]))
    write_tsv(os.path.join(output_dir, "canonical_seed_occurrence_audit.tsv"), canonical_occurrence_audit_rows, list(canonical_occurrence_audit_rows[0]))
    write_tsv(os.path.join(output_dir, "unified_array_annotations.tsv"), annotation_rows, list(annotation_rows[0]))
    write_tsv(os.path.join(output_dir, "unified_repeat_copies.tsv"), copy_rows, list(copy_rows[0]))
    boundary_columns = [
        "analysis_subgroup", "candidate_id", "label", "side", "alternative_type",
        "expected_start_relative_to_RT", "candidate_start_relative_to_RT", "offset_from_expected_nt",
        "profile_match_26", "seed_mismatches", "candidate_fixed26_sequence",
        "group_edge_candidate_loci", "group_majority_threshold", "status", "promoted_copy_id", "reason",
    ]
    write_tsv(os.path.join(output_dir, "boundary_alternatives.tsv"), boundary_alternative_rows, boundary_columns)
    write_tsv(os.path.join(output_dir, "edge_scan_audit.tsv"), edge_scan_rows, list(edge_scan_rows[0]))
    write_tsv(os.path.join(output_dir, "unified_units.tsv"), units, list(units[0]))
    write_fasta(os.path.join(output_dir, "unified_unit_sequences.fna"), unit_fastas)
    write_tsv(os.path.join(output_dir, "unit_pairwise_identity.tsv"), unit_pair_rows, list(unit_pair_rows[0]))
    write_tsv(os.path.join(output_dir, "unit_free_best_hits.tsv"), best_hit_rows, list(best_hit_rows[0]))
    write_tsv(os.path.join(output_dir, "supported_unit_correspondences.tsv"), supported_pair_rows, list(supported_pair_rows[0]))
    write_tsv(os.path.join(output_dir, "locus_pair_unit_correspondence_summary.tsv"), locus_correspondence_rows, list(locus_correspondence_rows[0]))
    write_tsv(os.path.join(output_dir, "cross_group_medoid_copy_mapping.tsv"), cross_copy_mapping_rows, list(cross_copy_mapping_rows[0]))
    write_tsv(os.path.join(output_dir, "cross_group_medoid_unit_boundary_mapping.tsv"), cross_unit_boundary_rows, list(cross_unit_boundary_rows[0]))
    write_tsv(os.path.join(output_dir, "cross_group_medoid_array_alignment_blocks.tsv"), cross_medoid_block_rows, list(cross_medoid_block_rows[0]))
    write_tsv(os.path.join(output_dir, "leave_one_locus_out_recovery.tsv"), loo_rows, list(loo_rows[0]))

    summary = {
        "analysis": "P2 unified ten-locus ART array annotation and comparison",
        "fixed_loci": len(loci),
        "groups": {},
        "total_copies": len(copy_rows),
        "total_units": len(units),
        "seed_chain_copies": sum(row["copy_class"] != "group_supported_degenerate_edge" for row in copy_rows),
        "seed_chain_units": sum(int(row["seed_chain_copy_count"]) - 1 for row in annotation_rows),
        "medoids": {group: medoid["label"] for group, medoid in medoids.items()},
        "cross_group_medoid_anchor_through_RT_global_identity": cross_medoid_evidence["identity"],
        "LOO_rows": len(loo_rows),
        "LOO_recovered_within_10nt": sum(bool(row["recovered_within_10nt"]) for row in loo_rows),
        "parameters": "scripts/p2_analysis/ANALYSIS_PARAMETERS.md",
        "biopython_runtime_receipt": "tools/biopython_runtime_receipt.json",
        "input_build_receipt": "data/source_metadata/p2/fixed_10_locus_build_receipt.json",
        "unresolved_boundary_alternatives": sum(
            row.get("status") != "promoted_to_main_group_supported_degenerate_edge"
            for row in boundary_alternative_rows
        ),
        "promoted_group_supported_degenerate_edge_copies": len(promoted_edge_rows),
        "within_group_directed_free_best_rows": sum(row["query_group"] == row["target_group"] for row in best_hit_rows),
        "within_group_sequence_and_coordinate_supported_rows": sum(
            row["query_group"] == row["target_group"] and row["sequence_and_coordinate_supported_correspondence"]
            for row in best_hit_rows
        ),
        "cross_group_directed_free_best_rows": sum(row["query_group"] != row["target_group"] for row in best_hit_rows),
        "cross_group_sequence_and_coordinate_supported_rows": sum(
            row["query_group"] != row["target_group"] and row["sequence_and_coordinate_supported_correspondence"]
            for row in best_hit_rows
        ),
        "unique_supported_unit_pairs": {
            "P0_seed7_within": sum(
                row["relationship"] == "within_group" and row["left_group"] == "P0_seed7"
                for row in supported_pair_rows
            ),
            "PB50_near3_within": sum(
                row["relationship"] == "within_group" and row["left_group"] == "PB50_near3"
                for row in supported_pair_rows
            ),
            "cross_group": sum(row["relationship"] == "cross_group" for row in supported_pair_rows),
        },
    }
    for group_name, group_loci in groups.items():
        group_annotations = [row for row in annotation_rows if row["analysis_subgroup"] == group_name]
        summary["groups"][group_name] = {
            "loci": len(group_loci),
            "selected_seed_word": group_seed[group_name],
            "operational_copy_counts": {row["label"]: row["copy_count"] for row in group_annotations},
            "seed_chain_copy_counts": {row["label"]: row["seed_chain_copy_count"] for row in group_annotations},
            "promoted_low_confidence_edge_counts": {row["label"]: row["promoted_degenerate_edge_copy_count"] for row in group_annotations},
            "operational_unit_counts": {row["label"]: int(row["copy_count"]) - 1 for row in group_annotations},
            "repeat_block_consensus": group_annotations[0]["repeat_block_consensus"],
            "repeat_block_offsets": [group_annotations[0]["repeat_block_left_offset_from_seed"], group_annotations[0]["repeat_block_right_offset_exclusive_from_seed"]],
        }
    write_json(os.path.join(output_dir, "p2_comparison_summary.json"), summary)
    receipt = {
        "status": "pass",
        "analysis_date": "2026-09-24",
        "command": "project-local Python 3.12 scripts/p2_analysis/run_unified_locus_analysis.py",
        "input_locus_count": len(loci),
        "input_build_receipt": "data/source_metadata/p2/fixed_10_locus_build_receipt.json",
        "parameter_record": "scripts/p2_analysis/ANALYSIS_PARAMETERS.md",
        "implementation_correction": {
            "issue": "The first provisional aggregation left four P0 proximal edge hits in the alternative table even though they met the <=2-mismatch route and the 4/7 group-majority slot rule predeclared in the P2 parameter record.",
            "resolution": "Before final delivery, the four -61/-63 hits were promoted to low-confidence operational edge candidates and every affected unit, correspondence, LOO, summary, and figure table was regenerated.",
            "scientific_scope": "They remain low-confidence candidate copies (14/26 profile matches, two seed mismatches), not experimentally established biological copies; the core correspondence claim does not depend on them.",
        },
        "validation": {
            "fixed_loci": len(loci),
            "seed_chain_copies": summary["seed_chain_copies"],
            "promoted_low_confidence_edges": summary["promoted_group_supported_degenerate_edge_copies"],
            "operational_copies": summary["total_copies"],
            "operational_units": summary["total_units"],
            "LOO_recovered": summary["LOO_recovered_within_10nt"],
            "LOO_rows": summary["LOO_rows"],
            "unresolved_boundary_alternatives": summary["unresolved_boundary_alternatives"],
        },
        "output_files": sorted(
            name for name in os.listdir(output_dir)
            if os.path.isfile(os.path.join(output_dir, name)) and name != "run_receipt.json"
        ),
    }
    write_json(os.path.join(output_dir, "run_receipt.json"), receipt)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
