#!/usr/bin/env python3
"""Render the three P2 research figures from completed comparison tables."""

from __future__ import division

import argparse
import copy as copy_module
import csv
import math
import os
from collections import Counter, defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D


P0_COLOR = "#3B6FB6"
NEAR3_COLOR = "#D97835"
EDGE_COLOR = "#8E5BA6"
ALT_COLOR = "#B94A48"
TEXT = "#252525"
LIGHT = "#E8E8E8"


def configure():
    matplotlib.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.titlesize": 9.5,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
    })


def read_tsv(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save_figure(fig, figure_dir, stem):
    os.makedirs(figure_dir, exist_ok=True)
    paths = []
    for extension in ("pdf", "svg", "png"):
        path = os.path.join(figure_dir, stem + "." + extension)
        kwargs = {"facecolor": "white"}
        if extension == "png":
            kwargs["dpi"] = 300
        fig.savefig(path, **kwargs)
        paths.append(path)
    plt.close(fig)
    return paths


def plot_architecture(root, out_dir, figure_dir):
    loci = read_tsv(os.path.join(root, "data", "processed", "p2_loci", "fixed_10_loci.tsv"))
    annotations = read_tsv(os.path.join(out_dir, "unified_array_annotations.tsv"))
    copies = read_tsv(os.path.join(out_dir, "unified_repeat_copies.tsv"))
    units = read_tsv(os.path.join(out_dir, "unified_units.tsv"))
    alternatives = read_tsv(os.path.join(out_dir, "boundary_alternatives.tsv"))
    coding = read_tsv(os.path.join(root, "data", "processed", "p2_loci", "coding_context.tsv"))
    annotation_by_id = {row["candidate_id"]: row for row in annotations}
    copies_by_id = defaultdict(list)
    units_by_id = defaultdict(list)
    coding_by_id = defaultdict(list)
    for row in copies:
        copies_by_id[row["candidate_id"]].append(row)
    for row in units:
        units_by_id[row["candidate_id"]].append(row)
    for row in coding:
        coding_by_id[row["candidate_id"]].append(row)

    fig = plt.figure(figsize=(175 / 25.4, 132 / 25.4))
    grid = GridSpec(1, 2, width_ratios=[1.82, 1.0], left=0.105, right=0.92, bottom=0.19, top=0.93, wspace=0.23)
    ax = fig.add_subplot(grid[0, 0])
    heat = fig.add_subplot(grid[0, 1])

    labels = [row["label"] for row in loci]
    y_positions = {row["candidate_id"]: len(loci) - 1 - index for index, row in enumerate(loci)}
    ax.axhspan(2.5, 9.5, color=P0_COLOR, alpha=0.045, zorder=0)
    ax.axhspan(-0.5, 2.5, color=NEAR3_COLOR, alpha=0.055, zorder=0)
    for locus in loci:
        candidate_id = locus["candidate_id"]
        y = y_positions[candidate_id]
        group_color = P0_COLOR if locus["analysis_subgroup"] == "P0_seed7" else NEAR3_COLOR
        anchor_start = int(locus["distal_anchor_relative_start_0based"])
        anchor_end = int(locus["distal_anchor_relative_end_0based"])
        ax.plot([anchor_end, 0], [y, y], color="#A9A9A9", lw=1.0, zorder=1)
        ax.add_patch(patches.Rectangle((anchor_start, y - 0.18), anchor_end - anchor_start + 1, 0.36,
                                       facecolor="#B8B8B8", edgecolor="none", zorder=2))
        ax.add_patch(patches.Rectangle((0, y - 0.24), 70, 0.48,
                                       facecolor="#4A4A4A", edgecolor="none", zorder=2))
        for feature in coding_by_id[candidate_id]:
            if feature["feature_role"] != "intervening_array_annotation_conflict":
                continue
            start = int(feature["oriented_relative_start_0based"])
            end = int(feature["oriented_relative_end_0based"])
            ax.add_patch(patches.Rectangle((start, y - 0.27), end - start + 1, 0.54,
                                           facecolor="#D65F5F", alpha=0.20, edgecolor="#B94A48",
                                           linewidth=0.5, hatch="///", zorder=1))
        locus_copies = sorted(copies_by_id[candidate_id], key=lambda row: int(row["calibrated_copy_start_relative_to_RT"]))
        positions = [int(row["calibrated_copy_start_relative_to_RT"]) for row in locus_copies]
        if positions:
            ax.plot([positions[0], positions[-1]], [y, y], color=group_color, lw=1.5, zorder=3)
        for copy in locus_copies:
            x = int(copy["calibrated_copy_start_relative_to_RT"])
            if copy["copy_class"] == "group_supported_degenerate_edge":
                ax.scatter([x], [y], marker="^", s=38, facecolor="white", edgecolor=EDGE_COLOR,
                           linewidth=1.1, zorder=5)
            elif copy["seed_mismatches"] == "1":
                ax.scatter([x], [y], marker="D", s=22, facecolor=group_color, edgecolor="white",
                           linewidth=0.4, zorder=5)
            else:
                ax.scatter([x], [y], marker="o", s=24, facecolor=group_color, edgecolor="white",
                           linewidth=0.4, zorder=5)
        for alternative in alternatives:
            if alternative["candidate_id"] != candidate_id or alternative["status"].startswith("promoted"):
                continue
            x = int(alternative["candidate_start_relative_to_RT"])
            ax.scatter([x], [y], marker="x", s=32, color=ALT_COLOR, linewidth=1.2, zorder=6)

    ax.axvline(0, color="#333333", lw=0.7)
    ax.set_xlim(-1660, 270)
    ax.set_ylim(-0.65, 9.65)
    ax.set_yticks(range(10))
    ax.set_yticklabels(list(reversed(labels)))
    ax.set_xlabel("Position relative to RT 5′ base (nt)")
    ax.set_title("a  Unified locus annotation", loc="left", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    legend = [
        patches.Patch(facecolor="#B8B8B8", label="distal coding anchor"),
        patches.Patch(facecolor="#4A4A4A", label="RT (start shown)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=P0_COLOR, markeredgecolor="white", markersize=6, label="periodic-chain copy"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=P0_COLOR, markeredgecolor="white", markersize=5, label="1-mismatch seed copy"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor="white", markeredgecolor=EDGE_COLOR, markersize=6, label="low-confidence edge candidate"),
        Line2D([0], [0], marker="x", color=ALT_COLOR, linestyle="none", markersize=6, label="unresolved register"),
    ]
    fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.51, 0.01), fontsize=6.2, handletextpad=0.5, columnspacing=1.0)

    max_unit = 5
    matrix = np.full((len(loci), max_unit), np.nan)
    edge_units = set()
    for row_index, locus in enumerate(loci):
        for unit in units_by_id[locus["candidate_id"]]:
            column = int(unit["unit_ordinal_distal_to_proximal"]) - 1
            matrix[row_index, column] = float(unit["core_to_core_length_nt"])
            right_copy = next(row for row in copies_by_id[locus["candidate_id"]] if row["copy_id"] == unit["right_copy_id"])
            if right_copy["copy_class"] == "group_supported_degenerate_edge":
                edge_units.add((row_index, column))
    cmap = copy_module.copy(plt.cm.viridis)
    cmap.set_bad("#EEEEEE")
    image = heat.imshow(matrix, aspect="auto", cmap=cmap, vmin=165, vmax=245, interpolation="none")
    for row_index in range(len(loci)):
        for column in range(max_unit):
            if np.isnan(matrix[row_index, column]):
                continue
            heat.text(column, row_index, str(int(matrix[row_index, column])), ha="center", va="center",
                      fontsize=6.8, color="white" if matrix[row_index, column] < 195 else "#1A1A1A")
            if (row_index, column) in edge_units:
                heat.add_patch(patches.Rectangle((column - 0.47, row_index - 0.47), 0.94, 0.94,
                                                 fill=False, edgecolor=EDGE_COLOR, linewidth=1.4))
    heat.axhline(6.5, color="white", lw=2.0)
    heat.set_xticks(range(max_unit))
    heat.set_xticklabels(["U{}".format(index) for index in range(1, max_unit + 1)])
    heat.set_yticks(range(len(loci)))
    heat.set_yticklabels(labels)
    heat.set_xlabel("Core-to-core unit")
    heat.set_title("b  Core-to-core lengths (nt)", loc="left", fontweight="bold")
    for spine in heat.spines.values():
        spine.set_visible(False)
    colorbar = fig.colorbar(image, ax=heat, fraction=0.046, pad=0.035)
    colorbar.set_label("Length (nt)")
    return save_figure(fig, figure_dir, "figure1_unified_array_architecture")


def plot_correspondence(root, out_dir, figure_dir):
    loci = read_tsv(os.path.join(root, "data", "processed", "p2_loci", "fixed_10_loci.tsv"))
    pairwise = read_tsv(os.path.join(out_dir, "locus_pairwise_identity.tsv"))
    locus_corr = read_tsv(os.path.join(out_dir, "locus_pair_unit_correspondence_summary.tsv"))
    supported = read_tsv(os.path.join(out_dir, "supported_unit_correspondences.tsv"))
    labels = [row["label"] for row in loci]
    index = {label: idx for idx, label in enumerate(labels)}
    matrix = np.eye(len(labels))
    for row in pairwise:
        left = index[row["left_label"]]
        right = index[row["right_label"]]
        value = float(row["intergenic_global_identity"])
        matrix[left, right] = value
        matrix[right, left] = value

    fig = plt.figure(figsize=(175 / 25.4, 142 / 25.4))
    grid = GridSpec(2, 2, height_ratios=[1.0, 0.88], width_ratios=[0.94, 1.06],
                    left=0.145, right=0.96, bottom=0.105, top=0.94, hspace=0.36, wspace=0.46)
    ax1 = fig.add_subplot(grid[:, 0])
    ax2 = fig.add_subplot(grid[0, 1])
    ax3 = fig.add_subplot(grid[1, 1])

    image = ax1.imshow(matrix, cmap="Blues", vmin=0.5, vmax=1.0, interpolation="none")
    ax1.axhline(6.5, color="white", lw=2)
    ax1.axvline(6.5, color="white", lw=2)
    ax1.set_xticks(range(len(labels)))
    ax1.set_xticklabels(labels, rotation=55, ha="right")
    ax1.set_yticks(range(len(labels)))
    ax1.set_yticklabels(labels)
    ax1.set_title("a  Anchor-to-RT identity", loc="left", fontweight="bold")
    for spine in ax1.spines.values():
        spine.set_visible(False)
    colorbar = fig.colorbar(image, ax=ax1, fraction=0.046, pad=0.03)
    colorbar.set_label("Identity")

    p0 = labels[:7]
    near3 = labels[7:]
    cross_matrix = np.zeros((len(p0), len(near3)))
    text_matrix = [["" for _ in near3] for _ in p0]
    for row in locus_corr:
        if row["relationship"] != "cross_group":
            continue
        left_label = row["left_label"]
        right_label = row["right_label"]
        if left_label in p0:
            p0_label, near_label = left_label, right_label
            ordinal_pairs = row["supported_ordinal_pairs_left_to_right"]
        else:
            p0_label, near_label = right_label, left_label
            ordinal_pairs = row["supported_ordinal_pairs_left_to_right"]
        i = p0.index(p0_label)
        j = near3.index(near_label)
        count = int(row["supported_unit_pair_count"])
        cross_matrix[i, j] = count
        ordinals = []
        if ordinal_pairs:
            for pair in ordinal_pairs.split(";"):
                values = pair.split(":")
                if len(values) == 2 and values[0] == values[1]:
                    ordinals.append(values[0])
                else:
                    ordinals.append(pair.replace(":", "→"))
        text_matrix[i][j] = ",".join(ordinals)
    cross_cmap = ListedColormap(["#F0F0F0", "#C9D8EC", "#79A4D6", "#315F9A"])
    ax2.imshow(cross_matrix, cmap=cross_cmap, vmin=0, vmax=3, interpolation="none", aspect="auto")
    for i in range(len(p0)):
        for j in range(len(near3)):
            text = text_matrix[i][j] if text_matrix[i][j] else "–"
            ax2.text(j, i, text, ha="center", va="center", fontsize=7,
                     color="white" if cross_matrix[i, j] >= 3 else TEXT)
    ax2.set_xticks(range(len(near3)))
    ax2.set_xticklabels(near3)
    ax2.set_yticks(range(len(p0)))
    ax2.set_yticklabels(p0)
    ax2.set_xlabel("PB50-near group locus")
    ax2.set_title("b  Cross-group unit correspondence", loc="left", fontweight="bold", fontsize=8.6)
    for spine in ax2.spines.values():
        spine.set_visible(False)

    categories = [
        ("P0 within\n(n=83)", [row for row in supported if row["relationship"] == "within_group" and row["left_group"] == "P0_seed7"], P0_COLOR),
        ("PB50 group\n(n=9)", [row for row in supported if row["relationship"] == "within_group" and row["left_group"] == "PB50_near3"], NEAR3_COLOR),
        ("Cross-group\n(n=35)", [row for row in supported if row["relationship"] == "cross_group"], "#4F8F72"),
    ]
    for category_index, (label, rows, color) in enumerate(categories):
        core_values = [float(row["core_to_core_global_identity"]) for row in rows]
        spacer_values = [float(row["trimmed_spacer_global_identity"]) for row in rows]
        for point_index, (core, spacer) in enumerate(zip(core_values, spacer_values)):
            jitter = 0.075 * math.sin((point_index + 1) * 2.399)
            ax3.plot([category_index - 0.09 + jitter, category_index + 0.09 + jitter], [core, spacer],
                     color=color, alpha=0.16, lw=0.45, zorder=1)
            ax3.scatter([category_index - 0.09 + jitter], [core], s=10, color=color, alpha=0.50, edgecolor="none", zorder=2)
            ax3.scatter([category_index + 0.09 + jitter], [spacer], s=10, facecolor="white", edgecolor=color,
                        linewidth=0.55, alpha=0.65, zorder=2)
        for shift, values in ((-0.09, core_values), (0.09, spacer_values)):
            median = float(np.median(values))
            ax3.plot([category_index + shift - 0.09, category_index + shift + 0.09], [median, median],
                     color="#111111", lw=1.2, zorder=3)
    ax3.set_xticks(range(3))
    ax3.set_xticklabels([item[0] for item in categories])
    ax3.set_ylim(0.45, 1.01)
    ax3.set_ylabel("Global nucleotide identity")
    ax3.set_title("c  Supported-pair identities", loc="left", fontweight="bold", fontsize=8.6)
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    ax3.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#666666", markersize=5, label="core-to-core unit"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="#666666", markersize=5, label="trimmed spacer"),
    ], frameon=False, loc="lower left", fontsize=6.6)
    return save_figure(fig, figure_dir, "figure2_locus_and_unit_correspondence")


def plot_local_differences(root, out_dir, figure_dir):
    copies = read_tsv(os.path.join(out_dir, "unified_repeat_copies.tsv"))
    edge_scan = read_tsv(os.path.join(out_dir, "edge_scan_audit.tsv"))
    copy_map = read_tsv(os.path.join(out_dir, "cross_group_medoid_copy_mapping.tsv"))
    blocks = read_tsv(os.path.join(out_dir, "cross_group_medoid_array_alignment_blocks.tsv"))

    fig = plt.figure(figsize=(175 / 25.4, 148 / 25.4))
    grid = GridSpec(2, 2, height_ratios=[1.12, 0.88], width_ratios=[1.15, 0.85],
                    left=0.19, right=0.96, bottom=0.105, top=0.94, hspace=0.37, wspace=0.34)
    ax1 = fig.add_subplot(grid[0, :])
    ax2 = fig.add_subplot(grid[1, 0])
    ax3 = fig.add_subplot(grid[1, 1])

    near_copies = [row for row in copies if row["analysis_subgroup"] == "PB50_near3"]
    order = {"AH12": 0, "Machias": 1, "PB50": 2}
    near_copies.sort(key=lambda row: (order[row["label"]], int(row["copy_ordinal_distal_to_proximal"])))
    sequences = [row["fixed26_sequence"] for row in near_copies]
    consensus = "".join(
        sorted(Counter(sequence[column] for sequence in sequences), key=lambda base: (-Counter(sequence[column] for sequence in sequences)[base], base))[0]
        for column in range(26)
    )
    mismatch = np.array([[int(base != consensus[column]) for column, base in enumerate(sequence)] for sequence in sequences])
    ax1.imshow(mismatch, aspect="auto", cmap=ListedColormap(["#F5E8D8", "#B6493A"]), vmin=0, vmax=1, interpolation="none")
    for row_index, sequence in enumerate(sequences):
        for column, base in enumerate(sequence):
            ax1.text(column, row_index, base, ha="center", va="center", family="monospace", fontsize=6.5,
                     color="white" if mismatch[row_index, column] else TEXT)
    ax1.add_patch(patches.Rectangle((5.5, -0.48), 10, len(sequences) - 0.04, fill=False, edgecolor="#2F2F2F", linewidth=1.0))
    ylabels = []
    for row in near_copies:
        short_label = {"AH12": "AH", "Machias": "MA", "PB50": "PB*"}[row["label"]]
        ylabels.append("{} C{} | {} | d14={}".format(
            short_label, row["copy_ordinal_distal_to_proximal"], row["calibrated_copy_start_relative_to_RT"],
            row["MarsHill14_best_mismatches_in_window61"],
        ))
    ax1.set_yticks(range(len(near_copies)))
    ax1.set_yticklabels(ylabels, fontsize=6.7)
    ax1.set_xticks(range(0, 26, 2))
    ax1.set_xticklabels([str(index - 6) for index in range(0, 26, 2)])
    ax1.set_xlabel("Position relative to selected 10-nt seed start (boxed: seed)")
    ax1.set_title("a  PB50-group copies retained by periodic context", loc="left", fontweight="bold")
    for spine in ax1.spines.values():
        spine.set_visible(False)

    p0_proximal = [row for row in edge_scan if row["analysis_subgroup"] == "P0_seed7" and row["side"] == "proximal"]
    label_order = ["SA1", "MarsHill", "Madawaska", "LY01", "S6", "PALS_2", "UFV_DC4"]
    p0_proximal.sort(key=lambda row: label_order.index(row["label"]))
    for row_index, row in enumerate(p0_proximal):
        accepted = row["meets_single_locus_candidate_floor"] == "true"
        score = float(row["profile_match_26"])
        marker = "^" if accepted else "o"
        color = EDGE_COLOR if accepted else "#9A9A9A"
        ax2.scatter([score], [row_index], marker=marker, s=44 if accepted else 26,
                    facecolor="white" if accepted else color, edgecolor=color, linewidth=1.0, zorder=3)
        ax2.text(score + 0.35, row_index, "{} mm; {} nt".format(row["seed_mismatches"], row["candidate_start_relative_to_RT"]),
                 va="center", fontsize=6.6, color=color)
    ax2.axvline(18, color="#555555", lw=0.8, linestyle="--")
    ax2.set_xlim(12.5, 24.5)
    ax2.set_ylim(len(p0_proximal) - 0.4, -0.8)
    ax2.set_yticks(range(len(p0_proximal)))
    ax2.set_yticklabels([row["label"] for row in p0_proximal])
    ax2.set_xlabel("Best 26-nt profile matches at expected proximal edge")
    ax2.set_title("b  One-pass P0 edge decision", loc="left", fontweight="bold")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.text(0.55, 0.89, "Triangles: ≤2-mismatch route\n4/7 same slot; all 14/26",
             transform=ax2.transAxes, fontsize=6.4, ha="left", va="top", color=TEXT)

    p0_native = [row for row in copy_map if row["source_label"] == "PALS_2"]
    pb_mapped = [row for row in copy_map if row["source_label"] == "PB50"]
    for y, label, color in ((1.0, "PALS_2 native", P0_COLOR), (0.0, "PB50 mapped into PALS_2", NEAR3_COLOR)):
        ax3.plot([-1160, -20], [y, y], color="#BDBDBD", lw=1.0)
        display_label = "PALS_2 native" if y == 1.0 else "PB50"
        ax3.text(-1175, y, display_label, ha="right", va="center", fontsize=6.8, color=color)
    p0_positions = {}
    for row in p0_native:
        position = int(row["source_copy_relative_start"])
        ordinal = int(row["source_copy_ordinal"])
        p0_positions[ordinal] = position
        is_edge = ordinal == max(int(item["source_copy_ordinal"]) for item in p0_native)
        ax3.scatter([position], [1.0], marker="^" if is_edge else "o", s=38 if is_edge else 28,
                    facecolor="white" if is_edge else P0_COLOR, edgecolor=EDGE_COLOR if is_edge else "white",
                    linewidth=1.0 if is_edge else 0.4, zorder=4)
        ax3.text(position, 1.17, "C{}".format(ordinal), ha="center", va="bottom", fontsize=6.4)
    for row in pb_mapped:
        if row["mapping_status"] != "aligned_base":
            continue
        position = int(row["mapped_target_relative_start"])
        ordinal = int(row["source_copy_ordinal"])
        ax3.scatter([position], [0.0], marker="s", s=28, facecolor=NEAR3_COLOR, edgecolor="white", linewidth=0.4, zorder=4)
        ax3.text(position, -0.17, "C{}".format(ordinal), ha="center", va="top", fontsize=6.4)
    for row in blocks:
        if row["block_type"] != "query_gap" or int(row["alignment_columns"]) < 10:
            continue
        left = int(row["target_relative_start"])
        right = int(row["target_relative_end_exclusive"])
        ax3.add_patch(patches.Rectangle((left, 0.37), right - left, 0.26, facecolor="#D65F5F", alpha=0.45, edgecolor="none"))
    gap_copy = next(row for row in p0_native if row["mapping_status"] == "maps_to_gap")
    gap_position = int(gap_copy["source_copy_relative_start"])
    ax3.scatter([gap_position], [1.0], marker="x", s=58, color=ALT_COLOR, linewidth=1.3, zorder=5)
    ax3.text(0.96, 0.53, "PB50 gap blocks\n14, 50, 61, 24 nt",
             transform=ax3.transAxes, ha="right", va="bottom", fontsize=6.2, color=ALT_COLOR)
    ax3.annotate("C4 maps to a PB50 gap", xy=(gap_position, 0.86), xytext=(-780, 1.47),
                 arrowprops=dict(arrowstyle="-", lw=0.7, color=ALT_COLOR), fontsize=6.6, color=ALT_COLOR)
    for ordinal in (1, 2, 3):
        pb = next((row for row in pb_mapped if int(row["source_copy_ordinal"]) == ordinal), None)
        if pb and pb["mapping_status"] == "aligned_base" and ordinal in p0_positions:
            ax3.plot([p0_positions[ordinal], int(pb["mapped_target_relative_start"])], [0.96, 0.04], color="#777777", lw=0.55, alpha=0.7)
    pb4 = next(row for row in pb_mapped if int(row["source_copy_ordinal"]) == 4)
    ax3.plot([p0_positions[5], int(pb4["mapped_target_relative_start"])], [0.96, 0.04], color="#777777", lw=0.55, alpha=0.7)
    ax3.set_xlim(-1210, 10)
    ax3.set_ylim(-0.42, 1.62)
    ax3.set_yticks([])
    ax3.set_xlabel("PALS_2 coordinate relative to RT (nt)")
    ax3.set_title("c  Medoid boundary map", loc="left", fontweight="bold")
    ax3.spines["left"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    ax3.spines["top"].set_visible(False)
    return save_figure(fig, figure_dir, "figure3_localized_annotation_differences")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    args = parser.parse_args()
    root = os.path.abspath(args.project_root)
    out_dir = os.path.join(root, "data", "processed", "p2_comparison")
    figure_dir = os.path.join(root, "figures", "p2")
    configure()
    outputs = []
    outputs.extend(plot_architecture(root, out_dir, figure_dir))
    outputs.extend(plot_correspondence(root, out_dir, figure_dir))
    outputs.extend(plot_local_differences(root, out_dir, figure_dir))
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
