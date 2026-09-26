#!/usr/bin/env python3
"""Render revised Figures 1-4 from the accepted P2/P3 tables."""

from __future__ import annotations

import copy
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

from common import (
    BLUE,
    GRAY,
    LIGHT_GRAY,
    MISSING,
    ORANGE,
    PURPLE,
    RED,
    ROOT,
    TEAL,
    TEXT,
    clean_axis,
    configure,
    panel_title,
    read_tsv,
    save_figure,
)


P2_DIR = ROOT / "data" / "processed" / "p2_comparison"
LOCI_DIR = ROOT / "data" / "processed" / "p2_loci"
P3_DIR = ROOT / "data" / "processed" / "p3_calibration"


def _median(values: list[float]) -> float:
    ordered = sorted(values)
    mid = len(ordered) // 2
    return ordered[mid] if len(ordered) % 2 else 0.5 * (ordered[mid - 1] + ordered[mid])


def _format_numeric_minus(value: object) -> str:
    """Use U+2212 for numeric display without changing ordinary hyphens."""
    return str(value).replace("-", "−")


def render_figure_1() -> list[Path]:
    loci = read_tsv(LOCI_DIR / "fixed_10_loci.tsv")
    annotations = read_tsv(P2_DIR / "unified_array_annotations.tsv")
    copies = read_tsv(P2_DIR / "unified_repeat_copies.tsv")
    units = read_tsv(P2_DIR / "unified_units.tsv")
    alternatives = read_tsv(P2_DIR / "boundary_alternatives.tsv")
    coding = read_tsv(LOCI_DIR / "coding_context.tsv")

    copies_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    units_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    coding_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in copies:
        copies_by_id[row["candidate_id"]].append(row)
    for row in units:
        units_by_id[row["candidate_id"]].append(row)
    for row in coding:
        coding_by_id[row["candidate_id"]].append(row)

    fig = plt.figure(figsize=(170 / 25.4, 190 / 25.4))
    grid = GridSpec(
        4,
        2,
        figure=fig,
        height_ratios=[1.48, 0.18, 0.72, 0.50],
        width_ratios=[1.72, 1.0],
        left=0.14,
        right=0.90,
        bottom=0.10,
        top=0.95,
        hspace=0.42,
        wspace=0.28,
    )
    ax = fig.add_subplot(grid[0, 0])
    heat = fig.add_subplot(grid[0, 1])
    legend_axis = fig.add_subplot(grid[1, :])
    anatomy = fig.add_subplot(grid[2, :])
    sa1 = fig.add_subplot(grid[3, :])

    labels = [row["label"] for row in loci]
    y_positions = {row["candidate_id"]: len(loci) - 1 - index for index, row in enumerate(loci)}
    ax.axhspan(2.5, 9.5, color=BLUE, alpha=0.045, zorder=0)
    ax.axhspan(-0.5, 2.5, color=ORANGE, alpha=0.055, zorder=0)
    for locus in loci:
        candidate_id = locus["candidate_id"]
        y = y_positions[candidate_id]
        group_color = BLUE if locus["analysis_subgroup"] == "P0_seed7" else ORANGE
        anchor_start = int(locus["distal_anchor_relative_start_0based"])
        anchor_end = int(locus["distal_anchor_relative_end_0based"])
        ax.plot([anchor_end, 0], [y, y], color="#A9A9A9", lw=0.9, zorder=1)
        ax.add_patch(
            patches.Rectangle(
                (anchor_start, y - 0.17), anchor_end - anchor_start + 1, 0.34,
                facecolor="#B8B8B8", edgecolor="none", zorder=2,
            )
        )
        ax.add_patch(
            patches.Rectangle((0, y - 0.22), 70, 0.44, facecolor="#4A4A4A", edgecolor="none", zorder=2)
        )
        for feature in coding_by_id[candidate_id]:
            if feature["feature_role"] == "intervening_array_annotation_conflict":
                start = int(feature["oriented_relative_start_0based"])
                end = int(feature["oriented_relative_end_0based"])
                ax.add_patch(
                    patches.Rectangle(
                        (start, y - 0.24), end - start + 1, 0.48,
                        facecolor=RED, alpha=0.16, edgecolor=RED, linewidth=0.45,
                        hatch="///", zorder=1,
                    )
                )
        locus_copies = sorted(copies_by_id[candidate_id], key=lambda row: int(row["calibrated_copy_start_relative_to_RT"]))
        positions = [int(row["calibrated_copy_start_relative_to_RT"]) for row in locus_copies]
        if positions:
            ax.plot([positions[0], positions[-1]], [y, y], color=group_color, lw=1.4, zorder=3)
        for copy_row in locus_copies:
            x = int(copy_row["calibrated_copy_start_relative_to_RT"])
            if copy_row["copy_class"] == "group_supported_degenerate_edge":
                ax.scatter([x], [y], marker="^", s=35, facecolor="white", edgecolor=PURPLE, linewidth=1.0, zorder=5)
            elif copy_row["seed_mismatches"] == "1":
                ax.scatter([x], [y], marker="D", s=22, facecolor=group_color, edgecolor="white", linewidth=0.4, zorder=5)
            else:
                ax.scatter([x], [y], marker="o", s=23, facecolor=group_color, edgecolor="white", linewidth=0.4, zorder=5)
        for alternative in alternatives:
            if alternative["candidate_id"] == candidate_id and not alternative["status"].startswith("promoted"):
                ax.scatter(
                    [int(alternative["candidate_start_relative_to_RT"])], [y], marker="x", s=30,
                    color=RED, linewidth=1.1, zorder=6,
                )

    ax.axvline(0, color="#333333", lw=0.7)
    ax.set_xlim(-1660, 270)
    ax.set_ylim(-0.65, 9.65)
    ax.set_xticks([-1500, -1000, -500, 0])
    ax.set_yticks(range(10))
    ax.set_yticklabels(list(reversed(labels)))
    ax.set_xlabel("Position relative to RT 5′ base (nt)")
    panel_title(ax, "a", "Ten-locus array annotation")
    clean_axis(ax)
    handles = [
        patches.Patch(facecolor="#B8B8B8", label="distal coding anchor"),
        patches.Patch(facecolor="#4A4A4A", label="RT start"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BLUE, markeredgecolor="white", markersize=5.5, label="periodic-chain copy"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=BLUE, markeredgecolor="white", markersize=5.2, label="1-mismatch seed copy"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor="white", markeredgecolor=PURPLE, markersize=5.5, label="weak edge candidate"),
        Line2D([0], [0], marker="x", color=RED, linestyle="none", markersize=5.5, label="unresolved register"),
    ]
    legend_axis.set_axis_off()
    legend_axis.legend(
        handles=handles,
        loc="center",
        ncol=3,
        frameon=False,
        handletextpad=0.45,
        columnspacing=0.8,
    )

    max_unit = 5
    matrix = np.full((len(loci), max_unit), np.nan)
    edge_units: set[tuple[int, int]] = set()
    for row_index, locus in enumerate(loci):
        for unit in units_by_id[locus["candidate_id"]]:
            column = int(unit["unit_ordinal_distal_to_proximal"]) - 1
            matrix[row_index, column] = float(unit["core_to_core_length_nt"])
            right_copy = next(row for row in copies_by_id[locus["candidate_id"]] if row["copy_id"] == unit["right_copy_id"])
            if right_copy["copy_class"] == "group_supported_degenerate_edge":
                edge_units.add((row_index, column))
    cmap = copy.copy(plt.cm.viridis)
    cmap.set_bad("#EEEEEE")
    image = heat.imshow(matrix, aspect="auto", cmap=cmap, vmin=165, vmax=245, interpolation="none")
    for row_index in range(len(loci)):
        for column in range(max_unit):
            if np.isnan(matrix[row_index, column]):
                continue
            heat.text(
                column, row_index, str(int(matrix[row_index, column])), ha="center", va="center",
                fontsize=8, color="white" if matrix[row_index, column] < 195 else "#1A1A1A",
            )
            if (row_index, column) in edge_units:
                heat.add_patch(
                    patches.Rectangle((column - 0.47, row_index - 0.47), 0.94, 0.94, fill=False, edgecolor=PURPLE, linewidth=1.3)
                )
    heat.axhline(6.5, color="white", lw=2.0)
    heat.set_xticks(range(max_unit))
    heat.set_xticklabels([f"U{i}" for i in range(1, max_unit + 1)])
    heat.set_yticks(range(len(loci)))
    heat.set_yticklabels(labels)
    heat.set_xlabel("Core-to-core unit")
    panel_title(heat, "b", "Unit lengths (nt)")
    for spine in heat.spines.values():
        spine.set_visible(False)
    colorbar = fig.colorbar(image, ax=heat, fraction=0.046, pad=0.035)
    colorbar.set_label("Length (nt)")

    # A schematic 212-nt unit makes the nested operational intervals explicit.
    anatomy.set_xlim(-30, 252)
    anatomy.set_ylim(0.0, 5.8)
    anatomy.set_yticks([])
    anatomy.set_xticks([])
    panel_title(anatomy, "c", "Unit interval definitions (212 nt)")
    anatomy.add_patch(patches.Rectangle((0, 3.72), 212, 0.80, facecolor="#E8F0F8", edgecolor=BLUE, linewidth=0.65))
    anatomy.text(106, 4.12, "core-to-core unit [left start, right start)", ha="center", va="center")
    for start in (0, 212):
        anatomy.axvline(start, color="#555555", lw=0.55, ls="--", ymin=0.10, ymax=0.78)
        anatomy.add_patch(patches.Rectangle((start - 26, 3.15), 20, 0.38, facecolor="#D7DDE5", edgecolor="none"))
        anatomy.add_patch(patches.Rectangle((start - 6, 3.15), 26, 0.38, facecolor=TEAL, alpha=0.70, edgecolor="none"))
        anatomy.add_patch(patches.Rectangle((start + 20, 3.15), 20, 0.38, facecolor="#D7DDE5", edgecolor="none"))
        anatomy.add_patch(patches.Rectangle((start - 3, 3.07), 17, 0.54, fill=False, edgecolor="#174A7A", linewidth=1.0))
        anatomy.add_patch(patches.Rectangle((start - 4, 2.57), 22, 0.30, fill=False, edgecolor=ORANGE, linewidth=1.05))
    anatomy_handles = [
        patches.Patch(facecolor="#D7DDE5", edgecolor="none", label="20-nt flank"),
        patches.Patch(facecolor=TEAL, edgecolor="none", alpha=0.70, label="fixed 26-nt window"),
        patches.Patch(facecolor="none", edgecolor="#174A7A", linewidth=1.0, label="selected 17-nt block"),
        patches.Patch(facecolor="none", edgecolor=ORANGE, linewidth=1.05, label="22-nt folding repeat"),
    ]
    anatomy.legend(
        handles=anatomy_handles,
        loc="center",
        bbox_to_anchor=(0.5, 0.86),
        ncol=4,
        frameon=False,
        handlelength=1.35,
        handletextpad=0.38,
        columnspacing=0.85,
    )
    def range_bar(start: float, end: float, y: float, label: str, color: str) -> None:
        anatomy.plot([start, end], [y, y], color=color, lw=0.85)
        anatomy.plot([start, start], [y - 0.18, y + 0.18], color=color, lw=0.85)
        anatomy.plot([end, end], [y - 0.18, y + 0.18], color=color, lw=0.85)
        anatomy.text((start + end) / 2, y + 0.12, label, ha="center", va="bottom", color=color, fontsize=8)

    range_bar(20, 206, 2.30, "between-window spacer [left+20, right−6)", TEAL)
    range_bar(40, 186, 1.60, "trimmed spacer [left+40, right−26)", PURPLE)
    range_bar(-4, 208, 0.90, "one-repeat context", "#555555")
    range_bar(-4, 230, 0.20, "two-repeat context", "#111111")
    for spine in anatomy.spines.values():
        spine.set_visible(False)

    sa1.set_xlim(9400, 12210)
    sa1.set_ylim(0, 1.75)
    sa1.set_yticks([])
    sa1.set_xticks([9447, 9597, 9809, 9980, 10181, 10357, 10646, 12152])
    sa1.set_xticklabels(["9447", "9597", "9809", "9980", "10181", "10357", "10646", "12152"], rotation=35, ha="right")
    sa1.tick_params(axis="x", pad=1)
    panel_title(sa1, "d", "SA1 coordinate context")
    segments = [
        (9447, 9597, "distal residual", "#BDBDBD", "//"),
        (9597, 9809, "U1", "#0072B2", None),
        (9809, 9980, "U2", "#009E73", None),
        (9980, 10181, "U3", "#CC79A7", None),
        (10181, 10357, "U4", "#D55E00", None),
        (10357, 10646, "terminal residual", "#666666", "\\\\"),
        (10646, 12152, "RT", "#333333", None),
    ]
    for start, end, label, color, hatch in segments:
        sa1.add_patch(patches.Rectangle((start, 0.45), end - start, 0.44, facecolor=color, edgecolor="white", linewidth=0.6, hatch=hatch))
        if label in {"U1", "U2", "U3", "U4", "RT"}:
            sa1.text((start + end) / 2, 0.67, label, ha="center", va="center", color="white", fontsize=8)
    sa1.annotate("distal residual", xy=((9447 + 9597) / 2, 0.88), xytext=(9522, 1.00), ha="center", va="bottom", arrowprops=dict(arrowstyle="-", lw=0.7, color="#666666"), fontsize=8)
    sa1.annotate("terminal residual", xy=((10357 + 10646) / 2, 0.88), xytext=(10502, 1.00), ha="center", va="bottom", arrowprops=dict(arrowstyle="-", lw=0.7, color="#444444"), fontsize=8)
    sa1.plot([9597, 10357], [0.33, 0.33], color=TEXT, lw=0.8)
    sa1.plot([9597, 9597], [0.29, 0.37], color=TEXT, lw=0.8)
    sa1.plot([10357, 10357], [0.29, 0.37], color=TEXT, lw=0.8)
    sa1.text((9597 + 10357) / 2, 0.08, "core-to-core U1–U4", ha="center", va="bottom", fontsize=8, color=TEXT)
    sa1.plot([9447, 10646], [1.34, 1.34], color=TEXT, lw=0.9)
    sa1.plot([9447, 9447], [1.26, 1.42], color=TEXT, lw=0.9)
    sa1.plot([10646, 10646], [1.26, 1.42], color=TEXT, lw=0.9)
    sa1.text((9447 + 10646) / 2, 1.46, "published whole array [9447, 10646)", ha="center", va="bottom", fontsize=8)
    sa1.text(10646, 0.12, "RT start", ha="center", va="bottom", fontweight="bold")
    sa1.set_xlabel("SA1 genome coordinate (0-based, half-open intervals)")
    for spine in sa1.spines.values():
        spine.set_visible(False)

    return save_figure(fig, "Figure_1")


def render_figure_2() -> list[Path]:
    loci = read_tsv(LOCI_DIR / "fixed_10_loci.tsv")
    pairwise = read_tsv(P2_DIR / "locus_pairwise_identity.tsv")
    locus_corr = read_tsv(P2_DIR / "locus_pair_unit_correspondence_summary.tsv")
    supported = read_tsv(P2_DIR / "supported_unit_correspondences.tsv")
    labels = [row["label"] for row in loci]
    index = {label: idx for idx, label in enumerate(labels)}
    matrix = np.full((len(labels), len(labels)), np.nan)
    np.fill_diagonal(matrix, 1.0)
    for row in pairwise:
        left = index[row["left_label"]]
        right = index[row["right_label"]]
        value = float(row["intergenic_global_identity"])
        matrix[left, right] = value
        matrix[right, left] = value

    fig = plt.figure(figsize=(170 / 25.4, 154 / 25.4))
    grid = GridSpec(
        2,
        2,
        figure=fig,
        height_ratios=[1.0, 0.78],
        width_ratios=[1.06, 0.94],
        left=0.125,
        right=0.965,
        bottom=0.125,
        top=0.925,
        hspace=0.66,
        wspace=0.48,
    )
    ax1 = fig.add_subplot(grid[0, 0])
    ax2 = fig.add_subplot(grid[0, 1])
    ax3 = fig.add_subplot(grid[1, :])

    identity_cmap = copy.copy(plt.cm.Blues)
    identity_cmap.set_bad(MISSING)
    image = ax1.imshow(np.ma.masked_invalid(matrix), cmap=identity_cmap, vmin=0.0, vmax=1.0, interpolation="none")
    ax1.axhline(6.5, color="white", lw=2)
    ax1.axvline(6.5, color="white", lw=2)
    ax1.set_xticks(range(len(labels)))
    ax1.set_xticklabels(labels, rotation=48, ha="right")
    ax1.set_yticks(range(len(labels)))
    ax1.set_yticklabels(labels)
    panel_title(ax1, "a", "Intervening-interval identity", pad=16)
    for spine in ax1.spines.values():
        spine.set_visible(False)
    colorbar = fig.colorbar(image, ax=ax1, fraction=0.046, pad=0.035, ticks=[0.0, 0.25, 0.5, 0.75, 1.0])
    colorbar.ax.set_title("Identity", fontsize=8, pad=5)

    reference = labels[:7]
    pb_group = labels[7:]
    cross_matrix = np.zeros((len(reference), len(pb_group)))
    text_matrix = [["" for _ in pb_group] for _ in reference]
    for row in locus_corr:
        if row["relationship"] != "cross_group":
            continue
        if row["left_label"] in reference:
            ref_label, pb_label = row["left_label"], row["right_label"]
        else:
            ref_label, pb_label = row["right_label"], row["left_label"]
        i = reference.index(ref_label)
        j = pb_group.index(pb_label)
        cross_matrix[i, j] = int(row["supported_unit_pair_count"])
        ordinals: list[str] = []
        for pair in filter(None, row["supported_ordinal_pairs_left_to_right"].split(";")):
            values = pair.split(":")
            if len(values) == 2 and values[0] == values[1]:
                ordinals.append(f"U{values[0]}")
            elif len(values) == 2:
                ordinals.append(f"U{values[0]}→U{values[1]}")
            else:
                ordinals.append(f"U{pair}")
        text_matrix[i][j] = ", ".join(ordinals)
    cross_cmap = ListedColormap(["#F2F2F2", "#C9D8EC", "#79A4D6", "#315F9A"])
    count_image = ax2.imshow(cross_matrix, cmap=cross_cmap, vmin=-0.5, vmax=3.5, interpolation="none", aspect="auto")
    for i in range(len(reference)):
        for j in range(len(pb_group)):
            text = text_matrix[i][j] if text_matrix[i][j] else "–"
            ax2.text(j, i, text, ha="center", va="center", fontsize=8, color="white" if cross_matrix[i, j] >= 3 else TEXT)
    ax2.add_patch(patches.Rectangle((-0.49, -0.49), 2.98, 0.98, fill=False, edgecolor=RED, linewidth=1.4))
    ax2.set_xticks(range(len(pb_group)))
    ax2.set_xticklabels(pb_group)
    ax2.set_yticks(range(len(reference)))
    ax2.set_yticklabels(reference)
    ax2.set_xlabel("PB50-related locus")
    panel_title(ax2, "b", "Supported unit ordinals", pad=16)
    ax2.text(
        0.0,
        1.02,
        "SA1: zero in all three columns",
        transform=ax2.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        fontweight="normal",
        color=TEXT,
        clip_on=False,
    )
    for spine in ax2.spines.values():
        spine.set_visible(False)
    count_bar = fig.colorbar(count_image, ax=ax2, fraction=0.050, pad=0.040, ticks=[0, 1, 2, 3])
    count_bar.ax.set_title("Count", fontsize=8, pad=5)

    categories = [
        ("Reference group\n(within; n=83)", [r for r in supported if r["relationship"] == "within_group" and r["left_group"] == "P0_seed7"], BLUE),
        ("PB50-related\n(within; n=9)", [r for r in supported if r["relationship"] == "within_group" and r["left_group"] == "PB50_near3"], ORANGE),
        ("Cross-group\n(n=35)", [r for r in supported if r["relationship"] == "cross_group"], TEAL),
    ]
    for category_index, (_, rows, color) in enumerate(categories):
        core_values = [float(row["core_to_core_global_identity"]) for row in rows]
        spacer_values = [float(row["trimmed_spacer_global_identity"]) for row in rows]
        for point_index, (core, spacer) in enumerate(zip(core_values, spacer_values)):
            jitter = 0.075 * math.sin((point_index + 1) * 2.399)
            ax3.plot([category_index - 0.09 + jitter, category_index + 0.09 + jitter], [core, spacer], color=color, alpha=0.16, lw=0.45, zorder=1)
            ax3.scatter([category_index - 0.09 + jitter], [core], s=10, color=color, alpha=0.50, edgecolor="none", zorder=2)
            ax3.scatter([category_index + 0.09 + jitter], [spacer], s=10, facecolor="white", edgecolor=color, linewidth=0.55, alpha=0.70, zorder=2)
        for shift, values in ((-0.09, core_values), (0.09, spacer_values)):
            med = _median(values)
            ax3.plot([category_index + shift - 0.09, category_index + shift + 0.09], [med, med], color="#111111", lw=1.2, zorder=3)
    ax3.set_xticks(range(3))
    ax3.set_xticklabels([item[0] for item in categories])
    ax3.set_ylim(0.45, 1.01)
    ax3.set_ylabel("Global nucleotide identity")
    panel_title(ax3, "c", "Supported-pair identities")
    clean_axis(ax3)
    ax3.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="#666666", markeredgecolor="none", markersize=5, label="core-to-core unit"),
            Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="white", markeredgecolor="#666666", markeredgewidth=0.55, markersize=5, label="trimmed spacer"),
        ],
        frameon=False,
        loc="lower right",
        bbox_to_anchor=(1.0, 1.055),
        ncol=2,
        borderaxespad=0.2,
    )
    fig.text(0.50, 0.020, "All 35 cross-group pairs came from the other six reference loci.", ha="center", va="bottom", fontsize=8, color=TEXT)
    return save_figure(fig, "Figure_2")


def _jitter(index: int, amplitude: float = 0.13) -> float:
    return amplitude * math.sin((index + 1) * 2.399963229728653)


def _plot_null_row(axis, values: list[float], observed: float, median: float, y: float = 0.0) -> None:
    for index, value in enumerate(values):
        axis.scatter(value, y + _jitter(index), s=17, color=GRAY, alpha=0.72, edgecolor="white", linewidth=0.25, zorder=2)
    axis.plot([min(values), max(values)], [y - 0.25, y - 0.25], color=TEXT, lw=1.0, zorder=3)
    axis.plot([median, median], [y - 0.33, y - 0.17], color=TEXT, lw=1.3, zorder=3)
    axis.scatter(observed, y + 0.30, marker="D", s=42, color=RED, edgecolor="white", linewidth=0.45, zorder=4)


def render_figure_3() -> list[Path]:
    rows = read_tsv(P3_DIR / "null_replicates.tsv")
    summary = json.loads((P3_DIR / "p3_calibration_summary.json").read_text(encoding="utf-8"))
    observed = summary["observed_reused_from_P2"]
    distributions = summary["null_distributions"]
    counts = [int(row["supported_pair_count"]) for row in rows]
    u1 = [int(row["u1_same_position_locus_pair_coverage"]) for row in rows]
    u2 = [int(row["u2_same_position_locus_pair_coverage"]) for row in rows]
    identities = [float(row["supported_trimmed_spacer_identity_median"]) for row in rows]

    fig = plt.figure(figsize=(170 / 25.4, 108 / 25.4))
    grid = GridSpec(2, 2, figure=fig, height_ratios=[0.95, 1.05], width_ratios=[1.02, 0.98], left=0.11, right=0.97, bottom=0.18, top=0.90, hspace=0.62, wspace=0.32)
    ax1 = fig.add_subplot(grid[0, :])
    ax2 = fig.add_subplot(grid[1, 0])
    ax3 = fig.add_subplot(grid[1, 1])

    count_dist = distributions["supported_pair_count"]
    _plot_null_row(ax1, counts, observed["supported_pair_count"], count_dist["median"])
    ax1.set_xlim(0, 37)
    ax1.set_ylim(-0.48, 0.55)
    ax1.set_yticks([])
    ax1.set_xlabel("Deduplicated cross-group supported unit pairs per full panel")
    panel_title(ax1, "a", "Supported cross-group pair count")
    ax1.text(0.01, 0.91, "Null median 8; range 4–14; p95 12", transform=ax1.transAxes, ha="left", va="top")
    ax1.text(0.99, 0.91, "Observed 35  |  0/64 null panels reached it", transform=ax1.transAxes, ha="right", va="top", color=RED, fontweight="normal")
    clean_axis(ax1)

    for row_index, (label, values, observed_value, median_value, range_text) in enumerate(
        [
            ("U1 ↔ U1", u1, observed["u1_same_position_locus_pair_coverage"], distributions["u1_same_position_locus_pair_coverage"]["median"], "null 1–7"),
            ("U2 ↔ U2", u2, observed["u2_same_position_locus_pair_coverage"], distributions["u2_same_position_locus_pair_coverage"]["median"], "null 0–4"),
        ]
    ):
        y = 1 - row_index
        for index, value in enumerate(values):
            ax2.scatter(value, y + _jitter(index, 0.10), s=16, color=GRAY, alpha=0.70, edgecolor="white", linewidth=0.25, zorder=2)
        ax2.plot([min(values), max(values)], [y - 0.20, y - 0.20], color=TEXT, lw=0.9)
        ax2.plot([median_value, median_value], [y - 0.27, y - 0.13], color=TEXT, lw=1.2)
        ax2.scatter(observed_value, y, marker="D", s=40, color=RED, edgecolor="white", linewidth=0.45, zorder=4)
        ax2.text(20.8, y + 0.18, f"observed {observed_value} / 21", ha="right", va="bottom", color=RED, fontweight="normal")
        ax2.text(20.8, y - 0.18, range_text, ha="right", va="top", color=TEXT)
    ax2.set_xlim(-0.5, 21.5)
    ax2.set_ylim(-0.48, 1.48)
    ax2.set_yticks([0, 1])
    ax2.set_yticklabels(["U2 ↔ U2", "U1 ↔ U1"])
    ax2.set_xticks([0, 5, 10, 15, 20])
    ax2.set_xlabel("Cross-group locus-pair coverage (of 21)")
    panel_title(ax2, "b", "Same-position coverage")
    clean_axis(ax2)

    identity_dist = distributions["supported_trimmed_spacer_identity_median"]
    _plot_null_row(ax3, identities, observed["supported_trimmed_spacer_identity_median"], identity_dist["median"])
    ax3.set_xlim(0.465, 0.568)
    ax3.set_ylim(-0.48, 0.55)
    ax3.set_yticks([])
    ax3.set_xlabel("Median identity (displayed 0.465–0.568)")
    panel_title(ax3, "c", "Trimmed-spacer identity")
    ax3.text(0.02, 0.91, "Null median 0.504\nrange 0.478–0.526", transform=ax3.transAxes, ha="left", va="top", color=TEXT)
    ax3.text(0.98, 0.91, "Observed 0.557", transform=ax3.transAxes, ha="right", va="top", color=RED, fontweight="normal")
    clean_axis(ax3)
    null_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=GRAY, markeredgecolor="white", markeredgewidth=0.25, markersize=5.2, label="64 null panels"),
        Line2D([0], [0], marker="D", linestyle="none", markerfacecolor=RED, markeredgecolor="white", markeredgewidth=0.45, markersize=5.5, label="observed"),
        Line2D([0, 1], [0, 0], color=TEXT, lw=1.0, marker="|", markersize=7, markeredgewidth=1.2, label="null min–max / median"),
    ]
    fig.legend(handles=null_handles, loc="lower center", bbox_to_anchor=(0.5, 0.018), ncol=3, frameon=False, handletextpad=0.45, columnspacing=1.05)
    return save_figure(fig, "Figure_3")


def render_figure_4() -> list[Path]:
    copies = read_tsv(P2_DIR / "unified_repeat_copies.tsv")
    edge_scan = read_tsv(P2_DIR / "edge_scan_audit.tsv")
    copy_map = read_tsv(P2_DIR / "cross_group_medoid_copy_mapping.tsv")
    blocks = read_tsv(P2_DIR / "cross_group_medoid_array_alignment_blocks.tsv")

    fig = plt.figure(figsize=(170 / 25.4, 138 / 25.4))
    grid = GridSpec(
        2,
        2,
        figure=fig,
        height_ratios=[0.80, 1.0],
        width_ratios=[0.88, 1.22],
        left=0.29,
        right=0.97,
        bottom=0.105,
        top=0.93,
        hspace=0.50,
        wspace=0.42,
    )
    ax1 = fig.add_subplot(grid[0, :])
    ax2 = fig.add_subplot(grid[1, 0])
    ax3 = fig.add_subplot(grid[1, 1])
    # Reserve physical left clearance for the longest profile label while
    # keeping the panel-a right boundary and the lower-panel topology fixed.
    ax1_position = ax1.get_position()
    ax1.set_position(
        [ax1_position.x0 + 0.012, ax1_position.y0, ax1_position.width - 0.012, ax1_position.height]
    )

    near_copies = [row for row in copies if row["analysis_subgroup"] == "PB50_near3"]
    all_sequences = [row["fixed26_sequence"] for row in near_copies]
    consensus = "".join(Counter(sequence[column] for sequence in all_sequences).most_common(1)[0][0] for column in range(26))
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in near_copies:
        grouped[row["fixed26_sequence"]].append(row)
    groups = sorted(grouped.items(), key=lambda item: (min(int(r["calibrated_copy_start_relative_to_RT"]) for r in item[1]), item[0]))
    sequences = [item[0] for item in groups]
    mismatch = np.array([[int(base != consensus[column]) for column, base in enumerate(sequence)] for sequence in sequences])
    ax1.imshow(mismatch, aspect="auto", cmap=ListedColormap(["#F5E8D8", RED]), vmin=0, vmax=1, interpolation="none")
    for row_index, sequence in enumerate(sequences):
        for column, base in enumerate(sequence):
            ax1.text(column, row_index, base, ha="center", va="center", family="monospace", fontsize=8, color="white" if mismatch[row_index, column] else TEXT)
    ax1.add_patch(patches.Rectangle((5.5, -0.48), 10, len(sequences) - 0.04, fill=False, edgecolor="#2F2F2F", linewidth=1.0))
    group_labels: list[str] = []
    abbreviations = {"AH12": "AH", "Machias": "Ma", "PB50": "PB"}
    for _, rows in groups:
        members = "/".join(abbreviations[r] for r in sorted({r["label"] for r in rows}))
        copies_here = "/".join(sorted({f"C{r['copy_ordinal_distal_to_proximal']}" for r in rows}))
        starts = "/".join(_format_numeric_minus(v) for v in sorted({int(r["calibrated_copy_start_relative_to_RT"]) for r in rows}))
        group_labels.append(f"{copies_here} {starts} | {members} (n={len(rows)})")
    ax1.set_yticks(range(len(groups)))
    ax1.set_yticklabels(group_labels)
    ax1.set_xticks(range(0, 26, 2))
    ax1.set_xticklabels([_format_numeric_minus(index - 6) for index in range(0, 26, 2)])
    ax1.set_xlabel("Position relative to selected 10-nt seed start (boxed: seed)")
    panel_title(ax1, "a", "PB50-related repeat windows", pad=14)
    ax1.text(1.0, 1.025, "5 unique profiles; 12 copies", transform=ax1.transAxes, ha="right", va="bottom", fontsize=8, color=TEXT)
    for spine in ax1.spines.values():
        spine.set_visible(False)

    ref_proximal = [row for row in edge_scan if row["analysis_subgroup"] == "P0_seed7" and row["side"] == "proximal"]
    label_order = ["SA1", "MarsHill", "Madawaska", "LY01", "S6", "PALS_2", "UFV_DC4"]
    ref_proximal.sort(key=lambda row: label_order.index(row["label"]))
    for row_index, row in enumerate(ref_proximal):
        accepted = row["meets_single_locus_candidate_floor"] == "true"
        score = float(row["profile_match_26"])
        marker = "^" if accepted else "o"
        color = PURPLE if accepted else "#8F8F8F"
        ax2.scatter([score], [row_index], marker=marker, s=43 if accepted else 27, facecolor="white" if accepted else color, edgecolor=color, linewidth=1.0, zorder=3)
        mismatch_count = int(row["seed_mismatches"])
        ax2.text(
            18.35,
            row_index,
            f"{mismatch_count}; {_format_numeric_minus(row['candidate_start_relative_to_RT'])} nt",
            ha="left",
            va="center",
            color=color,
            fontsize=8,
            zorder=3,
        )
    ax2.axvline(18, color="#555555", lw=0.8, linestyle="--", zorder=0)
    ax2.text(18.35, -0.62, "mismatches; start", ha="left", va="center", fontsize=8, color=TEXT)
    ax2.set_xlim(12.5, 27.0)
    ax2.set_ylim(len(ref_proximal) - 0.4, -0.8)
    ax2.set_yticks(range(len(ref_proximal)))
    ax2.set_yticklabels([row["label"] for row in ref_proximal])
    ax2.set_xlabel("Best 26-nt profile-match count")
    panel_title(ax2, "b", "Proximal-edge screen", pad=16)
    clean_axis(ax2)

    ref_native = [row for row in copy_map if row["source_label"] == "PALS_2"]
    pb_mapped = [row for row in copy_map if row["source_label"] == "PB50"]
    for y, display_label, color in ((1.0, "PALS_2 native", BLUE), (0.0, "PB50 mapped", ORANGE)):
        ax3.plot([-1160, -20], [y, y], color="#BDBDBD", lw=1.0)
        ax3.text(-1175, y, display_label, ha="right", va="center", color=color, fontsize=8)
    ref_positions: dict[int, int] = {}
    for row in ref_native:
        position = int(row["source_copy_relative_start"])
        ordinal = int(row["source_copy_ordinal"])
        ref_positions[ordinal] = position
        is_edge = ordinal == max(int(item["source_copy_ordinal"]) for item in ref_native)
        ax3.scatter([position], [1.0], marker="^" if is_edge else "o", s=37 if is_edge else 28, facecolor="white" if is_edge else BLUE, edgecolor=PURPLE if is_edge else "white", linewidth=1.0 if is_edge else 0.4, zorder=4)
        ax3.text(position, 1.17, f"C{ordinal}", ha="center", va="bottom", fontsize=8)
    for row in pb_mapped:
        if row["mapping_status"] != "aligned_base":
            continue
        position = int(row["mapped_target_relative_start"])
        ordinal = int(row["source_copy_ordinal"])
        ax3.scatter([position], [0.0], marker="s", s=28, facecolor=ORANGE, edgecolor="white", linewidth=0.4, zorder=4)
        ax3.text(position, -0.17, f"C{ordinal}", ha="center", va="top", fontsize=8)
    for row in blocks:
        if row["block_type"] == "query_gap" and int(row["alignment_columns"]) >= 10:
            left = int(row["target_relative_start"])
            right = int(row["target_relative_end_exclusive"])
            ax3.add_patch(patches.Rectangle((left, 0.37), right - left, 0.26, facecolor=RED, alpha=0.36, edgecolor="none"))
    gap_copy = next(row for row in ref_native if row["mapping_status"] == "maps_to_gap")
    gap_position = int(gap_copy["source_copy_relative_start"])
    ax3.scatter([gap_position], [1.0], marker="x", s=58, color=RED, linewidth=1.3, zorder=5)
    ax3.text(
        0.96,
        0.53,
        "PB50 gap blocks\n14, 50, 61, 24 nt",
        transform=ax3.transAxes,
        ha="right",
        va="bottom",
        color=RED,
        fontsize=8,
        zorder=7,
        bbox=dict(facecolor="white", edgecolor="none", pad=1.2),
    )
    ax3.annotate("C4 maps to gap", xy=(gap_position, 0.86), xytext=(-790, 1.47), arrowprops=dict(arrowstyle="-", lw=0.7, color=RED), fontsize=8, color=RED)
    for ordinal in (1, 2, 3):
        pb = next((row for row in pb_mapped if int(row["source_copy_ordinal"]) == ordinal), None)
        if pb and pb["mapping_status"] == "aligned_base" and ordinal in ref_positions:
            ax3.plot([ref_positions[ordinal], int(pb["mapped_target_relative_start"])], [0.96, 0.04], color="#777777", lw=0.55, alpha=0.7)
    pb4 = next(row for row in pb_mapped if int(row["source_copy_ordinal"]) == 4)
    ax3.plot([ref_positions[5], int(pb4["mapped_target_relative_start"])], [0.96, 0.04], color="#777777", lw=0.55, alpha=0.7)
    ax3.set_xlim(-1210, 10)
    ax3.set_ylim(-0.42, 1.62)
    ax3.set_yticks([])
    ax3.set_xlabel("Coordinate relative to RT (nt)")
    panel_title(ax3, "c", "PALS_2–PB50 boundary map", pad=16)
    ax3.spines["left"].set_visible(False)
    clean_axis(ax3)
    return save_figure(fig, "Figure_4")


def main() -> None:
    configure()
    outputs: list[Path] = []
    outputs.extend(render_figure_1())
    outputs.extend(render_figure_2())
    outputs.extend(render_figure_3())
    outputs.extend(render_figure_4())
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
