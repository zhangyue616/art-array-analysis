#!/usr/bin/env python3
"""Render revised Figure 5, Figure 6, and Supplementary Figure S1.

All panels consume saved P6 tables. The identity-range panels use saved
exploratory diagnostic tables.
"""

from __future__ import annotations

from collections import defaultdict
import json
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

from common import (
    BLUE,
    GRAY,
    GRID,
    ORANGE,
    RED,
    ROOT,
    TEXT,
    clean_axis,
    configure,
    panel_title,
    read_tsv,
    save_figure,
)


RNA_DIR = ROOT / "data" / "processed" / "p6_rna"
DIAG_DIR = ROOT / "data" / "processed" / "p6_revision_20260925"
PAIR_AUDIT = DIAG_DIR / "s1_two_repeat_pair_audit.tsv"
LOCUS_DIAGNOSTIC = DIAG_DIR / "s1_locus_pair_identity_range_diagnostic.tsv"

NULL = "#8A96A3"
PAIR_COLORS = {
    "spacer_spacer": "#D06B32",
    "repeat_repeat": "#3D6FA6",
    "mixed": "#8C5EA8",
}
PB_COLORS = {"AH12": "#3B78A7", "Machias": "#2F8F6B", "PB50": "#D37A32"}


def _median(values: list[float]) -> float:
    values = sorted(values)
    middle = len(values) // 2
    return values[middle] if len(values) % 2 else 0.5 * (values[middle - 1] + values[middle])


def _jitter(index: int, amplitude: float = 0.08) -> float:
    return amplitude * math.sin((index + 1) * 2.399963229728653)


def _panel_header(axis, letter: str, title: str, detail: str, letter_offset_pt: float = -13) -> None:
    """Place a two-line panel header in reserved space above an axis."""
    panel_title(axis, letter, title, pad=18, letter_offset_pt=letter_offset_pt)
    axis.text(
        0.0,
        1.035,
        detail,
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        color=TEXT,
        clip_on=False,
    )


def render_figure_5() -> list:
    pairs = read_tsv(PAIR_AUDIT)
    diagnostics = [
        row
        for row in read_tsv(LOCUS_DIAGNOSTIC)
        if row["identity_field"] == "context_global_identity"
        and row["range_scope"] == "global_supported_min_max"
        and row["classification"] == "evaluable"
    ]
    if len(pairs) != 252 or len(diagnostics) != 18:
        raise ValueError(f"Expected 252 pair rows and 18 evaluable loci; got {len(pairs)} and {len(diagnostics)}")

    supported = [row for row in pairs if row["supported_correspondence"] == "true"]
    other = [row for row in pairs if row["supported_correspondence"] == "false"]
    other_range = [row for row in other if row["other_in_global_supported_range__context_global_identity"] == "true"]
    if (len(supported), len(other), len(other_range)) != (35, 217, 86):
        raise ValueError("Figure 5 diagnostic counts do not match the accepted science receipt")

    fig = plt.figure(figsize=(170 / 25.4, 170 / 25.4))
    grid = GridSpec(
        2,
        2,
        figure=fig,
        height_ratios=(0.88, 1.12),
        left=0.18,
        right=0.96,
        bottom=0.09,
        top=0.88,
        hspace=0.76,
        wspace=0.40,
    )
    ax_all = fig.add_subplot(grid[0, 0])
    ax_range = fig.add_subplot(grid[0, 1])
    ax_locus = fig.add_subplot(grid[1, :])
    locus_position = ax_locus.get_position()
    ax_locus.set_position(
        [locus_position.x0 + 0.025, locus_position.y0, locus_position.width - 0.025, locus_position.height]
    )

    # a: all 252 two-repeat pair combinations.
    ax_all.scatter(
        [float(row["context_global_identity"]) for row in other],
        [float(row["spacer_spacer_bpp_overlap"]) for row in other],
        s=16,
        color="#B9BFC5",
        alpha=0.68,
        edgecolors="none",
        label="Other (n=217)",
        zorder=2,
    )
    ax_all.scatter(
        [float(row["context_global_identity"]) for row in supported],
        [float(row["spacer_spacer_bpp_overlap"]) for row in supported],
        s=31,
        marker="D",
        color=BLUE,
        alpha=0.92,
        edgecolors="white",
        linewidths=0.45,
        label="Supported (n=35)",
        zorder=3,
    )
    ax_all.set_xlabel("Context global identity")
    ax_all.set_ylabel("Spacer–spacer BPP overlap")
    _panel_header(ax_all, "a", "All pair combinations", "Spearman ρ = 0.615 (n = 252)")
    ax_all.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.20),
        ncol=2,
        frameon=False,
        handletextpad=0.35,
        columnspacing=0.9,
    )
    clean_axis(ax_all, "both")

    # b: raw outcomes after restricting other pairs to the supported identity range.
    range_groups = [("Supported\n(n=35)", supported, BLUE), ("Other in range\n(n=86)", other_range, GRAY)]
    for x_position, (label, rows, color) in enumerate(range_groups, 1):
        values = [float(row["spacer_spacer_bpp_overlap"]) for row in rows]
        for index, value in enumerate(values):
            ax_range.scatter(
                x_position + _jitter(index, 0.12),
                value,
                s=17,
                color=color,
                alpha=0.65,
                edgecolors="none",
                zorder=2,
            )
        median = _median(values)
        ax_range.plot([x_position - 0.23, x_position + 0.23], [median, median], color=TEXT, lw=1.35, zorder=4)
    ax_range.set_xlim(0.55, 2.45)
    ax_range.set_xticks([1, 2])
    ax_range.set_xticklabels([item[0] for item in range_groups])
    ax_range.set_ylabel("Spacer–spacer BPP overlap")
    _panel_header(
        ax_range,
        "b",
        "Shared-range restriction",
        "0.588–0.667; medians 0.085 / 0.104",
    )
    clean_axis(ax_range, "y")

    # c: locus-level effects under the same global range restriction.
    ordered = sorted(diagnostics, key=lambda row: float(row["supported_minus_range_other_bpp_median"]))
    positions = list(range(len(ordered)))
    effects = [float(row["supported_minus_range_other_bpp_median"]) for row in ordered]
    colors = [PB_COLORS[row["pb_label"]] for row in ordered]
    labels = [row["locus_pair_id"].replace("__", "–") for row in ordered]
    ax_locus.axvline(0, color="#707070", lw=0.85, zorder=1)
    for y_position, effect, color in zip(positions, effects, colors):
        ax_locus.plot([0, effect], [y_position, y_position], color="#C6C6C6", lw=0.8, zorder=1)
        ax_locus.scatter(effect, y_position, s=34, color=color, edgecolor="white", linewidth=0.45, zorder=3)
    ax_locus.set_yticks(positions)
    ax_locus.set_yticklabels(labels)
    ax_locus.set_xlabel("Supported median − identity-range other median (BPP overlap)")
    _panel_header(
        ax_locus,
        "c",
        "Locus-pair shared-range effects",
        "18 evaluable pairs; 9/18 > 0; median −0.00183",
    )
    locus_color_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=color, markeredgecolor="white", markeredgewidth=0.45, markersize=5.5, label=label)
        for label, color in PB_COLORS.items()
    ]
    ax_locus.legend(
        handles=locus_color_handles,
        loc="lower right",
        bbox_to_anchor=(1.0, 1.018),
        ncol=3,
        frameon=False,
        handletextpad=0.30,
        columnspacing=0.75,
        borderaxespad=0.0,
    )
    clean_axis(ax_locus, "x")
    fig.text(
        0.50,
        0.012,
        "Exploratory diagnostic.",
        ha="center",
        va="bottom",
        fontsize=8,
        color=TEXT,
    )
    return save_figure(fig, "Figure_5")


def _plot_null(axis, distributions, categories, ylim, annotation):
    for x_position, (label, field, observed) in enumerate(categories, 1):
        values = [float(row[field]) for row in distributions]
        for index, value in enumerate(values):
            axis.scatter(x_position + _jitter(index, 0.13), value, s=14, color=NULL, alpha=0.58, edgecolor="none", clip_on=False, zorder=2)
        axis.plot([x_position - 0.24, x_position + 0.24], [_median(values), _median(values)], color=TEXT, lw=1.1, zorder=3)
        axis.plot([x_position, x_position], [min(values), max(values)], color=TEXT, lw=0.65, zorder=1)
        axis.scatter(x_position, observed, marker="D", s=39, color=RED, edgecolor="white", linewidth=0.45, clip_on=False, zorder=4)
    axis.set_xlim(0.55, len(categories) + 0.45)
    axis.set_ylim(*ylim)
    axis.set_xticks(list(range(1, len(categories) + 1)))
    axis.set_xticklabels([item[0] for item in categories])
    axis.set_ylabel("BPP overlap")
    axis.text(
        0.98,
        1.025,
        annotation,
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color=RED,
        clip_on=False,
        zorder=6,
    )
    clean_axis(axis, "y")


def render_figure_6() -> list:
    rows = read_tsv(RNA_DIR / "dinucleotide_null_replicates.tsv")
    summary = json.loads((RNA_DIR / "p6_rna_summary.json").read_text(encoding="utf-8"))
    two = [row for row in rows if row["context"] == "two_repeat"]
    one = [row for row in rows if row["context"] == "one_repeat"]
    obs_two = summary["observed"]["two_repeat"]
    obs_one = summary["observed"]["one_repeat"]

    fig = plt.figure(figsize=(170 / 25.4, 136 / 25.4))
    grid = GridSpec(2, 2, figure=fig, left=0.10, right=0.97, bottom=0.15, top=0.88, hspace=0.72, wspace=0.32)
    axes = [fig.add_subplot(grid[index // 2, index % 2]) for index in range(4)]
    _plot_null(
        axes[0],
        two,
        [
            ("All 35", "spacer_spacer_bpp_overlap_median", obs_two["supported_spacer_spacer_bpp_overlap_median"]),
            ("U1 (15)", "u1_spacer_spacer_bpp_overlap_median", obs_two["u1_spacer_spacer_bpp_overlap_median"]),
            ("U2 (18)", "u2_spacer_spacer_bpp_overlap_median", obs_two["u2_spacer_spacer_bpp_overlap_median"]),
        ],
        (0.0, 0.105),
        "0/64 null ≥ observed",
    )
    _plot_null(
        axes[1],
        one,
        [
            ("All 35", "spacer_spacer_bpp_overlap_median", obs_one["supported_spacer_spacer_bpp_overlap_median"]),
            ("U1 (15)", "u1_spacer_spacer_bpp_overlap_median", obs_one["u1_spacer_spacer_bpp_overlap_median"]),
            ("U2 (18)", "u2_spacer_spacer_bpp_overlap_median", obs_one["u2_spacer_spacer_bpp_overlap_median"]),
        ],
        (0.0, 0.105),
        "0/64 null ≥ observed",
    )
    _plot_null(
        axes[2],
        two,
        [("Mixed", "mixed_bpp_overlap_median", obs_two["supported_mixed_bpp_overlap_median"])],
        (0.0, 0.105),
        "0/64 null ≥ observed",
    )
    _plot_null(
        axes[3],
        two,
        [("Repeat", "repeat_repeat_bpp_overlap_median", obs_two["supported_repeat_repeat_bpp_overlap_median"])],
        (0.15, 0.98),
        "20/64 null ≥ observed",
    )
    titles = (
        "Two-repeat spacer",
        "One-repeat sensitivity",
        "Repeat–spacer mixed",
        "Repeat–repeat (different scale)",
    )
    for letter, axis, title in zip("abcd", axes, titles):
        panel_title(axis, letter, title, pad=16)
    shared_ticks = [0.0, 0.02, 0.04, 0.06, 0.08, 0.10]
    for axis in axes[:3]:
        axis.set_ylim(0.0, 0.105)
        axis.set_yticks(shared_ticks)
        axis.set_yticklabels(["0", "0.02", "0.04", "0.06", "0.08", "0.10"])
    null_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=NULL, markeredgecolor="none", markersize=5.2, label="64 null replicates"),
        Line2D([0], [0], marker="D", linestyle="none", markerfacecolor=RED, markeredgecolor="white", markeredgewidth=0.45, markersize=5.5, label="observed"),
        Line2D([0, 1], [0, 0], color=TEXT, lw=1.0, marker="|", markersize=7, markeredgewidth=1.2, label="null min–max / median"),
    ]
    fig.legend(handles=null_handles, loc="lower center", bbox_to_anchor=(0.5, 0.018), ncol=3, frameon=False, handletextpad=0.45, columnspacing=1.05)
    return save_figure(fig, "Figure_6")


def _parse_ranges(value: str) -> list[tuple[int, int]]:
    return [tuple(int(part) for part in item.split("-")) for item in value.split(";")]


def _probability_size(probability: float) -> float:
    return 2.5 + 34.0 * math.sqrt(probability)


def render_figure_s1() -> list:
    edges = read_tsv(RNA_DIR / "representative_alignment_bpp.tsv")
    metadata = read_tsv(RNA_DIR / "representative_alignment_metadata.tsv")
    by_ordinal: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in edges:
        by_ordinal[row["ordinal"]].append(row)

    fig = plt.figure(figsize=(170 / 25.4, 122 / 25.4))
    grid = GridSpec(1, 2, figure=fig, left=0.085, right=0.985, bottom=0.22, top=0.78, wspace=0.24)
    axes = []
    for index, meta in enumerate(metadata):
        axis = fig.add_subplot(grid[0, index])
        axes.append(axis)
        alignment_length = int(meta["alignment_length"])
        for row in by_ordinal[meta["ordinal"]]:
            left = int(row["left_alignment_column_1based"])
            right = int(row["right_alignment_column_1based"])
            probability = float(row["probability"])
            x_value, y_value = (left, right) if row["side"] == "PALS_2" else (right, left)
            axis.scatter(
                x_value,
                y_value,
                s=_probability_size(probability),
                color=PAIR_COLORS[row["pair_class"]],
                alpha=0.55,
                edgecolor="none",
                zorder=2,
            )
        axis.plot([1, alignment_length], [1, alignment_length], color="#6E6E6E", lw=0.65, zorder=1)
        for start, end in _parse_ranges(meta["left_repeat_columns_1based"]):
            axis.plot([start, end], [alignment_length + 4, alignment_length + 4], color=PAIR_COLORS["repeat_repeat"], lw=3.0, solid_capstyle="butt", clip_on=False)
        for start, end in _parse_ranges(meta["right_repeat_columns_1based"]):
            axis.plot([start, end], [-4, -4], color=PAIR_COLORS["repeat_repeat"], lw=3.0, solid_capstyle="butt", clip_on=False)
        axis.set_xlim(0, alignment_length + 1)
        axis.set_ylim(0, alignment_length + 1)
        axis.set_aspect("equal", adjustable="box")
        axis.set_xlabel("Aligned context column")
        if index == 0:
            axis.set_ylabel("Aligned context column")
        panel_title(axis, chr(ord("a") + index), f"{meta['ordinal']} supported pair", pad=18)
        axis.text(
            0.0,
            1.025,
            f"BPP overlap: spacer {float(meta['spacer_spacer_bpp_overlap']):.3f}; repeat {float(meta['repeat_repeat_bpp_overlap']):.3f}",
            transform=axis.transAxes,
            ha="left",
            va="bottom",
            fontsize=8,
            color=TEXT,
            clip_on=False,
        )
        clean_axis(axis)

    class_handles = [
        Line2D([0], [0], marker="o", lw=0, markersize=5.5, markerfacecolor=PAIR_COLORS[key], markeredgecolor="none", label=label)
        for key, label in (("spacer_spacer", "spacer–spacer"), ("repeat_repeat", "repeat–repeat"), ("mixed", "mixed"))
    ]
    size_handles = [
        plt.scatter([], [], s=_probability_size(value), color="#666666", alpha=0.65, edgecolor="none", label=f"{value:.2f}")
        for value in (0.01, 0.10, 0.50)
    ]
    fig.legend(handles=class_handles, title="Pair class", loc="lower left", ncol=3, frameon=False, bbox_to_anchor=(0.08, 0.055), handletextpad=0.4, columnspacing=1.10)
    fig.legend(handles=size_handles, title="Pair probability", loc="lower right", ncol=3, frameon=False, bbox_to_anchor=(0.96, 0.055), handletextpad=0.35, columnspacing=0.90)
    fig.text(0.50, 0.955, "Triangle source: PALS_2 above diagonal  |  PB50 below diagonal", ha="center", va="top", fontsize=8, color=TEXT)
    fig.text(0.50, 0.015, "Displayed points: probability ≥ 0.01. Overlap calculations used saved edges at probability ≥ 10⁻⁶.", ha="center", va="bottom", fontsize=8, color=TEXT)
    return save_figure(fig, "Figure_S1")


def main() -> None:
    configure()
    outputs = []
    outputs.extend(render_figure_5())
    outputs.extend(render_figure_6())
    outputs.extend(render_figure_s1())
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
