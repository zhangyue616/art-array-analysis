#!/usr/bin/env python3
"""Render the three P6 RNA ensemble figures from completed tables."""

from __future__ import division

import csv
import json
import math
import os
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D


TEXT = "#252525"
GRID = "#D9D9D9"
NULL = "#8A96A3"
OBSERVED = "#B5473C"
PANEL_COLORS = {
    "AH12": "#3B78A7",
    "Machias": "#2F8F6B",
    "PB50": "#D37A32",
}
PAIR_COLORS = {
    "spacer_spacer": "#D06B32",
    "repeat_repeat": "#3D6FA6",
    "mixed": "#8C5EA8",
}


def configure():
    matplotlib.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.titlesize": 9.5,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "legend.fontsize": 8.0,
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
    if not os.path.isdir(figure_dir):
        os.makedirs(figure_dir)
    outputs = []
    for extension in ("pdf", "svg", "png"):
        path = os.path.join(figure_dir, stem + "." + extension)
        kwargs = {"facecolor": "white"}
        if extension == "png":
            kwargs["dpi"] = 300
        fig.savefig(path, **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def deterministic_jitter(index, amplitude=0.08):
    return amplitude * math.sin((index + 1) * 2.399963229728653)


def clean_axis(axis, grid_axis=None):
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    if grid_axis:
        axis.grid(True, axis=grid_axis, color=GRID, lw=0.45, alpha=0.65, zorder=0)


def panel_letter(axis, text):
    axis.text(-0.12, 1.045, text, transform=axis.transAxes,
              ha="left", va="bottom", fontweight="bold", fontsize=10)


def float_value(row, key):
    return float(row[key])


def plot_locus_scatter(axis, rows, context, title):
    subset = [row for row in rows
              if row["context"] == context and int(row["supported_pairs"]) > 0]
    values = []
    for row in subset:
        x = float_value(row, "misaligned_spacer_bpp_overlap_median")
        y = float_value(row, "supported_spacer_bpp_overlap_median")
        values.append((x, y))
        axis.scatter(x, y, s=32, color=PANEL_COLORS[row["pb_label"]],
                     edgecolor="white", linewidth=0.45, zorder=3)
        if y <= x:
            axis.annotate(row["locus_pair_id"].replace("__", "–"), (x, y),
                          xytext=(4, -7), textcoords="offset points", fontsize=8.0,
                          color=TEXT)
    upper = max(max(item) for item in values) * 1.10
    upper = max(0.13, upper)
    axis.plot([0, upper], [0, upper], color="#626262", lw=0.8, ls="--", zorder=1)
    axis.set_xlim(0, upper)
    axis.set_ylim(0, upper)
    axis.set_xlabel("Other unit combinations: median overlap")
    axis.set_ylabel("Supported-unit median overlap")
    axis.set_title(title, loc="left", fontweight="bold", fontsize=8.9)
    clean_axis(axis, "both")


def figure_specificity(data_dir, figure_dir):
    locus_rows = read_tsv(os.path.join(data_dir, "locus_pair_structure_summary.tsv"))
    pair_rows = read_tsv(os.path.join(data_dir, "unit_pair_structure_similarity.tsv"))

    fig = plt.figure(figsize=(170 / 25.4, 124 / 25.4))
    grid = GridSpec(2, 2, height_ratios=[1.0, 0.78], left=0.095, right=0.985,
                    bottom=0.15, top=0.925, hspace=0.52, wspace=0.31)
    ax1 = fig.add_subplot(grid[0, 0])
    ax2 = fig.add_subplot(grid[0, 1])
    ax3 = fig.add_subplot(grid[1, :])

    plot_locus_scatter(ax1, locus_rows, "two_repeat",
                       "Two-repeat: 18/18 locus pairs positive")
    plot_locus_scatter(ax2, locus_rows, "one_repeat",
                       "One-repeat: 17/18 positive")
    panel_letter(ax1, "a")
    panel_letter(ax2, "b")

    supported = [row for row in pair_rows
                 if row["context"] == "two_repeat" and row["supported_correspondence"] == "true"]
    ordinal_rows = defaultdict(list)
    for row in supported:
        ordinal_rows[int(row["p0_ordinal"])].append(row)
    for ordinal in (1, 2, 3):
        rows = ordinal_rows[ordinal]
        values = sorted(float_value(row, "spacer_spacer_bpp_overlap") for row in rows)
        for index, row in enumerate(rows):
            x = ordinal + deterministic_jitter(index, 0.11)
            ax3.scatter(x, float_value(row, "spacer_spacer_bpp_overlap"), s=27,
                        color=PANEL_COLORS[row["pb_label"]], edgecolor="white",
                        linewidth=0.40, alpha=0.88, zorder=3)
        midpoint = values[len(values) // 2] if len(values) % 2 else 0.5 * (values[len(values)//2 - 1] + values[len(values)//2])
        ax3.plot([ordinal - 0.22, ordinal + 0.22], [midpoint, midpoint],
                 color=TEXT, lw=1.3, zorder=4)
        ax3.text(ordinal, 0.238, "n={}".format(len(values)), ha="center", va="top", fontsize=8.0)
    ax3.set_xlim(0.55, 3.45)
    ax3.set_ylim(-0.005, 0.245)
    ax3.set_xticks([1, 2, 3])
    ax3.set_xticklabels(["U1", "U2", "U3"])
    ax3.set_ylabel("Spacer–spacer BPP overlap")
    ax3.set_title("Supported unit pairs in the primary two-repeat context",
                  loc="left", fontweight="bold")
    clean_axis(ax3, "y")
    panel_letter(ax3, "c")

    handles = [Line2D([0], [0], marker="o", lw=0, markersize=5.5,
                      markerfacecolor=PANEL_COLORS[name], markeredgecolor="white", label=name)
               for name in ("AH12", "Machias", "PB50")]
    ax1.legend(handles=handles, title="PB50-near locus", loc="upper left",
               frameon=False, borderpad=0.1, handletextpad=0.4)
    fig.text(0.50, 0.024,
             "Panels a–b: locus-pair medians. Panel c: 35 sequence-supported pairs.\nPair points share loci and are descriptive.",
             ha="center", va="bottom", fontsize=8.0, color=TEXT, linespacing=1.25)
    return save_figure(fig, figure_dir, "figure_p6_1_supported_vs_wrong_units")


def plot_null_rows(axis, distributions, observed, categories, ylim, title):
    for x_index, (label, field, observed_value) in enumerate(categories, 1):
        values = [float(row[field]) for row in distributions]
        for index, value in enumerate(values):
            axis.scatter(x_index + deterministic_jitter(index, 0.13), value,
                         s=14, color=NULL, alpha=0.58, edgecolor="none", zorder=2)
        sorted_values = sorted(values)
        median = 0.5 * (sorted_values[31] + sorted_values[32])
        axis.plot([x_index - 0.24, x_index + 0.24], [median, median], color=TEXT, lw=1.1, zorder=3)
        axis.plot([x_index, x_index], [min(values), max(values)], color=TEXT, lw=0.65, zorder=1)
        axis.scatter(x_index, observed_value, marker="D", s=39, color=OBSERVED,
                     edgecolor="white", linewidth=0.45, zorder=4)
    axis.set_xlim(0.55, len(categories) + 0.45)
    axis.set_ylim(*ylim)
    axis.set_xticks(range(1, len(categories) + 1))
    axis.set_xticklabels([item[0] for item in categories])
    axis.set_ylabel("BPP overlap")
    axis.set_title(title, loc="left", fontweight="bold")
    clean_axis(axis, "y")


def figure_null(data_dir, figure_dir):
    rows = read_tsv(os.path.join(data_dir, "dinucleotide_null_replicates.tsv"))
    with open(os.path.join(data_dir, "p6_rna_summary.json"), "r", encoding="utf-8") as handle:
        summary = json.load(handle)
    two = [row for row in rows if row["context"] == "two_repeat"]
    one = [row for row in rows if row["context"] == "one_repeat"]
    obs_two = summary["observed"]["two_repeat"]
    obs_one = summary["observed"]["one_repeat"]

    fig = plt.figure(figsize=(170 / 25.4, 126 / 25.4))
    grid = GridSpec(2, 2, left=0.095, right=0.98, bottom=0.145, top=0.93,
                    hspace=0.50, wspace=0.30)
    ax1 = fig.add_subplot(grid[0, 0])
    ax2 = fig.add_subplot(grid[0, 1])
    ax3 = fig.add_subplot(grid[1, 0])
    ax4 = fig.add_subplot(grid[1, 1])

    plot_null_rows(
        ax1, two, obs_two,
        [
            ("All 35", "spacer_spacer_bpp_overlap_median", obs_two["supported_spacer_spacer_bpp_overlap_median"]),
            ("U1 (15)", "u1_spacer_spacer_bpp_overlap_median", obs_two["u1_spacer_spacer_bpp_overlap_median"]),
            ("U2 (18)", "u2_spacer_spacer_bpp_overlap_median", obs_two["u2_spacer_spacer_bpp_overlap_median"]),
        ],
        (-0.003, 0.105), "Two-repeat spacer signal",
    )
    plot_null_rows(
        ax2, one, obs_one,
        [
            ("All 35", "spacer_spacer_bpp_overlap_median", obs_one["supported_spacer_spacer_bpp_overlap_median"]),
            ("U1 (15)", "u1_spacer_spacer_bpp_overlap_median", obs_one["u1_spacer_spacer_bpp_overlap_median"]),
            ("U2 (18)", "u2_spacer_spacer_bpp_overlap_median", obs_one["u2_spacer_spacer_bpp_overlap_median"]),
        ],
        (-0.003, 0.095), "One-repeat boundary sensitivity",
    )
    plot_null_rows(
        ax3, two, obs_two,
        [("Mixed", "mixed_bpp_overlap_median", obs_two["supported_mixed_bpp_overlap_median"])],
        (-0.004, 0.09), "Repeat–spacer mixed pairing",
    )
    plot_null_rows(
        ax4, two, obs_two,
        [("Repeat", "repeat_repeat_bpp_overlap_median", obs_two["supported_repeat_repeat_bpp_overlap_median"])],
        (0.15, 0.98), "Repeat–repeat pairing",
    )
    for letter, axis in zip(("a", "b", "c", "d"), (ax1, ax2, ax3, ax4)):
        panel_letter(axis, letter)
    ax1.text(0.98, 0.94, "0/64 null panels ≥ observed", transform=ax1.transAxes,
             ha="right", va="top", fontsize=8.0, color=OBSERVED)
    ax2.text(0.98, 0.94, "0/64 null panels ≥ observed", transform=ax2.transAxes,
             ha="right", va="top", fontsize=8.0, color=OBSERVED)
    ax3.text(0.98, 0.94, "0/64", transform=ax3.transAxes,
             ha="right", va="top", fontsize=8.0, color=OBSERVED)
    ax4.text(0.98, 0.94, "20/64", transform=ax4.transAxes,
             ha="right", va="top", fontsize=8.0, color=OBSERVED)
    fig.text(0.50, 0.022,
             "Gray: 64 whole-panel dinucleotide-preserving shuffles; ◆ observed median.\nBlack lines: null range and median.",
             ha="center", va="bottom", fontsize=8.0, color=TEXT, linespacing=1.25)
    return save_figure(fig, figure_dir, "figure_p6_2_dinucleotide_panel_calibration")


def parse_ranges(value):
    ranges = []
    for item in value.split(";"):
        left, right = item.split("-")
        ranges.append((int(left), int(right)))
    return ranges


def figure_representative(data_dir, figure_dir):
    edges = read_tsv(os.path.join(data_dir, "representative_alignment_bpp.tsv"))
    metadata = read_tsv(os.path.join(data_dir, "representative_alignment_metadata.tsv"))
    edges_by_ordinal = defaultdict(list)
    for row in edges:
        edges_by_ordinal[row["ordinal"]].append(row)

    fig = plt.figure(figsize=(170 / 25.4, 100 / 25.4))
    grid = GridSpec(1, 2, left=0.075, right=0.985, bottom=0.16, top=0.82, wspace=0.24)
    for index, meta in enumerate(metadata):
        axis = fig.add_subplot(grid[0, index])
        alignment_length = int(meta["alignment_length"])
        for row in edges_by_ordinal[meta["ordinal"]]:
            left = int(row["left_alignment_column_1based"])
            right = int(row["right_alignment_column_1based"])
            probability = float(row["probability"])
            color = PAIR_COLORS[row["pair_class"]]
            size = 2.5 + 34.0 * math.sqrt(probability)
            if row["side"] == "PALS_2":
                x, y = left, right
            else:
                x, y = right, left
            axis.scatter(x, y, s=size, color=color, alpha=0.55,
                         edgecolor="none", rasterized=False, zorder=2)
        axis.plot([1, alignment_length], [1, alignment_length], color="#6E6E6E", lw=0.65, zorder=1)
        for start, end in parse_ranges(meta["left_repeat_columns_1based"]):
            axis.plot([start, end], [alignment_length + 4, alignment_length + 4],
                      color=PAIR_COLORS["repeat_repeat"], lw=3.0, solid_capstyle="butt", clip_on=False)
        for start, end in parse_ranges(meta["right_repeat_columns_1based"]):
            axis.plot([start, end], [-4, -4], color=PAIR_COLORS["repeat_repeat"],
                      lw=3.0, solid_capstyle="butt", clip_on=False)
        axis.set_xlim(0, alignment_length + 1)
        axis.set_ylim(0, alignment_length + 1)
        axis.set_aspect("equal", adjustable="box")
        axis.set_xlabel("Aligned context column")
        if index == 0:
            axis.set_ylabel("Aligned context column")
        axis.set_title(
            "{}  |  spacer {:.3f}  |  repeat {:.3f}".format(
                meta["ordinal"], float(meta["spacer_spacer_bpp_overlap"]),
                float(meta["repeat_repeat_bpp_overlap"])),
            loc="left", fontweight="bold", fontsize=8.8,
        )
        axis.text(0.03, 0.96, "PALS_2 (upper)", transform=axis.transAxes,
                  ha="left", va="top", fontsize=8.0)
        axis.text(0.97, 0.04, "PB50 (lower)", transform=axis.transAxes,
                  ha="right", va="bottom", fontsize=8.0)
        clean_axis(axis)
        panel_letter(axis, chr(ord("a") + index))

    handles = [Line2D([0], [0], marker="o", lw=0, markersize=5.5,
                      markerfacecolor=PAIR_COLORS[key], markeredgecolor="none", label=label)
               for key, label in [
                   ("spacer_spacer", "spacer–spacer"),
                   ("repeat_repeat", "repeat–repeat"),
                   ("mixed", "mixed"),
               ]]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 0.005), handletextpad=0.4, columnspacing=1.6)
    fig.text(0.50, 0.965,
             "Representative supported PALS_2–PB50 ensemble maps (pair probability ≥ 0.01)",
             ha="center", va="top", fontweight="bold", fontsize=9.5)
    return save_figure(fig, figure_dir, "figure_p6_3_representative_pair_maps")


def main():
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    data_dir = os.path.join(root, "data", "processed", "p6_rna")
    figure_dir = os.path.join(root, "figures", "p6_rna")
    configure()
    outputs = []
    outputs.extend(figure_specificity(data_dir, figure_dir))
    outputs.extend(figure_null(data_dir, figure_dir))
    outputs.extend(figure_representative(data_dir, figure_dir))
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
