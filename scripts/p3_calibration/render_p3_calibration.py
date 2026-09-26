#!/usr/bin/env python3
"""Render the single P3 panel-level calibration figure."""

from __future__ import division

import csv
import json
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


NULL_COLOR = "#7F8C99"
OBSERVED_COLOR = "#B6493A"
MEDIAN_COLOR = "#242424"
TEXT_COLOR = "#252525"


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


def jitter(index, amplitude=0.13):
    return amplitude * math.sin((index + 1) * 2.399963229728653)


def style_axis(axis):
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_visible(False)
    axis.tick_params(axis="y", length=0)


def plot_single_row(axis, values, observed, median, y=0.0):
    for index, value in enumerate(values):
        axis.scatter(value, y + jitter(index), s=17, color=NULL_COLOR,
                     alpha=0.72, edgecolor="white", linewidth=0.25, zorder=2)
    axis.plot([min(values), max(values)], [y - 0.25, y - 0.25],
              color=MEDIAN_COLOR, lw=1.0, zorder=3)
    axis.plot([median, median], [y - 0.33, y - 0.17],
              color=MEDIAN_COLOR, lw=1.3, zorder=3)
    axis.scatter(observed, y + 0.30, marker="D", s=42, color=OBSERVED_COLOR,
                 edgecolor="white", linewidth=0.45, zorder=4)


def main():
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    data_dir = os.path.join(root, "data", "processed", "p3_calibration")
    figure_dir = os.path.join(root, "figures", "p3")
    os.makedirs(figure_dir, exist_ok=True)
    rows = read_tsv(os.path.join(data_dir, "null_replicates.tsv"))
    with open(os.path.join(data_dir, "p3_calibration_summary.json"), "r", encoding="utf-8") as handle:
        summary = json.load(handle)
    observed = summary["observed_reused_from_P2"]
    distributions = summary["null_distributions"]

    counts = [int(row["supported_pair_count"]) for row in rows]
    u1 = [int(row["u1_same_position_locus_pair_coverage"]) for row in rows]
    u2 = [int(row["u2_same_position_locus_pair_coverage"]) for row in rows]
    identities = [float(row["supported_trimmed_spacer_identity_median"]) for row in rows]

    configure()
    fig = plt.figure(figsize=(175 / 25.4, 105 / 25.4))
    grid = GridSpec(2, 2, height_ratios=[0.95, 1.05], width_ratios=[1.02, 0.98],
                    left=0.08, right=0.97, bottom=0.16, top=0.92,
                    hspace=0.52, wspace=0.32)
    ax1 = fig.add_subplot(grid[0, :])
    ax2 = fig.add_subplot(grid[1, 0])
    ax3 = fig.add_subplot(grid[1, 1])

    count_dist = distributions["supported_pair_count"]
    plot_single_row(ax1, counts, observed["supported_pair_count"], count_dist["median"])
    ax1.set_xlim(0, 37)
    ax1.set_ylim(-0.48, 0.55)
    ax1.set_yticks([])
    ax1.set_xlabel("Deduplicated cross-group supported unit pairs per full panel")
    ax1.set_title("a  Whole-panel null under fixed 47-copy / 37-unit geometry",
                  loc="left", fontweight="bold")
    ax1.text(0.01, 0.91, "Null: median 8, range 4–14, p95 12",
             transform=ax1.transAxes, ha="left", va="top", color=TEXT_COLOR)
    ax1.text(0.99, 0.91, "Observed 35  |  0/64 null panels reached it",
             transform=ax1.transAxes, ha="right", va="top", color=OBSERVED_COLOR,
             fontweight="bold")
    style_axis(ax1)

    for row_index, (label, values, observed_value, median_value, range_text) in enumerate([
        ("U1 ↔ U1", u1, observed["u1_same_position_locus_pair_coverage"],
         distributions["u1_same_position_locus_pair_coverage"]["median"], "null 1–7"),
        ("U2 ↔ U2", u2, observed["u2_same_position_locus_pair_coverage"],
         distributions["u2_same_position_locus_pair_coverage"]["median"], "null 0–4"),
    ]):
        y = 1 - row_index
        for index, value in enumerate(values):
            ax2.scatter(value, y + jitter(index, 0.10), s=16, color=NULL_COLOR,
                        alpha=0.70, edgecolor="white", linewidth=0.25, zorder=2)
        ax2.plot([min(values), max(values)], [y - 0.20, y - 0.20], color=MEDIAN_COLOR, lw=0.9)
        ax2.plot([median_value, median_value], [y - 0.27, y - 0.13], color=MEDIAN_COLOR, lw=1.2)
        ax2.scatter(observed_value, y, marker="D", s=40, color=OBSERVED_COLOR,
                    edgecolor="white", linewidth=0.45, zorder=4)
        ax2.text(20.8, y + 0.18, "observed {} / 21".format(observed_value),
                 ha="right", va="bottom", fontsize=7.2, color=OBSERVED_COLOR,
                 fontweight="bold")
        ax2.text(20.8, y - 0.18, range_text, ha="right", va="top",
                 fontsize=6.8, color=TEXT_COLOR)
    ax2.set_xlim(-0.5, 21.5)
    ax2.set_ylim(-0.48, 1.48)
    ax2.set_yticks([0, 1])
    ax2.set_yticklabels(["U2 ↔ U2", "U1 ↔ U1"])
    ax2.set_xticks([0, 5, 10, 15, 20])
    ax2.set_xlabel("Cross-group locus-pair coverage (of 21)")
    ax2.set_title("b  Distal positional coverage", loc="left", fontweight="bold")
    style_axis(ax2)

    identity_dist = distributions["supported_trimmed_spacer_identity_median"]
    plot_single_row(
        ax3,
        identities,
        observed["supported_trimmed_spacer_identity_median"],
        identity_dist["median"],
    )
    ax3.set_xlim(0.465, 0.568)
    ax3.set_ylim(-0.48, 0.55)
    ax3.set_yticks([])
    ax3.set_xlabel("Median trimmed-spacer identity")
    ax3.set_title("c  Sequence identity after selection", loc="left", fontweight="bold")
    ax3.text(0.02, 0.91, "Null median 0.504\nrange 0.478–0.526",
             transform=ax3.transAxes, ha="left", va="top", color=TEXT_COLOR)
    ax3.text(0.98, 0.91, "Observed 0.557",
             transform=ax3.transAxes, ha="right", va="top", color=OBSERVED_COLOR,
             fontweight="bold")
    style_axis(ax3)

    fig.text(0.5, 0.025,
             "Each gray circle: one shared ten-locus null panel   ◆ observed P2 chain-only result   — null range with median tick",
             ha="center", va="bottom", fontsize=7.0, color=TEXT_COLOR)

    outputs = []
    stem = "figure4_spacer_order_calibration"
    for extension in ("pdf", "svg", "png"):
        path = os.path.join(figure_dir, stem + "." + extension)
        kwargs = {"facecolor": "white"}
        if extension == "png":
            kwargs["dpi"] = 300
        fig.savefig(path, **kwargs)
        outputs.append(path)
    plt.close(fig)
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
