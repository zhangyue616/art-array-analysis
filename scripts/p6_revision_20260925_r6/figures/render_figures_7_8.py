#!/usr/bin/env python3
"""Render revised Figures 7-8 from accepted P6 RNA-seq aggregate tables."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from common import ROOT, TEXT, UNIT_COLORS, clean_axis, configure, panel_label, panel_title, read_tsv, save_figure


DATA_DIR = ROOT / "data" / "processed" / "p6_rnaseq" / "aggregate"
UNITS = ("U1", "U2", "U3", "U4")
CULTURES = (1, 2, 3)
POST_TIMES = (5, 15, 55)
CULTURE_MARKERS = {1: "o", 2: "s", 3: "^"}
CULTURE_STYLES = {1: "-", 2: "--", 3: ":"}


def _median(values: list[float]) -> float:
    values = sorted(values)
    mid = len(values) // 2
    return values[mid] if len(values) % 2 else 0.5 * (values[mid - 1] + values[mid])


def _unit_panel_title(axis: plt.Axes, letter: str, unit: str) -> None:
    """Use the shared title geometry while retaining the unit color identity."""
    first_new_text = len(axis.texts)
    panel_title(axis, letter, unit, pad=8)
    axis.texts[first_new_text].set_color(UNIT_COLORS[unit])


def render_figure_7() -> list[Path]:
    unit_rows = read_tsv(DATA_DIR / "unit_metrics.tsv")
    paired_rows = read_tsv(DATA_DIR / "paired_contrasts.tsv")
    unit_index = {
        (int(row["time_min"]), int(row["culture"]), row["unit"]): row
        for row in unit_rows
    }
    primary = {
        (row["unit"], int(row["culture"])): row
        for row in paired_rows
        if row["comparison"] == "5_to_55" and row["endpoint"] == "spacer_midpoint_density_share"
    }

    fig = plt.figure(figsize=(170 / 25.4, 145 / 25.4))
    grid = GridSpec(
        2,
        4,
        figure=fig,
        height_ratios=(1.0, 0.86),
        left=0.09,
        right=0.98,
        bottom=0.14,
        top=0.94,
        hspace=0.46,
        wspace=0.28,
    )
    top_axes = [fig.add_subplot(grid[0, index]) for index in range(4)]
    delta_axis = fig.add_subplot(grid[1, :])
    # Use elapsed minutes as the actual x coordinate so the 5->15 and 15->55
    # intervals retain their 10:40 spacing rather than appearing equidistant.
    x_time = list(POST_TIMES)

    for panel_index, (axis, unit) in enumerate(zip(top_axes, UNITS)):
        for culture in CULTURES:
            values = [float(unit_index[(time, culture, unit)]["spacer_midpoint_density_share"]) for time in POST_TIMES]
            axis.plot(
                x_time,
                values,
                color=UNIT_COLORS[unit],
                marker=CULTURE_MARKERS[culture],
                linestyle=CULTURE_STYLES[culture],
                markersize=4.2,
                markeredgecolor="white",
                markeredgewidth=0.45,
                alpha=0.90,
            )
        axis.set_xlim(0, 60)
        axis.set_ylim(0, 1)
        axis.set_xticks(x_time)
        axis.set_xticklabels([str(value) for value in POST_TIMES])
        axis.set_xlabel("Minutes")
        if panel_index == 0:
            axis.set_ylabel("Absolute density share (0–1)")
        else:
            axis.set_yticklabels([])
        clean_axis(axis, "y")
        _unit_panel_title(axis, chr(ord("a") + panel_index), unit)

    x_units = list(range(len(UNITS)))
    jitter = {1: -0.12, 2: 0.0, 3: 0.12}
    for unit_index_value, unit in enumerate(UNITS):
        values_pp: list[float] = []
        for culture in CULTURES:
            row = primary[(unit, culture)]
            value_pp = 100.0 * float(row["delta"])
            values_pp.append(value_pp)
            delta_axis.scatter(
                [unit_index_value + jitter[culture]],
                [value_pp],
                s=36,
                marker=CULTURE_MARKERS[culture],
                facecolor=UNIT_COLORS[unit],
                edgecolor="white",
                linewidth=0.55,
                zorder=3,
            )
        med = _median(values_pp)
        delta_axis.plot([unit_index_value - 0.22, unit_index_value + 0.22], [med, med], color="#111111", lw=1.2, zorder=4)
    delta_axis.axhline(0, color="#777777", lw=0.8)
    for threshold in (-5, 5):
        delta_axis.axhline(threshold, color="#8A8A8A", lw=0.8, ls="--", zorder=0)
    delta_axis.set_xlim(-0.5, 3.5)
    # Preserve a symmetric blank band beyond the +/-5-point reference lines so
    # the rule annotation never competes with observations or the threshold.
    delta_axis.set_ylim(-7.0, 7.0)
    delta_axis.set_xticks(x_units)
    delta_axis.set_xticklabels(UNITS)
    delta_axis.set_ylabel("5→55 min change (percentage points)")
    panel_title(delta_axis, "e", "Paired change by culture", pad=8)
    delta_axis.text(
        0.99,
        0.96,
        "Prespecified magnitude scale: |Δ| = 5 percentage points",
        transform=delta_axis.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        fontweight="normal",
        color=TEXT,
        bbox=dict(facecolor="white", edgecolor="none", pad=1.0),
        zorder=5,
    )
    clean_axis(delta_axis, "y")

    culture_handles = [
        Line2D(
            [0],
            [0],
            marker=CULTURE_MARKERS[culture],
            color="#555555",
            linestyle=CULTURE_STYLES[culture],
            linewidth=1.0,
            markersize=5.5,
            label=f"Culture {culture}",
        )
        for culture in CULTURES
    ]
    fig.legend(handles=culture_handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.012))
    return save_figure(fig, "Figure_7")


def render_figure_8() -> list[Path]:
    paired_rows = read_tsv(DATA_DIR / "paired_contrasts.tsv")
    partition_rows = read_tsv(DATA_DIR / "whole_array_partitions.tsv")
    paired_index = {
        (row["endpoint"], row["unit"], int(row["culture"])): row
        for row in paired_rows
        if row["comparison"] == "5_to_55"
    }

    confirm_points: list[tuple[str, int, float, float]] = []
    sensitivity_points: list[tuple[str, int, float, float]] = []
    for unit in UNITS:
        for culture in CULTURES:
            primary = paired_index[("spacer_midpoint_density_share", unit, culture)]
            confirm = paired_index[("spacer_depth_share", unit, culture)]
            core = paired_index[("core_midpoint_density_share", unit, culture)]
            if primary["pair_status"] == "complete" and confirm["pair_status"] == "complete" and primary["both_libraries_estimable"] == "true" and confirm["both_libraries_estimable"] == "true":
                confirm_points.append((unit, culture, float(primary["delta"]), float(confirm["delta"])))
            if primary["pair_status"] == "complete" and core["pair_status"] == "complete" and primary["both_libraries_estimable"] == "true" and core["both_libraries_estimable"] == "true":
                sensitivity_points.append((unit, culture, float(primary["delta"]), float(core["delta"])))

    fig = plt.figure(figsize=(170 / 25.4, 145 / 25.4))
    grid = GridSpec(2, 2, figure=fig, height_ratios=(1.0, 1.08), left=0.08, right=0.98, bottom=0.24, top=0.90, hspace=0.54, wspace=0.30)
    ax_confirm = fig.add_subplot(grid[0, 0])
    ax_sensitivity = fig.add_subplot(grid[0, 1])
    ax_partition = fig.add_subplot(grid[1, :])

    for axis, points in ((ax_confirm, confirm_points), (ax_sensitivity, sensitivity_points)):
        for unit, culture, x_value, y_value in points:
            axis.scatter(
                [x_value],
                [y_value],
                color=UNIT_COLORS[unit],
                marker=CULTURE_MARKERS[culture],
                s=27,
                linewidths=0.45,
                edgecolors="white",
                zorder=3,
            )
    all_deltas = [abs(value) for _, _, x_value, y_value in confirm_points + sensitivity_points for value in (x_value, y_value)]
    limit = max(0.05, max(all_deltas) * 1.12)
    for axis, point_count in ((ax_confirm, len(confirm_points)), (ax_sensitivity, len(sensitivity_points))):
        axis.plot([-limit, limit], [-limit, limit], color="#777777", linewidth=0.8, linestyle="--", zorder=0)
        axis.axhline(0, color="#BDBDBD", linewidth=0.6, zorder=0)
        axis.axvline(0, color="#BDBDBD", linewidth=0.6, zorder=0)
        axis.set_xlim(-limit, limit)
        axis.set_ylim(-limit, limit)
        axis.set_aspect("equal", adjustable="box")
        clean_axis(axis, "both")
        axis.text(0.02, 0.98, f"evaluable pairs: {point_count}/12", transform=axis.transAxes, ha="left", va="top", fontsize=8)
        axis.text(
            0.96,
            0.91,
            "y = x",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=8,
            color="#666666",
            bbox=dict(facecolor="white", edgecolor="none", pad=0.5),
        )
    ax_confirm.set_xlabel("Δ spacer density share")
    ax_confirm.set_ylabel("Δ spacer depth share")
    ax_sensitivity.set_xlabel("Δ spacer density share")
    ax_sensitivity.set_ylabel("Δ core-to-core density share")
    panel_title(ax_confirm, "a", "Density vs fragment depth", pad=8)
    panel_title(ax_sensitivity, "b", "Density vs core-to-core", pad=8)
    unit_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=UNIT_COLORS[unit], markeredgecolor="white", markersize=5.5, label=unit)
        for unit in UNITS
    ]
    culture_handles = [
        Line2D([0], [0], marker=CULTURE_MARKERS[culture], color="#4D4D4D", linestyle="None", markersize=5.5, label=f"Culture {culture}")
        for culture in CULTURES
    ]
    fig.text(0.10, 0.968, "Unit", ha="left", va="center", fontsize=8.5, color=TEXT)
    fig.legend(
        handles=unit_handles,
        loc="upper left",
        bbox_to_anchor=(0.145, 0.995),
        ncol=4,
        frameon=False,
        columnspacing=0.75,
        handletextpad=0.35,
    )
    fig.text(0.545, 0.968, "Culture", ha="left", va="center", fontsize=8.5, color=TEXT)
    fig.legend(
        handles=culture_handles,
        loc="upper left",
        bbox_to_anchor=(0.601, 0.995),
        ncol=3,
        frameon=False,
        columnspacing=0.75,
        handletextpad=0.35,
    )

    rows_by_run: dict[tuple[int, int, str], dict[str, dict[str, str]]] = defaultdict(dict)
    for row in partition_rows:
        key = (int(row["time_min"]), int(row["culture"]), row["run_accession"])
        rows_by_run[key][row["partition"]] = row
    ordered_keys = sorted(rows_by_run)
    x = list(range(len(ordered_keys)))
    bottoms = [0.0] * len(ordered_keys)
    partitions = ("distal_residual", "U1", "U2", "U3", "U4", "terminal_residual")
    partition_colours = {
        "distal_residual": "#BDBDBD",
        "U1": UNIT_COLORS["U1"],
        "U2": UNIT_COLORS["U2"],
        "U3": UNIT_COLORS["U3"],
        "U4": UNIT_COLORS["U4"],
        "terminal_residual": "#666666",
    }
    partition_hatches = {"distal_residual": "//", "terminal_residual": "\\\\"}
    partition_labels = {
        "distal_residual": "distal residual",
        "U1": "U1",
        "U2": "U2",
        "U3": "U3",
        "U4": "U4",
        "terminal_residual": "terminal residual (RT-adjacent)",
    }
    for partition in partitions:
        values = [float(rows_by_run[key][partition]["whole_array_midpoint_share"]) for key in ordered_keys]
        container = ax_partition.bar(
            x,
            values,
            bottom=bottoms,
            width=0.72,
            color=partition_colours[partition],
            edgecolor="white",
            linewidth=0.35,
            hatch=partition_hatches.get(partition),
            label=partition_labels[partition],
        )
        for rect, key in zip(container.patches, ordered_keys):
            if key[0] == 0:
                rect.set_alpha(0.48)
        bottoms = [bottom + value for bottom, value in zip(bottoms, values)]

    labels = [f"{time}-{culture}" for time, culture, _ in ordered_keys]
    ax_partition.set_xticks(x)
    ax_partition.set_xticklabels(labels, rotation=45, ha="right")
    ax_partition.set_ylim(0, 1)
    ax_partition.set_ylabel("Whole-array midpoint share")
    ax_partition.set_xlabel("Time (min) – culture")
    ax_partition.axvline(2.5, color="#4D4D4D", linewidth=0.8, linestyle="--")
    ax_partition.text(1, 1.14, "0 min: low-count background", ha="center", va="bottom", fontsize=8, clip_on=False)
    ax_partition.text(7, 1.14, "post-infection", ha="center", va="bottom", fontsize=8, clip_on=False)
    whole_counts: list[int] = []
    for key in ordered_keys:
        total = sum(int(rows_by_run[key][partition]["sense_midpoint_fragments"]) for partition in partitions)
        whole_counts.append(total)
    for index, key in enumerate(ordered_keys[:3]):
        ax_partition.text(index, 1.025, f"n={whole_counts[index]}", ha="center", va="bottom", fontsize=8, color=TEXT, clip_on=False)
    ax_partition.text(0.99, 1.02, "Distal residual: 0–2 midpoint fragments per library", transform=ax_partition.transAxes, ha="right", va="bottom", fontsize=8, color=TEXT, clip_on=False)
    clean_axis(ax_partition, "y")
    panel_label(ax_partition, "c", x=-0.075, y=1.22)
    # Build independent, fully opaque partition handles: the first plotted bars
    # are the reduced-opacity 0-min samples and must not fade the partition key.
    partition_handles = [
        Patch(
            facecolor=partition_colours[partition],
            edgecolor="white",
            linewidth=0.35,
            hatch=partition_hatches.get(partition),
            label=partition_labels[partition],
        )
        for partition in partitions
    ]
    baseline_handle = Patch(
        facecolor="#777777",
        edgecolor="#555555",
        linewidth=0.35,
        alpha=0.48,
        label="0 min (reduced opacity): low-count background",
    )
    # Matplotlib fills legend columns top-to-bottom. This ordering displays the
    # six components left-to-right in their actual array partition order.
    partition_legend_order = [0, 3, 1, 4, 2, 5]
    partition_legend = fig.legend(
        handles=[partition_handles[index] for index in partition_legend_order],
        title="Whole-array partitions (U1–U4 core-to-core)",
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.045),
        columnspacing=1.05,
        handletextpad=0.45,
    )
    partition_legend.get_title().set_fontsize(8)
    fig.legend(
        handles=[baseline_handle],
        loc="lower center",
        ncol=1,
        frameon=False,
        bbox_to_anchor=(0.5, 0.005),
    )
    return save_figure(fig, "Figure_8")


def main() -> None:
    configure()
    outputs: list[Path] = []
    outputs.extend(render_figure_7())
    outputs.extend(render_figure_8())
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
