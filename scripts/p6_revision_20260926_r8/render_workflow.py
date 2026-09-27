#!/usr/bin/env python3
"""Render Figure 1, the object-led analysis workflow, from saved tables.

The renderer reads only saved, derived tables.  It does not rerun sequence
matching, RNA folding, null generation, read alignment, or RNA-seq aggregation.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter
from pathlib import Path

# Prefer an optional repository-local dependency directory when present. Keep
# the change process-local and leave global Python configuration untouched.
PROJECT_ROOT = Path(__file__).resolve().parents[2]
LOCAL_PYTHON_PACKAGES = PROJECT_ROOT / "tools" / "p6_rnaseq" / "python_packages"
if LOCAL_PYTHON_PACKAGES.exists():
    sys.path.insert(0, str(LOCAL_PYTHON_PACKAGES))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle


WIDTH_MM = 170.0
HEIGHT_MM = 115.0

INK = "#1F2C36"
MUTED = "#5D6A73"
FAINT = "#D9E0E5"
PALE = "#EEF2F4"
SEQUENCE = "#176CA4"
PAIRING = "#20786A"
TIME = "#A95B20"
UNIT_COLORS = {
    "U1": "#0072B2",
    "U2": "#009E73",
    "U3": "#CC79A7",
    "U4": "#D55E00",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load_saved_objects(root: Path) -> dict[str, object]:
    data_root = root / "data" / "processed"
    object_rows = read_tsv(data_root / "p6_rna" / "unit_sequence_objects.tsv")
    primary = [row for row in object_rows if row["context"] == "one_repeat"]
    require(len({row["unit_id"] for row in primary}) == 37,
            "Expected 37 unique primary one-repeat units")

    group_counts = Counter(row["analysis_subgroup"] for row in primary)
    require(group_counts == Counter({"P0_seed7": 28, "PB50_near3": 9}),
            f"Unexpected primary group sizes: {group_counts}")
    all_cross_group = group_counts["P0_seed7"] * group_counts["PB50_near3"]
    require(all_cross_group == 252, "Expected 28 x 9 = 252 cross-group combinations")

    correspondence_rows = read_tsv(
        data_root / "p2_comparison" / "supported_unit_correspondences.tsv"
    )
    selected = [row for row in correspondence_rows if row["relationship"] == "cross_group"]
    require(len(selected) == 35, "Expected 35 selected cross-group correspondences")
    other = all_cross_group - len(selected)
    require(other == 217, "Expected 217 non-selected cross-group combinations")

    representative = [
        row for row in selected
        if {row["left_label"], row["right_label"]} == {"PALS_2", "PB50"}
    ]
    require(len(representative) == 2, "Expected two PALS_2-PB50 correspondences")
    representative_ordinals = {
        (int(row["left_ordinal"]), int(row["right_ordinal"])) for row in representative
    }
    require(representative_ordinals == {(1, 1), (2, 2)},
            f"Unexpected PALS_2-PB50 ordinal map: {representative_ordinals}")

    tracks: dict[str, list[dict[str, str]]] = {}
    for label in ("PALS_2", "PB50"):
        rows = sorted(
            [row for row in primary if row["label"] == label],
            key=lambda row: int(row["unit_ordinal_distal_to_proximal"]),
        )
        require(rows, f"Missing primary unit objects for {label}")
        tracks[label] = rows

    metadata_rows = read_tsv(data_root / "p6_rna" / "representative_alignment_metadata.tsv")
    metadata = {row["ordinal"]: row for row in metadata_rows}
    require(set(metadata) == {"U1", "U2"},
            f"Unexpected representative alignment metadata: {set(metadata)}")
    bpp_rows = read_tsv(data_root / "p6_rna" / "representative_alignment_bpp.tsv")
    bpp_u1 = [row for row in bpp_rows if row["ordinal"] == "U1"]
    require(len(bpp_u1) > 0, "No saved U1 BPP edges")
    require({row["side"] for row in bpp_u1} == {"PALS_2", "PB50"},
            "Expected saved U1 BPP edges from both representative units")

    window_rows = read_tsv(data_root / "p6_rnaseq" / "sa1_windows.tsv")
    window_by_id = {row["window_id"]: row for row in window_rows}
    needed_windows = {
        "array_whole", "array_distal_residual", "array_terminal_residual",
        "U1_spacer", "U2_spacer", "U3_spacer", "U4_spacer",
        "repeat01_mask", "repeat02_mask", "repeat03_mask",
        "repeat04_mask", "repeat05_mask",
    }
    require(needed_windows <= set(window_by_id),
            f"Missing SA1 windows: {sorted(needed_windows - set(window_by_id))}")

    temporal_rows = read_tsv(
        data_root / "p6_revision_20260925" / "s2_temporal_absolute_relative_changes.tsv"
    )
    require({row["unit"] for row in temporal_rows} == set(UNIT_COLORS),
            "Unexpected units in saved temporal-change table")
    require({int(row["culture"]) for row in temporal_rows} == {1, 2, 3},
            "Expected three cultures in saved temporal-change table")

    library_receipt = None
    library_path = data_root / "p6_rnaseq" / "aggregate" / "library_metrics.tsv"
    if library_path.exists():
        library_rows = read_tsv(library_path)
        require(len(library_rows) == 12, "Expected 12 saved SA1 library rows")
        require({int(row["time_min"]) for row in library_rows} == {0, 5, 15, 55},
                "Expected SA1 time points 0, 5, 15, and 55 min")
        require({int(row["culture"]) for row in library_rows} == {1, 2, 3},
                "Expected three cultures in saved SA1 library rows")
        background_rows = [row for row in library_rows if int(row["time_min"]) == 0]
        background_counts = [int(row["spacer_sense_midpoint_fragments"])
                             for row in background_rows]
        require(all(row["unit_estimable_ge100"].lower() == "false"
                    for row in background_rows),
                "Expected all 0-min libraries to be below the unit-estimability threshold")
        library_receipt = {
            "rows": len(library_rows),
            "background_count_min": min(background_counts),
            "background_count_max": max(background_counts),
        }

    return {
        "tracks": tracks,
        "selected_count": len(selected),
        "other_count": other,
        "all_cross_group": all_cross_group,
        "bpp_u1": bpp_u1,
        "metadata": metadata,
        "windows": window_by_id,
        "library_receipt": library_receipt,
    }


def configure_matplotlib() -> None:
    matplotlib.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
        "axes.linewidth": 0.6,
    })


def add_text(ax, x: float, y: float, label: str, *, size: float = 8.3,
             color: str = INK, weight: str = "normal", ha: str = "left",
             va: str = "top", linespacing: float = 1.1, zorder: int = 10):
    require(size >= 8.0, f"Figure text below 8 pt: {label!r} at {size} pt")
    return ax.text(x, y, label, fontsize=size, color=color, weight=weight,
                   ha=ha, va=va, linespacing=linespacing, zorder=zorder)


def add_arrow(ax, start: tuple[float, float], end: tuple[float, float], *,
              color: str = MUTED, width: float = 0.7, style: str = "-|>",
              linestyle: str = "-") -> None:
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle=style, mutation_scale=7.0,
        linewidth=width, color=color, linestyle=linestyle,
        shrinkA=0, shrinkB=0, zorder=4,
    ))


def gene_arrow(ax, x0: float, x1: float, y: float, height: float,
               color: str, *, edgecolor: str | None = None) -> None:
    head = min(2.4, (x1 - x0) * 0.22)
    points = [
        (x0, y - height / 2), (x1 - head, y - height / 2),
        (x1 - head, y - height), (x1, y),
        (x1 - head, y + height), (x1 - head, y + height / 2),
        (x0, y + height / 2),
    ]
    ax.add_patch(Polygon(points, closed=True, facecolor=color,
                         edgecolor=edgecolor or color, linewidth=0.5, zorder=3))


def panel_heading(ax, x0: float, x1: float, y: float, letter: str,
                  title: str, color: str, note: str | None = None) -> None:
    ax.plot([x0, x1], [y, y], color=color, lw=1.15, solid_capstyle="butt", zorder=2)
    add_text(ax, x0, y + 1.5, letter, size=10.2, weight="bold", color=INK)
    add_text(ax, x0 + 5.1, y + 1.5, title, size=9.8, weight="bold", color=color)
    if note:
        add_text(ax, x1, y + 1.8, note, size=8.0, color=MUTED, ha="right")


def draw_shared_locus(ax) -> None:
    add_text(ax, 4.5, 3.2, "Shared ART-array objects", size=10.2,
             weight="bold", color=INK)
    add_text(ax, 109.5, 3.7, "10 loci · 37 primary units", size=8.8,
             weight="bold", color=SEQUENCE, ha="right")

    add_text(ax, 15.6, 10.1, "coding anchor", size=8.0, color=MUTED, ha="center")
    add_text(ax, 49.4, 10.1, "repeat / spacer array", size=8.0,
             color=MUTED, ha="center")
    add_text(ax, 84.8, 10.1, "RT", size=8.0, color=MUTED, ha="center")
    add_text(ax, 102.0, 10.1, "partner", size=8.0, color=MUTED, ha="center")
    gene_arrow(ax, 7.0, 25.2, 17.1, 3.7, "#9DA7AE")

    repeat_x = [29.0, 38.3, 48.2, 59.2, 71.2]
    spacer_widths = [8.0, 8.6, 9.7, 10.7]
    for idx, x in enumerate(repeat_x):
        ax.add_patch(Rectangle((x - 0.85, 14.4), 1.7, 5.4,
                               facecolor="#515C64", edgecolor="none", zorder=4))
        if idx < len(spacer_widths):
            x_next = repeat_x[idx + 1]
            ax.add_patch(Rectangle((x + 0.85, 15.2), x_next - x - 1.7, 3.8,
                                   facecolor="#DDE8ED", edgecolor="#91A5AF",
                                   linewidth=0.45, zorder=3))
    add_arrow(ax, (27.0, 17.1), (73.0, 17.1), color="#87949C",
              width=0.45, style="->")
    gene_arrow(ax, 77.0, 94.2, 17.1, 3.7, "#3F4A52")
    gene_arrow(ax, 96.3, 109.0, 17.1, 3.7, "#7D8991")
    add_text(ax, 58.8, 23.1, "Locus schematic",
             size=8.0, color=MUTED, ha="center")


def unit_x(position: int, x0: float, x1: float) -> float:
    low, high = -1120.0, -285.0
    return x0 + (position - low) / (high - low) * (x1 - x0)


def draw_unit_track(ax, rows: list[dict[str, str]], label: str,
                    y: float, x0: float = 27.0, x1: float = 76.0) -> dict[int, tuple[float, float]]:
    add_text(ax, x0 - 2.2, y, label, size=8.2, weight="bold", ha="right", va="center")
    ax.plot([x0, x1], [y, y], color="#9CA7AE", lw=0.6, zorder=1)
    boundaries: dict[int, tuple[float, float]] = {}
    for row in rows:
        ordinal = int(row["unit_ordinal_distal_to_proximal"])
        left = unit_x(int(row["left_calibrated_start_relative_to_RT"]), x0, x1)
        right = unit_x(int(row["right_calibrated_start_relative_to_RT"]), x0, x1)
        unit = f"U{ordinal}"
        color = UNIT_COLORS[unit]
        ax.add_patch(Rectangle((left, y - 2.05), right - left, 4.1,
                               facecolor=color, edgecolor="white", linewidth=0.45,
                               alpha=0.90, zorder=3))
        add_text(ax, (left + right) / 2, y + 0.05, unit, size=8.0,
                 weight="bold", color="white", ha="center", va="center")
        ax.add_patch(Rectangle((left - 0.55, y - 2.65), 1.1, 5.3,
                               facecolor="#424D54", edgecolor="none", zorder=5))
        boundaries[ordinal] = (left, right)
    final_right = max(value[1] for value in boundaries.values())
    ax.add_patch(Rectangle((final_right - 0.55, y - 2.65), 1.1, 5.3,
                           facecolor="#424D54", edgecolor="none", zorder=5))
    return boundaries


def draw_correspondence_panel(ax, tracks: dict[str, list[dict[str, str]]]) -> None:
    panel_heading(ax, 4.5, 110.0, 31.5, "a", "Unit correspondence",
                  SEQUENCE)
    pals = draw_unit_track(ax, tracks["PALS_2"], "PALS_2", 45.0)
    pb50 = draw_unit_track(ax, tracks["PB50"], "PB50", 54.5)

    for ordinal in (1, 2):
        color = UNIT_COLORS[f"U{ordinal}"]
        for p_x, b_x in zip(pals[ordinal], pb50[ordinal]):
            ax.plot([p_x, b_x], [47.1, 52.4], color=color, lw=0.85,
                    alpha=0.75, zorder=2)
    add_text(ax, 51.5, 58.4, "sequence correspondence + both mapped boundaries",
             size=8.0, color=MUTED, ha="center")
    add_arrow(ax, (77.8, 49.5), (98.4, 49.5), color=SEQUENCE, width=0.8)
    add_text(ax, 103.7, 43.2, "35 pairs", size=9.8, weight="bold",
             color=SEQUENCE, ha="center")
    add_text(ax, 103.7, 49.4, "overall", size=8.0,
             color=MUTED, ha="center")


def draw_bpp_map(ax, edges: list[dict[str, str]], meta: dict[str, str],
                 x0: float, y0: float, size: float) -> None:
    length = int(meta["alignment_length"])
    ax.add_patch(Rectangle((x0, y0), size, size, facecolor="#F7F9FA",
                           edgecolor="#97A4AC", linewidth=0.55, zorder=1))
    ax.plot([x0, x0 + size], [y0, y0 + size], color="#B0B9BF", lw=0.5, zorder=2)
    plotted = []
    for row in edges:
        left = int(row["left_alignment_column_1based"])
        right = int(row["right_alignment_column_1based"])
        if row["side"] == "PALS_2":
            x_pos, y_pos = left, right
        else:
            x_pos, y_pos = right, left
        x = x0 + (x_pos - 1) / (length - 1) * size
        y = y0 + (y_pos - 1) / (length - 1) * size
        probability = float(row["probability"])
        plotted.append((probability, x, y))
    plotted.sort()
    probabilities = [item[0] for item in plotted]
    xs = [item[1] for item in plotted]
    ys = [item[2] for item in plotted]
    point_sizes = [1.0 + 8.5 * math.sqrt(value) for value in probabilities]
    ax.scatter(xs, ys, s=point_sizes, color=PAIRING, alpha=0.40,
               edgecolors="none", zorder=3, rasterized=False)
    add_text(ax, x0 + size / 2, y0 - 2.9, "PALS_2 / PB50 · U1 BPP", size=8.0,
             color=PAIRING, weight="bold", ha="center", va="bottom")


def draw_sheet_stack(ax, x: float, y: float) -> None:
    for offset in (1.8, 0.9, 0.0):
        ax.add_patch(Rectangle((x + offset, y - offset), 8.0, 5.0,
                               facecolor="white", edgecolor=PAIRING,
                               linewidth=0.55, zorder=3))
        for row in range(2):
            ax.plot([x + offset + 1.2, x + offset + 6.8],
                    [y - offset + 1.6 + row * 1.3] * 2,
                    color="#A8B3B9", lw=0.4, zorder=4)


def draw_pairing_panel(ax, data: dict[str, object]) -> None:
    panel_heading(ax, 4.5, 110.0, 66.5, "b", "Pairing diagnostic",
                  PAIRING)
    add_text(ax, 17.0, 78.0, "28 × 9 = 252", size=10.0,
             weight="bold", color=INK)
    ax.plot([18.0, 20.5], [88.8, 88.8], color=SEQUENCE, lw=2.6,
            solid_capstyle="butt")
    add_text(ax, 22.0, 88.8, f"{data['selected_count']} selected",
             size=8.1, weight="bold", va="center")
    ax.plot([18.0, 20.5], [94.7, 94.7], color="#AAB4BA", lw=2.6,
            solid_capstyle="butt")
    add_text(ax, 22.0, 94.7, f"{data['other_count']} other",
             size=8.1, va="center")

    ax.plot([41.0, 44.0], [88.8, 88.8], color=SEQUENCE, lw=0.75)
    ax.plot([41.0, 44.0], [94.7, 94.7], color="#8B979E", lw=0.75)
    ax.plot([44.0, 44.0], [88.8, 94.7], color="#8B979E", lw=0.75)
    add_arrow(ax, (44.0, 91.75), (46.3, 91.75), color=PAIRING, width=0.8)

    meta = data["metadata"]["U1"]
    draw_bpp_map(ax, data["bpp_u1"], meta, 47.3, 78.7, 21.4)
    add_arrow(ax, (70.1, 89.4), (75.0, 89.4), color=PAIRING, width=0.8)
    add_text(ax, 77.0, 81.0, "BPP overlap", size=8.8,
             weight="bold", color=PAIRING)
    add_text(ax, 77.0, 86.2, "vs sequence", size=8.8,
             weight="bold", color=PAIRING)
    add_text(ax, 77.0, 91.4, "identity", size=8.8,
             weight="bold", color=PAIRING)

    # The calibration branch starts at the blue selected-set marker and routes
    # around the gray non-selected row; null panels are never drawn as units.
    ax.plot([18.0, 14.8, 14.8, 21.8], [88.8, 88.8, 105.5, 105.5],
            color=SEQUENCE, lw=0.65, ls="--", zorder=2)
    add_arrow(ax, (21.8, 105.5), (23.0, 105.5), color=SEQUENCE,
              width=0.65, linestyle="--")
    draw_sheet_stack(ax, 22.8, 103.2)
    add_text(ax, 34.0, 101.6, "64 dinucleotide", size=8.0,
             weight="bold", color=PAIRING)
    add_text(ax, 34.0, 105.5, "whole-panel shuffles", size=8.0,
             color=MUTED)


def culture_marker(ax, x: float, y: float, culture: int, *,
                   facecolor: str, edgecolor: str) -> None:
    if culture == 1:
        ax.add_patch(Circle((x, y), 1.15, facecolor=facecolor,
                            edgecolor=edgecolor, linewidth=0.7, zorder=5))
    elif culture == 2:
        ax.add_patch(Rectangle((x - 1.05, y - 1.05), 2.1, 2.1,
                               facecolor=facecolor, edgecolor=edgecolor,
                               linewidth=0.7, zorder=5))
    else:
        ax.add_patch(Polygon([(x, y - 1.25), (x - 1.2, y + 1.05),
                              (x + 1.2, y + 1.05)], closed=True,
                             facecolor=facecolor, edgecolor=edgecolor,
                             linewidth=0.7, zorder=5))


def draw_sa1_windows(ax, windows: dict[str, dict[str, str]]) -> None:
    whole = windows["array_whole"]
    start = int(whole["start0"])
    end = int(whole["end0"])
    x0, x1 = 118.0, 166.5

    def wx(coordinate: int) -> float:
        return x0 + (coordinate - start) / (end - start) * (x1 - x0)

    ax.add_patch(Rectangle((x0, 72.0), x1 - x0, 5.5, facecolor="white",
                           edgecolor="none", zorder=1))
    ordered = [
        ("array_distal_residual", "", "#C5CDD2"),
        ("array_terminal_residual", "", "#C5CDD2"),
        ("U1_spacer", "U1", UNIT_COLORS["U1"]),
        ("U2_spacer", "U2", UNIT_COLORS["U2"]),
        ("U3_spacer", "U3", UNIT_COLORS["U3"]),
        ("U4_spacer", "U4", UNIT_COLORS["U4"]),
    ]
    for window_id, label, color in ordered:
        row = windows[window_id]
        left, right = wx(int(row["start0"])), wx(int(row["end0"]))
        ax.add_patch(Rectangle((left, 72.0), right - left, 5.5,
                               facecolor=color, edgecolor="white",
                               linewidth=0.55, zorder=3))
        if label:
            add_text(ax, (left + right) / 2, 74.75, label, size=8.0,
                     color="white", weight="bold", ha="center", va="center")
    for index in range(1, 6):
        row = windows[f"repeat{index:02d}_mask"]
        left, right = wx(int(row["start0"])), wx(int(row["end0"]))
        ax.add_patch(Rectangle((left, 72.0), right - left, 5.5,
                               facecolor="#727D84", edgecolor="white",
                               linewidth=0.35, zorder=4))
    ax.add_patch(Rectangle((x0, 72.0), x1 - x0, 5.5, facecolor="none",
                           edgecolor="#7C898F", linewidth=0.55, zorder=4))
    add_text(ax, x0, 79.0, "Repeat-excluded spacers", size=8.0,
             color=MUTED)

    fragment_left, fragment_right = 132.0, 151.5
    fragment_mid = (fragment_left + fragment_right) / 2
    ax.plot([fragment_left, fragment_right], [65.2, 65.2], color=TIME,
            lw=1.4, solid_capstyle="round", zorder=4)
    ax.plot([fragment_left, fragment_left], [64.0, 66.4], color=TIME, lw=0.7)
    ax.plot([fragment_right, fragment_right], [64.0, 66.4], color=TIME, lw=0.7)
    ax.add_patch(Circle((fragment_mid, 65.2), 1.0, facecolor=TIME,
                        edgecolor="white", linewidth=0.55, zorder=5))
    add_text(ax, fragment_mid, 60.6, "Fragment",
             size=8.0, color=MUTED, ha="center")
    add_arrow(ax, (fragment_mid, 67.0), (fragment_mid, 71.0),
              color=TIME, width=0.65)


def draw_timecourse_panel(ax, data: dict[str, object]) -> None:
    ax.plot([113.2, 113.2], [3.0, 111.8], color=FAINT, lw=0.75)
    panel_heading(ax, 116.0, 168.0, 4.0, "c", "SA1 time course", TIME)
    add_text(ax, 116.0, 11.7, "12 libraries · 3 cultures", size=8.8,
             weight="bold", color=INK)
    add_text(ax, 168.0, 17.0, "sampling design · not time-scaled", size=8.0,
             color=MUTED, ha="right")

    x_positions = {0: 127.0, 5: 138.3, 15: 149.5, 55: 162.2}
    add_text(ax, 118.0, 22.0, "min", size=8.0, color=MUTED)
    for time, x in x_positions.items():
        add_text(ax, x, 22.0, str(time), size=8.0, color=INK,
                 weight="bold", ha="center")
    for culture, y in ((1, 30.5), (2, 37.2), (3, 43.9)):
        add_text(ax, 118.0, y, f"C{culture}", size=8.0,
                 color=MUTED, ha="left", va="center")
        ax.plot([x_positions[0], x_positions[55]], [y, y],
                color="#B8C1C6", lw=0.65, zorder=1)
        for time, x in x_positions.items():
            if time == 0:
                face, edge = "white", "#9BA6AC"
            else:
                face, edge = "#D47C36", TIME
            culture_marker(ax, x, y, culture, facecolor=face, edgecolor=edge)

    add_text(ax, 116.0, 49.0, "0 min: low-count background", size=8.0,
             color=MUTED)
    ax.plot([x_positions[5], x_positions[55]], [52.0, 52.0],
            color=TIME, lw=0.7)
    ax.plot([x_positions[5], x_positions[5]], [51.1, 52.9], color=TIME, lw=0.7)
    ax.plot([x_positions[55], x_positions[55]], [51.1, 52.9], color=TIME, lw=0.7)
    add_arrow(ax, ((x_positions[5] + x_positions[55]) / 2, 52.9),
              ((x_positions[5] + x_positions[55]) / 2, 57.2),
              color=TIME, width=0.75)
    add_text(ax, 116.0, 56.3, "Midpoint density", size=8.8,
             weight="bold", color=TIME)

    draw_sa1_windows(ax, data["windows"])
    add_arrow(ax, (142.2, 83.2), (142.2, 90.0), color=TIME, width=0.8)
    add_text(ax, 116.0, 89.8, "Unit shares", size=9.5,
             weight="bold", color=TIME)

    for index, unit in enumerate(("U1", "U2", "U3", "U4")):
        x = 119.5 + index * 11.2
        ax.add_patch(Rectangle((x, 96.8), 8.0, 5.9,
                               facecolor="white", edgecolor=UNIT_COLORS[unit],
                               linewidth=1.15, zorder=3))
        add_text(ax, x + 4.0, 99.75, unit, size=8.1,
                 weight="bold", color=UNIT_COLORS[unit],
                 ha="center", va="center")
    add_text(ax, 116.0, 105.2, "Midpoints / nt", size=8.0,
             color=INK)
    add_text(ax, 116.0, 109.3, "→ Normalize across U1–U4", size=8.0,
             color=INK)


def validate_canvas(fig, ax) -> None:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    canvas = fig.bbox
    for item in ax.texts:
        box = item.get_window_extent(renderer)
        require(canvas.x0 - 0.5 <= box.x0 and box.x1 <= canvas.x1 + 0.5 and
                canvas.y0 - 0.5 <= box.y0 and box.y1 <= canvas.y1 + 0.5,
                f"Text exceeds canvas: {item.get_text()!r}")


def resolve_output_dir(root: Path, requested: Path | None) -> Path:
    if requested is not None:
        return requested.resolve()
    return root / "results" / "reproduced" / "figures"


def render(root: Path, destination: Path) -> dict[str, object]:
    configure_matplotlib()
    data = load_saved_objects(root)
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, WIDTH_MM)
    ax.set_ylim(HEIGHT_MM, 0)
    ax.set_axis_off()

    draw_shared_locus(ax)
    # One shared input rail: the same 37-unit cohort enters both left analyses.
    ax.plot([8.0, 8.0], [25.9, 84.0], color="#A9B4BA", lw=0.65, zorder=1)
    add_arrow(ax, (8.0, 47.8), (14.1, 47.8), color="#8C989F", width=0.65)
    add_arrow(ax, (8.0, 84.0), (14.1, 84.0), color="#8C989F", width=0.65)
    draw_correspondence_panel(ax, data["tracks"])
    draw_pairing_panel(ax, data)
    draw_timecourse_panel(ax, data)

    validate_canvas(fig, ax)
    destination.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for extension in ("svg", "png", "pdf"):
        path = destination / f"Figure_1.{extension}"
        kwargs = {"facecolor": "white", "bbox_inches": None, "pad_inches": 0}
        if extension == "png":
            kwargs["dpi"] = 300
        fig.savefig(path, format=extension, **kwargs)
        outputs[extension] = path
    plt.close(fig)
    return {"data": data, "outputs": outputs}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="Repository root containing data/processed.",
    )
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="Output directory (default: results/reproduced/figures).")
    args = parser.parse_args()
    root = args.project_root.resolve()
    destination = resolve_output_dir(root, args.output_dir)
    receipt = render(root, destination)
    data = receipt["data"]
    print(f"Figure_1: {WIDTH_MM:.0f} x {HEIGHT_MM:.0f} mm")
    print("Shared primary units: 37 (28 x 9 = 252 cross-group combinations)")
    print(f"Selected/other combinations: {data['selected_count']}/{data['other_count']}")
    print(f"Saved representative U1 BPP edges: {len(data['bpp_u1'])}")
    if data["library_receipt"]:
        info = data["library_receipt"]
        print("Saved SA1 libraries: "
              f"{info['rows']}; 0-min spacer midpoint counts "
              f"{info['background_count_min']}-{info['background_count_max']}")
    for extension, path in receipt["outputs"].items():
        print(f"Wrote {extension.upper()}: {path}")


if __name__ == "__main__":
    main()
