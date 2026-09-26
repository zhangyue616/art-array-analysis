#!/usr/bin/env python3
"""Render revised Supplementary Figure S2 from accepted protein tables."""

from __future__ import annotations

import copy

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from common import BLUE, GRID, MISSING, ORANGE, ROOT, TEXT, clean_axis, configure, panel_title, read_tsv, save_figure


DATA_DIR = ROOT / "data" / "processed" / "p6_protein"
ORDER = ["SA1", "MarsHill", "Madawaska", "LY01", "S6", "PALS_2", "UFV_DC4", "AH12", "Machias", "PB50"]


def build_matrices(rows: list[dict[str, str]]) -> dict[str, np.ndarray]:
    keys = {
        "RT N-terminal segment": "rt_ntd_identity_per_alignment_column",
        "RT core segment": "rt_core_identity_per_alignment_column",
        "Partner protein": "partner_identity_per_alignment_column",
    }
    index = {label: position for position, label in enumerate(ORDER)}
    matrices = {title: np.full((len(ORDER), len(ORDER)), np.nan, dtype=float) for title in keys}
    for matrix in matrices.values():
        np.fill_diagonal(matrix, 1.0)
    for row in rows:
        left = index[row["label_a"]]
        right = index[row["label_b"]]
        for title, field in keys.items():
            value = float(row[field])
            matrices[title][left, right] = value
            matrices[title][right, left] = value
    return matrices


def main() -> None:
    configure()
    features = sorted(read_tsv(DATA_DIR / "protein_feature_summary.tsv"), key=lambda row: int(row["locus_order"]))
    pairs = read_tsv(DATA_DIR / "protein_pairwise_identity.tsv")
    if [row["label"] for row in features] != ORDER or len(pairs) != 45:
        raise ValueError("Supplementary Figure S2 requires the accepted ten-locus, 45-pair tables")
    matrices = build_matrices(pairs)
    off_diagonal = np.concatenate([matrix[np.triu_indices_from(matrix, 1)] for matrix in matrices.values()])
    if not np.all(np.isfinite(off_diagonal)):
        raise ValueError("Supplementary Figure S2 heatmaps contain non-finite identity values")
    heatmap_min = float(np.min(off_diagonal))
    heatmap_max = float(np.max(off_diagonal))
    if heatmap_min < 0.0 or heatmap_max > 1.0:
        raise ValueError(f"Supplementary Figure S2 identity range is outside [0, 1]: {heatmap_min} to {heatmap_max}")
    # Identity is a bounded quantity. All three panels therefore use the same
    # explicit semantic range. NaN remains masked and is not rendered as zero.
    normalizer = Normalize(vmin=0.0, vmax=1.0)
    # Match the identity color semantics used by Figure 2a.
    heatmap_cmap = copy.copy(plt.cm.Blues)
    heatmap_cmap.set_bad(color=MISSING)

    fig = plt.figure(figsize=(170 / 25.4, 178 / 25.4))
    grid = GridSpec(
        4,
        3,
        figure=fig,
        height_ratios=(0.10, 1.05, 0.95, 0.10),
        left=0.15,
        right=0.975,
        bottom=0.07,
        top=0.98,
        hspace=0.39,
        wspace=0.22,
    )
    legend_axis = fig.add_subplot(grid[0, :])
    legend_axis.axis("off")
    architecture_grid = GridSpecFromSubplotSpec(
        1, 2, subplot_spec=grid[1, :], width_ratios=(0.74, 0.26), wspace=0.04
    )
    ax_arch = fig.add_subplot(architecture_grid[0, 0])
    ax_metadata = fig.add_subplot(architecture_grid[0, 1], sharey=ax_arch)
    heat_axes = [fig.add_subplot(grid[2, column]) for column in range(3)]
    color_axis = fig.add_subplot(grid[3, :])

    positions = np.arange(len(features))
    ntd = np.array([int(row["mapped_ntd_length_aa"]) for row in features])
    core = np.array([int(row["mapped_core_length_aa"]) for row in features])
    total = ntd + core
    core_colors = [BLUE if row["subgroup"] == "seven_locus_reference" else ORANGE for row in features]
    ax_arch.barh(positions, ntd, color="#D7DDE5", edgecolor="#555555", linewidth=0.35, height=0.70)
    ax_arch.barh(positions, core, left=ntd, color=core_colors, edgecolor="#555555", linewidth=0.35, height=0.70)
    for y_position, row in zip(positions, features):
        motif_position = int(row["c_terminal_yxdd_position_1based"]) - 1
        ax_arch.plot(motif_position, y_position, marker="|", markersize=11, markeredgewidth=1.35, color="#111111", zorder=4)
        ax_metadata.text(
            0.10,
            y_position,
            row["c_terminal_yxdd"],
            ha="left",
            va="center",
            fontsize=8,
        )
        ax_metadata.text(
            0.92,
            y_position,
            row["partner_length_aa"],
            ha="right",
            va="center",
            fontsize=8,
        )
    ax_arch.axhline(6.5, color="#666666", lw=0.70, ls=(0, (3, 2)))
    ax_metadata.axhline(6.5, color="#666666", lw=0.70, ls=(0, (3, 2)))
    ax_arch.set_yticks(positions)
    ax_arch.set_yticklabels(ORDER)
    ax_arch.invert_yaxis()
    ax_arch.set_xlim(0, max(total) + 8)
    ax_arch.set_xlabel("RT residue position")
    panel_title(ax_arch, "a", "Operational RT partitions and partner lengths", pad=8)
    clean_axis(ax_arch, "x")
    ax_metadata.set_xlim(0, 1)
    ax_metadata.axis("off")
    ax_metadata.plot([0.02, 0.02], [-0.5, len(features) - 0.5], color=GRID, lw=0.45, clip_on=False)
    ax_metadata.text(0.10, 1.04, "Motif", transform=ax_metadata.transAxes, ha="left", va="bottom", fontsize=8, color=TEXT)
    ax_metadata.text(0.92, 1.04, "Partner (aa)", transform=ax_metadata.transAxes, ha="right", va="bottom", fontsize=8, color=TEXT)
    architecture_handles = [
        Patch(facecolor="#D7DDE5", edgecolor="#555555", label="Transferred N-terminal"),
        Patch(facecolor=BLUE, edgecolor="#555555", label="Reference core"),
        Patch(facecolor=ORANGE, edgecolor="#555555", label="PB50-related core"),
        Line2D([0], [0], marker="|", color="#111111", linestyle="None", markersize=10, markeredgewidth=1.35, label="YxDD motif"),
    ]
    legend_axis.legend(handles=architecture_handles, loc="center", ncol=4, frameon=False, columnspacing=1.1, handletextpad=0.45)

    image = None
    for letter, axis, (title, matrix) in zip("bcd", heat_axes, matrices.items()):
        image = axis.imshow(np.ma.masked_invalid(matrix), cmap=heatmap_cmap, norm=normalizer, interpolation="none", aspect="equal")
        axis.axhline(6.5, color="white", lw=0.75, alpha=0.90)
        axis.axvline(6.5, color="white", lw=0.75, alpha=0.90)
        axis.set_xticks(list(range(len(ORDER))))
        axis.set_xticklabels(ORDER, rotation=90)
        axis.set_yticks(list(range(len(ORDER))))
        axis.set_yticklabels(ORDER if letter == "b" else [""] * len(ORDER))
        panel_title(axis, letter, title, pad=7, title_size=9.0, letter_offset_pt=-12)
        axis.tick_params(length=0, pad=1.5)
    if image is None:
        raise RuntimeError("Protein heatmaps were not rendered")
    colorbar = fig.colorbar(image, cax=color_axis, orientation="horizontal")
    colorbar.set_ticks(np.linspace(0.0, 1.0, 6))
    colorbar.set_label("Global pairwise identity per alignment column")
    colorbar.ax.tick_params(labelsize=8)

    for path in save_figure(fig, "Figure_S2"):
        print(path)


if __name__ == "__main__":
    main()
