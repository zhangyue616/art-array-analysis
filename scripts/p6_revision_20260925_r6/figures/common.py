#!/usr/bin/env python3
"""Shared rendering helpers for the P6 article figure set."""

from __future__ import annotations

import csv
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation


ROOT = Path(
    os.environ.get("ART_ARRAY_PROJECT_ROOT", Path(__file__).resolve().parents[3])
).resolve()
OUTPUT_DIR = Path(
    os.environ.get("ART_ARRAY_FIGURE_OUTPUT", ROOT / "results" / "reproduced" / "figures")
).resolve()

TEXT = "#252525"
GRID = "#D9D9D9"
BLUE = "#3775BA"
ORANGE = "#D97932"
TEAL = "#2C8C7C"
PURPLE = "#9A4D8E"
RED = "#B5473C"
GRAY = "#7F8C99"
LIGHT_GRAY = "#E8E8E8"
MISSING = "#ECEFF1"

UNIT_COLORS = {
    "U1": "#0072B2",
    "U2": "#009E73",
    "U3": "#CC79A7",
    "U4": "#D55E00",
}


def configure() -> None:
    """Set final-size typography without leaking settings outside each process."""
    matplotlib.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.titlesize": 9.2,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "legend.fontsize": 8.0,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "lines.linewidth": 1.1,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.axisbelow": True,
            "legend.frameon": False,
        }
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save_figure(fig: plt.Figure, label: str) -> list[Path]:
    """Write vector PDF/SVG and a 300-dpi PNG from the same fixed-size canvas."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for suffix in ("pdf", "svg", "png"):
        path = OUTPUT_DIR / f"{label}.{suffix}"
        fig.savefig(path, dpi=300 if suffix == "png" else None, facecolor="white")
        outputs.append(path)
    plt.close(fig)
    return outputs


def panel_label(axis: plt.Axes, text: str, x: float = -0.12, y: float = 1.04) -> None:
    axis.text(
        x,
        y,
        text,
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        fontweight="bold",
        color=TEXT,
    )


def panel_title(axis: plt.Axes, letter: str, title: str, pad: float = 8,
                title_size: float = 9.2, letter_offset_pt: float = -13) -> None:
    """Align a bold panel letter and regular title at physical point offsets."""
    title_transform = axis.transAxes + ScaledTranslation(0, pad / 72, axis.figure.dpi_scale_trans)
    letter_transform = axis.transAxes + ScaledTranslation(
        letter_offset_pt / 72, pad / 72, axis.figure.dpi_scale_trans
    )
    axis.text(0, 1, title, transform=title_transform, ha="left", va="baseline",
              fontsize=title_size, fontweight="normal", color=TEXT, clip_on=False)
    axis.text(0, 1, letter, transform=letter_transform, ha="left", va="baseline",
              fontsize=10, fontweight="bold", color=TEXT, clip_on=False)


def clean_axis(axis: plt.Axes, grid_axis: str | None = None) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    if grid_axis:
        axis.grid(axis=grid_axis, color=GRID, linewidth=0.40, alpha=0.55, zorder=0)
