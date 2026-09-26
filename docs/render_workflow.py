#!/usr/bin/env python3
"""Draw the repository's analysis overview as editable SVG and a 300-dpi PNG."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

WIDTH, HEIGHT = 1040, 750
WIDTH_MM = 170
POINTS_PER_UNIT = WIDTH_MM / 25.4 * 72 / WIDTH
INK = "#1F2C36"
MUTED = "#52616C"


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none",
                         "savefig.facecolor": "white", "figure.facecolor": "white"})
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, WIDTH_MM / 25.4 * HEIGHT / WIDTH))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, WIDTH)
    ax.set_ylim(HEIGHT, 0)
    ax.set_axis_off()

    def text(x, y, label, size=19, bold=False, color=INK):
        return ax.text(x, y, label, fontsize=size * POINTS_PER_UNIT,
                       weight="bold" if bold else "normal", color=color,
                       va="top", ha="left", linespacing=1.32)

    def arrow(x1, x2, y):
        ax.add_patch(FancyArrowPatch((x1, y), (x2, y), arrowstyle="-|>",
                                    mutation_scale=8, linewidth=0.8, color=MUTED))

    text(32, 24, "Public data, three ART-array analyses", 31, True)
    text(32, 73, "Traceable inputs, fixed parameters and reproducible tables and figures.", 19, color=MUTED)
    for x, label in [(52, "PUBLIC INPUTS"), (360, "ANALYSES"), (790, "OUTPUTS")]:
        text(x, 114, label, 18, True, MUTED)

    lanes = [
        {"y": 147, "color": "#176CA4", "fill": "#F1F6FA",
         "name": "01  SEQUENCE",
         "input": ["10 public genome loci", "Oriented ART intervals", "Repeat / unit annotations"],
         "title": "Correspondence + calibration",
         "steps": ["Sequence + boundary rule", "64 composition-matched panels", "Protein identity / locus context"],
         "figures": "Fig. 1–4, S2", "outputs": ["Unit mappings", "Calibrated counts"]},
        {"y": 324, "color": "#20786A", "fill": "#F1F7F5",
         "name": "02  PREDICTED STRUCTURE",
         "input": ["Unit sequences", "One- / two-repeat contexts", "Aligned sequence pairs"],
         "title": "RNA folding ensembles",
         "steps": ["252 cross-group combinations", "BPP overlap vs sequence identity", "Dinucleotide-preserving null"],
         "figures": "Fig. 5–6, S1", "outputs": ["Pairing / identity", "diagnostics"]},
        {"y": 501, "color": "#A95B20", "fill": "#FAF5F0",
         "name": "03  RNA-SEQ",
         "input": ["PRJNA836150", "12 RNA-seq libraries", "3 cultures × 4 time points"],
         "title": "Alignments → fragment counts",
         "steps": ["Quality filtering + mapping", "Culture-paired density shares", "Whole-array residual partition"],
         "figures": "Fig. 7–8", "outputs": ["Unit composition", "Residual coverage"]},
    ]
    for lane in lanes:
        y = lane["y"]
        ax.add_patch(Rectangle((32, y), 976, 158, facecolor=lane["fill"], edgecolor="none"))
        ax.add_patch(Rectangle((32, y), 4, 158, facecolor=lane["color"], edgecolor="none"))
        text(52, y + 14, lane["name"], 18, True, lane["color"])
        for i, label in enumerate(lane["input"]):
            text(52, y + 52 + 29 * i, label, 21 if i == 0 else 18)
        text(360, y + 19, lane["title"], 20, True)
        for i, label in enumerate(lane["steps"]):
            text(360, y + 54 + 29 * i, label, 18)
        arrow(311, 342, y + 88)
        arrow(739, 772, y + 88)
        text(790, y + 41, lane["figures"], 24, True, lane["color"])
        for i, label in enumerate(lane["outputs"]):
            text(790, y + 80 + 29 * i, label, 19)

    text(32, 684, "Offline reproduction: saved tables → statistics + figures", 20, True)
    text(32, 720, "Full route: source accessions → sequences / reads → recorded analysis parameters", 18, color=MUTED)

    # Guard the shared geometry used by both exports; SVG text stays editable.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    canvas = fig.bbox
    for item in ax.texts:
        box = item.get_window_extent(renderer)
        if not (canvas.x0 <= box.x0 and box.x1 <= canvas.x1 and canvas.y0 <= box.y0 and box.y1 <= canvas.y1):
            raise ValueError("Text exceeds the overview canvas: " + item.get_text())
    destination = Path(__file__).resolve().parent
    fig.savefig(destination / "analysis_workflow.svg", format="svg")
    fig.savefig(destination / "analysis_workflow.png", dpi=300)
    plt.close(fig)
    print("Wrote docs/analysis_workflow.svg and docs/analysis_workflow.png")


if __name__ == "__main__":
    main()
