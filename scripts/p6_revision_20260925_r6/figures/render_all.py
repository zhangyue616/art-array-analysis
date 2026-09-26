#!/usr/bin/env python3
"""Render the complete or compact article figure set from saved tables."""

from __future__ import annotations

import argparse
import os
from pathlib import Path


def parse_args() -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=default_root)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--figures",
        choices=("all", "key"),
        default="all",
        help="Render all ten figures, or the compact key pair (Figures 5 and 7).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    output_dir = (args.output_dir or root / "results" / "reproduced" / "figures").resolve()
    os.environ["ART_ARRAY_PROJECT_ROOT"] = str(root)
    os.environ["ART_ARRAY_FIGURE_OUTPUT"] = str(output_dir)

    # Import after configuring paths because the renderers bind their input
    # directories from common.py at import time.
    import common
    import render_figure_s2
    import render_figures_1_4
    import render_figures_5_6_s1
    import render_figures_7_8

    common.configure()
    if args.figures == "key":
        outputs = []
        outputs.extend(render_figures_5_6_s1.render_figure_5())
        outputs.extend(render_figures_7_8.render_figure_7())
        for path in outputs:
            print(path)
        return

    render_figures_1_4.main()
    render_figures_5_6_s1.main()
    render_figures_7_8.main()
    render_figure_s2.main()


if __name__ == "__main__":
    main()
