#!/usr/bin/env python3
"""Render the current Figure 1-9 and Figure S1-S2 numbering.

The accepted data-figure renderers retain their original Figure 1-8 names for
provenance. This wrapper runs them in an isolated staging directory and maps
those outputs to current Figure 2-9, so no renderer writes over another figure.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


FORMATS = ("pdf", "svg", "png")


def parse_args() -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=default_root)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--figures",
        choices=("all", "key"),
        default="all",
        help="Render all 11 figures, or the current key pair (Figures 6 and 8).",
    )
    return parser.parse_args()


def copy_mapped(staging: Path, output: Path, source: str, target: str) -> list[Path]:
    written: list[Path] = []
    for extension in FORMATS:
        source_path = staging / f"Figure_{source}.{extension}"
        if not source_path.is_file():
            raise FileNotFoundError(f"Legacy renderer did not produce {source_path.name}")
        target_path = output / f"Figure_{target}.{extension}"
        shutil.copyfile(source_path, target_path)
        written.append(target_path)
    return written


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    old_renderer = root / "scripts/p6_revision_20260925_r6/figures/render_all.py"
    workflow_renderer = root / "scripts/p6_revision_20260926_r8/render_workflow.py"
    written: list[Path] = []

    with tempfile.TemporaryDirectory(prefix=".legacy_figure_names_", dir=output) as temporary:
        staging = Path(temporary)
        subprocess.run(
            [
                sys.executable,
                str(old_renderer),
                "--project-root",
                str(root),
                "--output-dir",
                str(staging),
                "--figures",
                args.figures,
            ],
            cwd=root,
            check=True,
        )
        mapping = (
            (("5", "6"), ("7", "8"))
            if args.figures == "key"
            else tuple((str(old), str(old + 1)) for old in range(1, 9))
            + (("S1", "S1"), ("S2", "S2"))
        )
        for source, target in mapping:
            written.extend(copy_mapped(staging, output, source, target))

    if args.figures == "all":
        subprocess.run(
            [
                sys.executable,
                str(workflow_renderer),
                "--project-root",
                str(root),
                "--output-dir",
                str(output),
            ],
            cwd=root,
            check=True,
        )
        written.extend(output / f"Figure_1.{extension}" for extension in FORMATS)

    for path in sorted(written):
        if not path.is_file():
            raise FileNotFoundError(path)
        print(path)


if __name__ == "__main__":
    main()
