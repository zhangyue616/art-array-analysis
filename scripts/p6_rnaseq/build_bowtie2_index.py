#!/usr/bin/env python3
"""Build the combined-reference Bowtie2 index required by the full raw route."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REFERENCE = ROOT / "data/processed/p6_rnaseq/MW218148.1_NZ_CP059679.1.fna"
DEFAULT_PREFIX = ROOT / "data/processed/p6_rnaseq/MW218148.1_NZ_CP059679.1"
DEFAULT_BUNDLED = ROOT / "tools/p6_rnaseq/bowtie2-2.5.5-mingw-x86_64/bowtie2-build-s.exe"


def resolve_executable(explicit: str | None) -> Path:
    candidates = [
        explicit,
        os.environ.get("ART_ARRAY_BOWTIE2_BUILD"),
        shutil.which("bowtie2-build"),
        shutil.which("bowtie2-build-s"),
        str(DEFAULT_BUNDLED),
    ]
    for candidate in candidates:
        if candidate and Path(candidate).expanduser().is_file():
            return Path(candidate).expanduser().resolve()
    raise SystemExit(
        "Cannot find bowtie2-build; use --bowtie2-build, set ART_ARRAY_BOWTIE2_BUILD, or add it to PATH"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bowtie2-build")
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--index-prefix", type=Path, default=DEFAULT_PREFIX)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    if args.threads < 1:
        raise SystemExit("--threads must be positive")

    executable = resolve_executable(args.bowtie2_build)
    reference = args.reference if args.reference.is_absolute() else ROOT / args.reference
    prefix = args.index_prefix if args.index_prefix.is_absolute() else ROOT / args.index_prefix
    if not reference.is_file():
        raise SystemExit(f"Combined reference absent: {reference}")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    command = [str(executable), "--threads", str(args.threads), str(reference), str(prefix)]
    subprocess.run(command, cwd=ROOT, check=True)

    suffix_sets = (
        (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
        (".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"),
    )
    outputs = next(
        ([Path(str(prefix) + suffix) for suffix in suffixes] for suffixes in suffix_sets
         if all(Path(str(prefix) + suffix).is_file() for suffix in suffixes)),
        None,
    )
    if outputs is None:
        raise RuntimeError("bowtie2-build exited successfully but the six expected index components are absent")
    print(json.dumps({
        "status": "complete",
        "reference": reference.relative_to(ROOT).as_posix(),
        "index_prefix": prefix.relative_to(ROOT).as_posix(),
        "components": [path.relative_to(ROOT).as_posix() for path in outputs],
    }, indent=2))


if __name__ == "__main__":
    main()
