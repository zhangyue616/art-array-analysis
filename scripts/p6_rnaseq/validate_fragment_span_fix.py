#!/usr/bin/env python3
"""Targeted checks for the production SAM outer-template-span helper."""

from __future__ import annotations

import json
from pathlib import Path

from run_library import outer_template_span


ROOT = Path(__file__).resolve().parents[2]
STATE = ROOT / "data" / "source_metadata" / "p6_rnaseq" / "orchestrator_state.json"
RECEIPT = ROOT / "data" / "source_metadata" / "p6_rnaseq" / "fragment_span_fix_receipt.json"


def main() -> None:
    with STATE.open(encoding="utf-8") as handle:
        state = json.load(handle)
    if state.get("stage") != "downloading" or state.get("completed_runs"):
        raise RuntimeError("Targeted fix validation expected downloading state with no completed run")

    cases = [
        {
            "name": "positive_TLEN_standard",
            "position0": 100,
            "mate_position1": 221,
            "template_length": 171,
            "expected": [100, 271],
        },
        {
            "name": "negative_TLEN_standard",
            "position0": 220,
            "mate_position1": 101,
            "template_length": -171,
            "expected": [100, 271],
        },
        {
            "name": "read2_contained_in_mate_reported_example",
            "position0": 150,
            "mate_position1": 101,
            "template_length": -151,
            "expected": [100, 251],
        },
        {
            "name": "mate_contained_in_read2",
            "position0": 100,
            "mate_position1": 151,
            "template_length": 151,
            "expected": [100, 251],
        },
    ]
    for case in cases:
        observed = outer_template_span(
            case["position0"],
            case["mate_position1"],
            case["template_length"],
            1000,
        )
        case["observed"] = list(observed)
        case["pass"] = list(observed) == case["expected"]

    invalid_pnext_rejected = False
    try:
        outer_template_span(100, 0, 151, 1000)
    except ValueError:
        invalid_pnext_rejected = True

    receipt = {
        "status": "pass" if all(case["pass"] for case in cases) and invalid_pnext_rejected else "fail",
        "production_function": "scripts/p6_rnaseq/run_library.py::outer_template_span",
        "formula": "start0=min(POS-1,PNEXT-1); end0=start0+abs(TLEN)",
        "pre_run_state": {
            "stage": state.get("stage"),
            "completed_runs": state.get("completed_runs"),
            "current_run": state.get("current_run"),
            "orchestrator_pid": state.get("orchestrator_pid"),
        },
        "cases": cases,
        "invalid_pnext_rejected": invalid_pnext_rejected,
        "scope": "Coordinate helper only; no FASTQ, alignment, fold, protein, or complete-library rerun.",
        "references": [
            "https://samtools.github.io/hts-specs/SAMv1.pdf (p. 8, POS/PNEXT/TLEN)",
            "https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml (overlapping/containing mates)",
        ],
    }
    with RECEIPT.open("w", encoding="utf-8") as handle:
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if receipt["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
