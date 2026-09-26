#!/usr/bin/env python3
"""Resumable downloader and complete-library runner for all 12 P6 RNA-seq runs."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "p6_rnaseq"
PROCESSED = ROOT / "data" / "processed" / "p6_rnaseq"
META = ROOT / "data" / "source_metadata" / "p6_rnaseq"
TOOLS = ROOT / "tools" / "p6_rnaseq"
ARIA2 = TOOLS / "aria2-1.37.0-win-64bit-build1" / "aria2c.exe"
PYTHON = Path(sys.executable)
RUNNER = ROOT / "scripts" / "p6_rnaseq" / "run_library.py"
ENA_META = META / "ena_PRJNA836150_read_run.json"
STATE_PATH = META / "orchestrator_state.json"
STOP_PATH = META / "STOP_ORCHESTRATOR"
MASTER_LOG = META / "orchestrator.log"
FIRST_LIBRARY = "SRR19152327"

RUN_ORDER = [
    "SRR19152327",  # 15 min, culture 1: benchmark first
    "SRR19152326",
    "SRR19152325",
    "SRR19152330",
    "SRR19152329",
    "SRR19152328",
    "SRR19152324",
    "SRR19152333",
    "SRR19152332",
    "SRR19152335",
    "SRR19152334",
    "SRR19152331",
]


def resolve_aria2(explicit: str | None) -> Path:
    candidates = [explicit, os.environ.get("ART_ARRAY_ARIA2"), shutil.which("aria2c"), str(ARIA2)]
    for candidate in candidates:
        if candidate and Path(candidate).expanduser().is_file():
            return Path(candidate).expanduser().resolve()
    raise SystemExit("Cannot find aria2c; use --aria2, set ART_ARRAY_ARIA2, or add aria2c to PATH")


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    # A reader can briefly hold the destination open on Windows.  Preserve the
    # atomic replace, but tolerate that transient sharing violation.
    for attempt in range(10):
        try:
            temporary.replace(path)
            return
        except PermissionError:
            if attempt == 9:
                raise
            time.sleep(0.1)


def append_log(message: str) -> None:
    with MASTER_LOG.open("a", encoding="utf-8") as handle:
        handle.write(f"{now()} {message}\n")


def load_metadata() -> dict[str, dict[str, str]]:
    with ENA_META.open(encoding="utf-8") as handle:
        return {row["run_accession"]: row for row in json.load(handle)}


def completed_runs() -> list[str]:
    complete = []
    for run in RUN_ORDER:
        status_path = PROCESSED / "runs" / run / "run_status.json"
        if status_path.exists():
            try:
                with status_path.open(encoding="utf-8") as handle:
                    if json.load(handle).get("status") == "complete":
                        complete.append(run)
            except (OSError, json.JSONDecodeError):
                pass
    return complete


def state(stage: str, **extra: object) -> None:
    value: dict[str, object] = {
        "orchestrator_pid": os.getpid(),
        "python": str(PYTHON),
        "stage": stage,
        "updated_utc": now(),
        "run_order": RUN_ORDER,
        "completed_runs": completed_runs(),
        "stop_file": str(STOP_PATH.relative_to(ROOT)),
    }
    value.update(extra)
    write_json(STATE_PATH, value)


def check_stop() -> None:
    if STOP_PATH.exists():
        state("stopped_by_stop_file")
        append_log("stop file observed; exiting between atomic operations")
        raise SystemExit(0)


def download_pair(run: str, row: dict[str, str], aria2: Path) -> None:
    urls = ["https://" + item for item in row["fastq_ftp"].split(";")]
    names = [Path(item).name for item in row["fastq_ftp"].split(";")]
    expected_bytes = [int(item) for item in row["fastq_bytes"].split(";")]
    if len(urls) != 2:
        raise RuntimeError(f"Expected paired FASTQ URLs for {run}")

    processes: list[subprocess.Popen[bytes]] = []
    log_handles = []
    commands = []
    try:
        for mate, (url, name) in enumerate(zip(urls, names), 1):
            log_path = META / f"{run}_R{mate}_download.log"
            log_handle = log_path.open("ab")
            log_handles.append(log_handle)
            command = [
                str(aria2),
                "--continue=true",
                "--allow-overwrite=true",
                "--auto-file-renaming=false",
                "--max-connection-per-server=8",
                "--split=8",
                "--min-split-size=4M",
                "--max-tries=20",
                "--retry-wait=5",
                "--timeout=60",
                "--lowest-speed-limit=1K",
                "--summary-interval=60",
                "--console-log-level=notice",
                "--file-allocation=none",
                f"--dir={RAW}",
                f"--out={name}",
                url,
            ]
            commands.append(command)
            creation_flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
            processes.append(
                subprocess.Popen(
                    command,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    creationflags=creation_flags,
                )
            )

        append_log(f"started paired download for {run}: pids={[process.pid for process in processes]}")
        started = time.monotonic()
        while any(process.poll() is None for process in processes):
            state(
                "downloading",
                current_run=run,
                download_pids=[process.pid for process in processes],
                download_returncodes=[process.poll() for process in processes],
                expected_fastq_bytes=expected_bytes,
                fastq_paths=[str((RAW / name).relative_to(ROOT)) for name in names],
                aria2_control_files=[str((RAW / (name + ".aria2")).relative_to(ROOT)) for name in names],
                elapsed_seconds=time.monotonic() - started,
                commands=commands,
            )
            time.sleep(60)
        returncodes = [process.wait() for process in processes]
        if any(code != 0 for code in returncodes):
            raise RuntimeError(f"aria2 download failure for {run}: {returncodes}")
        for name, expected in zip(names, expected_bytes):
            observed = (RAW / name).stat().st_size
            if observed != expected:
                raise RuntimeError(f"FASTQ byte mismatch after download for {name}: {observed} != {expected}")
        append_log(f"completed paired download for {run}")
    finally:
        for process in processes:
            if process.poll() is None:
                process.terminate()
        for handle in log_handles:
            handle.close()


def run_library(
    run: str,
    threads: int,
    fastp: str | None,
    bowtie2: str | None,
    bowtie2_index: Path | None,
) -> None:
    out_dir = PROCESSED / "runs" / run
    out_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = out_dir / "runner.stdout.log"
    stderr_path = out_dir / "runner.stderr.log"
    command = [str(PYTHON), str(RUNNER), run, "--threads", str(threads)]
    if fastp:
        command.extend(["--fastp", fastp])
    if bowtie2:
        command.extend(["--bowtie2", bowtie2])
    if bowtie2_index:
        command.extend(["--bowtie2-index", str(bowtie2_index)])
    creation_flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    with stdout_path.open("ab") as stdout_handle, stderr_path.open("ab") as stderr_handle:
        process = subprocess.Popen(
            command,
            stdout=stdout_handle,
            stderr=stderr_handle,
            creationflags=creation_flags,
        )
        append_log(f"started complete-library analysis for {run}: pid={process.pid}")
        started = time.monotonic()
        while process.poll() is None:
            state(
                "analyzing",
                current_run=run,
                analysis_pid=process.pid,
                elapsed_seconds=time.monotonic() - started,
                command=command,
            )
            time.sleep(60)
        returncode = process.wait()
    if returncode != 0:
        raise RuntimeError(f"complete-library analysis failed for {run}: return code {returncode}")
    append_log(f"completed analysis for {run}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--aria2", help="aria2c executable (default: ART_ARRAY_ARIA2, PATH, bundled fallback)")
    parser.add_argument("--fastp", help="forwarded to run_library.py")
    parser.add_argument("--bowtie2", help="forwarded to run_library.py")
    parser.add_argument("--bowtie2-index", type=Path, help="forwarded to run_library.py")
    parser.add_argument(
        "--stop-after-first-library",
        action="store_true",
        help="Stop cleanly after the first complete library for an optional user-controlled check.",
    )
    args = parser.parse_args()
    if args.threads < 1 or args.threads > 12:
        raise SystemExit("--threads must be between 1 and 12")
    aria2 = resolve_aria2(args.aria2)

    META.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    metadata = load_metadata()
    if set(metadata) != set(RUN_ORDER):
        raise SystemExit("Frozen ENA metadata does not contain exactly the 12 planned runs")
    state("starting")
    append_log(f"orchestrator started pid={os.getpid()}")
    try:
        for run in RUN_ORDER:
            check_stop()
            status_path = PROCESSED / "runs" / run / "run_status.json"
            run_complete = False
            if status_path.exists():
                with status_path.open(encoding="utf-8") as handle:
                    prior = json.load(handle)
                if prior.get("status") == "complete":
                    append_log(f"skipping already-complete run {run}")
                    run_complete = True
            if not run_complete:
                download_pair(run, metadata[run], aria2)
                check_stop()
                run_library(run, args.threads, args.fastp, args.bowtie2, args.bowtie2_index)
            if run == FIRST_LIBRARY and args.stop_after_first_library:
                state("stopped_after_first_library", current_run=FIRST_LIBRARY)
                append_log("stopped after the first complete library by user request")
                return
        state("all_runs_complete")
        append_log("all 12 complete libraries finished")
    except BaseException as exc:
        state("failed", error_type=type(exc).__name__, error=str(exc))
        append_log(f"orchestrator failed: {type(exc).__name__}: {exc}")
        raise


if __name__ == "__main__":
    main()
