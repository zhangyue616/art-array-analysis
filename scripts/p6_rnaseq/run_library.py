#!/usr/bin/env python3
"""Run one complete PRJNA836150 library and aggregate SA1 fragments in-stream."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from array import array
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "p6_rnaseq"
PROCESSED = ROOT / "data" / "processed" / "p6_rnaseq"
META = ROOT / "data" / "source_metadata" / "p6_rnaseq"
TOOLS = ROOT / "tools" / "p6_rnaseq"

FASTP = TOOLS / "fastp-1.3.3-windows-ucrt64" / "fastp.exe"
BOWTIE2 = TOOLS / "bowtie2-2.5.5-mingw-x86_64" / "bowtie2-align-s.exe"
INDEX = PROCESSED / "MW218148.1_NZ_CP059679.1"
FEATURES = PROCESSED / "sa1_features.tsv"
WINDOWS = PROCESSED / "sa1_windows.tsv"
ENA_META = META / "ena_PRJNA836150_read_run.json"
PHAGE_NAME = b"MW218148.1"
PHAGE_LENGTH = 260727
COVERAGE_START0 = 9000
COVERAGE_END0 = 13200
CIGAR_TOKEN = re.compile(rb"(\d+)([MIDNSHP=X])")


def resolve_executable(
    explicit: str | None,
    environment_name: str,
    path_names: tuple[str, ...],
    bundled_fallback: Path,
) -> Path:
    """Resolve a command from CLI, environment, PATH, then the historical bundle."""
    candidates = []
    if explicit:
        candidates.append(explicit)
    if os.environ.get(environment_name):
        candidates.append(os.environ[environment_name])
    for name in path_names:
        discovered = shutil.which(name)
        if discovered:
            candidates.append(discovered)
    candidates.append(str(bundled_fallback))
    for candidate in candidates:
        path = Path(candidate).expanduser()
        if path.is_file():
            return path.resolve()
    raise SystemExit(
        f"Cannot find {path_names[0]}; use the corresponding CLI option, "
        f"set {environment_name}, or add it to PATH"
    )


def iso_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def file_digest(path: Path, algorithm: str = "md5") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def reference_length(cigar: bytes) -> int:
    if cigar == b"*":
        return 0
    return sum(int(length) for length, op in CIGAR_TOKEN.findall(cigar) if op in b"MDN=X")


def outer_template_span(
    position0: int,
    mate_position1: int,
    template_length: int,
    reference_size: int,
) -> tuple[int, int]:
    """Return the SAM outer template span as a 0-based half-open interval.

    POS and PNEXT are the leftmost aligned positions of the current read and
    its mate.  TLEN is the signed outer template length.  Using the current
    read's CIGAR right edge is incorrect when one mate is contained within the
    other alignment.
    """

    mate_position0 = mate_position1 - 1
    if position0 < 0 or mate_position0 < 0 or template_length == 0:
        raise ValueError("invalid POS, PNEXT, or zero TLEN")
    fragment_start0 = min(position0, mate_position0)
    fragment_end0 = fragment_start0 + abs(template_length)
    if fragment_end0 <= fragment_start0 or fragment_end0 > reference_size:
        raise ValueError("template span lies outside the reference")
    return fragment_start0, fragment_end0


def parse_metadata(run_accession: str) -> dict[str, str]:
    with ENA_META.open(encoding="utf-8") as handle:
        rows = json.load(handle)
    try:
        return next(row for row in rows if row["run_accession"] == run_accession)
    except StopIteration as exc:
        raise RuntimeError(f"Run absent from frozen ENA metadata: {run_accession}") from exc


def parse_time_culture(sample_alias: str) -> tuple[int, int]:
    match = re.search(r"_(\d+)-(\d+)$", sample_alias)
    if not match:
        raise RuntimeError(f"Cannot parse time/culture from sample alias {sample_alias!r}")
    return int(match.group(1)), int(match.group(2))


def load_windows() -> list[dict[str, object]]:
    with WINDOWS.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [
        {
            **row,
            "start0": int(row["start0"]),
            "end0": int(row["end0"]),
            "length_nt": int(row["length_nt"]),
        }
        for row in rows
    ]


def load_features() -> list[dict[str, object]]:
    with FEATURES.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    parsed = [
        {
            **row,
            "feature_ordinal": int(row["feature_ordinal"]),
            "start0": int(row["start0"]),
            "end0": int(row["end0"]),
            "length_nt": int(row["length_nt"]),
        }
        for row in rows
    ]
    if len(parsed) != 258:
        raise RuntimeError(f"Expected 258 features, found {len(parsed)}")
    return parsed


def build_feature_lookup(features: list[dict[str, object]]) -> tuple[array, dict[int, tuple[int, ...]]]:
    lookup = array("h", [-1]) * PHAGE_LENGTH
    overlaps: dict[int, tuple[int, ...]] = {}
    for index, feature in enumerate(features):
        for position in range(int(feature["start0"]), int(feature["end0"])):
            current = lookup[position]
            if current == -1:
                lookup[position] = index
            elif current >= 0:
                overlaps[position] = (current, index)
                lookup[position] = -2
            else:
                overlaps[position] = overlaps[position] + (index,)
    return lookup, overlaps


def write_status(path: Path, content: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(content, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_accession")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--fastp", help="fastp executable (default: ART_ARRAY_FASTP, PATH, bundled fallback)")
    parser.add_argument(
        "--bowtie2",
        help="Bowtie2 executable (default: ART_ARRAY_BOWTIE2, PATH, bundled fallback)",
    )
    parser.add_argument(
        "--bowtie2-index",
        type=Path,
        default=INDEX,
        help="Bowtie2 index prefix (default: published combined reference under data/processed)",
    )
    args = parser.parse_args()

    run = args.run_accession
    if not re.fullmatch(r"SRR\d+", run):
        raise SystemExit(f"Invalid run accession: {run}")
    if args.threads < 1 or args.threads > 12:
        raise SystemExit("--threads must be between 1 and 12")

    fastp_executable = resolve_executable(args.fastp, "ART_ARRAY_FASTP", ("fastp",), FASTP)
    bowtie2_executable = resolve_executable(
        args.bowtie2,
        "ART_ARRAY_BOWTIE2",
        ("bowtie2", "bowtie2-align-s"),
        BOWTIE2,
    )
    index_prefix = args.bowtie2_index
    if not index_prefix.is_absolute():
        index_prefix = ROOT / index_prefix

    for path in (FEATURES, WINDOWS, ENA_META):
        if not path.exists():
            raise SystemExit(f"Required input absent: {path}")
    for suffix in (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"):
        if not Path(str(index_prefix) + suffix).exists():
            raise SystemExit(f"Bowtie2 index component absent: {index_prefix}{suffix}")

    run_meta = parse_metadata(run)
    time_min, culture = parse_time_culture(run_meta["sample_alias"])
    ftp_names = [Path(item).name for item in run_meta["fastq_ftp"].split(";")]
    expected_md5 = run_meta["fastq_md5"].split(";")
    expected_bytes = [int(item) for item in run_meta["fastq_bytes"].split(";")]
    if len(ftp_names) != 2 or len(expected_md5) != 2 or len(expected_bytes) != 2:
        raise SystemExit(f"Expected two FASTQ files in ENA row for {run}")
    read1, read2 = (RAW / ftp_names[0], RAW / ftp_names[1])
    for index, path in enumerate((read1, read2)):
        if not path.exists():
            raise SystemExit(f"FASTQ absent: {path}")
        if path.stat().st_size != expected_bytes[index]:
            raise SystemExit(
                f"FASTQ byte mismatch for {path.name}: {path.stat().st_size} != {expected_bytes[index]}"
            )

    out_dir = PROCESSED / "runs" / run
    out_dir.mkdir(parents=True, exist_ok=True)
    status_path = out_dir / "run_status.json"
    if status_path.exists() and not args.force:
        with status_path.open(encoding="utf-8") as handle:
            old_status = json.load(handle)
        if old_status.get("status") == "complete":
            print(f"{run} already complete; use --force to rerun")
            return

    started_iso = iso_now()
    started_clock = time.perf_counter()
    status: dict[str, object] = {
        "run_accession": run,
        "sample_alias": run_meta["sample_alias"],
        "time_min": time_min,
        "culture": culture,
        "status": "validating_fastq",
        "started_utc": started_iso,
    }
    write_status(status_path, status)

    observed_md5 = [file_digest(read1), file_digest(read2)]
    if observed_md5 != expected_md5:
        status.update({"status": "failed_fastq_md5", "observed_md5": observed_md5, "expected_md5": expected_md5})
        write_status(status_path, status)
        raise SystemExit(f"FASTQ MD5 mismatch for {run}")
    status.update(
        {
            "fastq_validation": "pass",
            "expected_fastq_bytes": expected_bytes,
            "observed_fastq_bytes": [read1.stat().st_size, read2.stat().st_size],
            "expected_fastq_md5": expected_md5,
            "observed_fastq_md5": observed_md5,
        }
    )

    windows = load_windows()
    features = load_features()
    feature_lookup, overlapping_features = build_feature_lookup(features)
    feature_sense = [0] * len(features)
    feature_antisense = [0] * len(features)
    window_stats = {
        str(window["window_id"]): {
            "sense_midpoint_fragments": 0,
            "antisense_midpoint_fragments": 0,
            "sense_span_bases": 0,
            "antisense_span_bases": 0,
        }
        for window in windows
    }
    plus_difference = array("q", [0]) * (COVERAGE_END0 - COVERAGE_START0 + 1)
    minus_difference = array("q", [0]) * (COVERAGE_END0 - COVERAGE_START0 + 1)
    counters: Counter[str] = Counter()

    fastp_json = out_dir / "fastp.json"
    fastp_html = out_dir / "fastp.html"
    fastp_log = out_dir / "fastp.stderr.log"
    bowtie_log = out_dir / "bowtie2.stderr.log"
    parser_log = out_dir / "parser.log"
    fastp_threads = min(2, args.threads)
    # fastp 1.3.3 for Windows copies argv into its JSON report without escaping
    # backslashes. Forward-slash paths keep that report valid JSON.
    fastp_command = [
        fastp_executable.as_posix(),
        "--in1",
        read1.as_posix(),
        "--in2",
        read2.as_posix(),
        "--stdout",
        "--length_required",
        "30",
        "--thread",
        str(fastp_threads),
        "--json",
        fastp_json.as_posix(),
        "--html",
        fastp_html.as_posix(),
        "--report_title",
        f"P6 {run}",
    ]
    bowtie_command = [
        str(bowtie2_executable),
        "--very-sensitive",
        "-X",
        "1000",
        "--no-unal",
        "--no-head",
        "--no-sq",
        "-p",
        str(args.threads),
        "-x",
        str(index_prefix),
        "--interleaved",
        "-",
    ]
    status.update(
        {
            "status": "running",
            "expected_md5": expected_md5,
            "observed_md5": observed_md5,
            "expected_fastq_bytes": expected_bytes,
            "observed_fastq_bytes": [read1.stat().st_size, read2.stat().st_size],
            "fastp_command": fastp_command,
            "bowtie2_command": bowtie_command,
            "sam_filter_mode": "inline exact RNAME match to MW218148.1",
        }
    )
    write_status(status_path, status)

    with fastp_log.open("wb") as fastp_stderr, bowtie_log.open("wb") as bowtie_stderr:
        fastp_process = subprocess.Popen(
            fastp_command,
            stdout=subprocess.PIPE,
            stderr=fastp_stderr,
            bufsize=1024 * 1024,
        )
        assert fastp_process.stdout is not None
        bowtie_process = subprocess.Popen(
            bowtie_command,
            stdin=fastp_process.stdout,
            stdout=subprocess.PIPE,
            stderr=bowtie_stderr,
            bufsize=1024 * 1024,
        )
        fastp_process.stdout.close()
        assert bowtie_process.stdout is not None
        with parser_log.open("w", encoding="utf-8") as progress_handle:
            for line_number, raw_line in enumerate(bowtie_process.stdout, 1):
                counters["sam_lines"] += 1
                if raw_line.startswith(b"@"):
                    continue
                fields = raw_line.rstrip(b"\r\n").split(b"\t")
                if len(fields) < 11:
                    counters["malformed_sam_lines"] += 1
                    continue
                if fields[2] != PHAGE_NAME:
                    continue
                counters["phage_alignment_records"] += 1
                try:
                    flag = int(fields[1])
                    position0 = int(fields[3]) - 1
                    mapq = int(fields[4])
                    mate_position1 = int(fields[7])
                    template_length = int(fields[8])
                except ValueError:
                    counters["malformed_numeric_sam_lines"] += 1
                    continue
                if not (flag & 0x80):
                    counters["not_read2"] += 1
                    continue
                if flag & (0x100 | 0x800):
                    counters["secondary_or_supplementary_read2"] += 1
                    continue
                counters["primary_read2_phage"] += 1
                if flag & (0x4 | 0x8):
                    counters["filtered_unmapped_mate"] += 1
                    continue
                if not (flag & 0x2):
                    counters["filtered_not_proper_pair"] += 1
                    continue
                if fields[6] not in (b"=", PHAGE_NAME):
                    counters["filtered_mate_other_reference"] += 1
                    continue
                if mapq < 10:
                    counters["filtered_mapq_lt10"] += 1
                    continue
                if template_length == 0:
                    counters["filtered_zero_template"] += 1
                    continue
                if abs(template_length) > 1500:
                    counters["filtered_template_gt1500"] += 1
                    continue
                cigar_reference_length = reference_length(fields[5])
                if cigar_reference_length <= 0:
                    counters["filtered_invalid_cigar"] += 1
                    continue
                if mate_position1 <= 0:
                    counters["filtered_invalid_pnext"] += 1
                    continue
                try:
                    fragment_start0, fragment_end0 = outer_template_span(
                        position0,
                        mate_position1,
                        template_length,
                        PHAGE_LENGTH,
                    )
                except ValueError:
                    counters["filtered_invalid_fragment_span"] += 1
                    continue

                counters["retained_phage_fragments"] += 1
                strand = "-" if flag & 0x10 else "+"
                counters[f"retained_phage_fragments_{'minus' if strand == '-' else 'plus'}"] += 1
                midpoint0 = (fragment_start0 + fragment_end0) // 2

                feature_code = feature_lookup[midpoint0]
                if feature_code >= 0:
                    feature_indices = (feature_code,)
                elif feature_code == -2:
                    feature_indices = overlapping_features[midpoint0]
                    counters["midpoints_in_overlapping_CDS"] += 1
                else:
                    feature_indices = ()
                if not feature_indices:
                    counters["midpoints_outside_CDS"] += 1
                for feature_index in feature_indices:
                    if strand == features[feature_index]["strand"]:
                        feature_sense[feature_index] += 1
                    else:
                        feature_antisense[feature_index] += 1

                if fragment_end0 > COVERAGE_START0 and fragment_start0 < COVERAGE_END0:
                    clipped_start = max(fragment_start0, COVERAGE_START0) - COVERAGE_START0
                    clipped_end = min(fragment_end0, COVERAGE_END0) - COVERAGE_START0
                    difference = plus_difference if strand == "+" else minus_difference
                    difference[clipped_start] += 1
                    difference[clipped_end] -= 1

                    for window in windows:
                        overlap = max(
                            0,
                            min(fragment_end0, int(window["end0"]))
                            - max(fragment_start0, int(window["start0"])),
                        )
                        if overlap:
                            key = "sense_span_bases" if strand == window["strand"] else "antisense_span_bases"
                            window_stats[str(window["window_id"])][key] += overlap

                for window in windows:
                    if int(window["start0"]) <= midpoint0 < int(window["end0"]):
                        key = "sense_midpoint_fragments" if strand == window["strand"] else "antisense_midpoint_fragments"
                        window_stats[str(window["window_id"])][key] += 1

                if counters["retained_phage_fragments"] % 1_000_000 == 0:
                    progress_handle.write(
                        f"{iso_now()} retained_phage_fragments={counters['retained_phage_fragments']} "
                        f"sam_lines={counters['sam_lines']}\n"
                    )
                    progress_handle.flush()

        bowtie_process.stdout.close()
        bowtie_return = bowtie_process.wait()
        fastp_return = fastp_process.wait()

    if fastp_return != 0 or bowtie_return != 0:
        status.update(
            {
                "status": "failed_pipeline",
                "fastp_returncode": fastp_return,
                "bowtie2_returncode": bowtie_return,
                "findstr_returncode": 0,
                "sam_filter_mode": "inline exact RNAME match to MW218148.1",
                "parser_counters": dict(counters),
                "finished_utc": iso_now(),
                "elapsed_seconds": time.perf_counter() - started_clock,
            }
        )
        write_status(status_path, status)
        raise SystemExit(
            f"Pipeline failed for {run}: fastp={fastp_return}, bowtie2={bowtie_return}"
        )

    plus_depth = array("q", [0]) * (COVERAGE_END0 - COVERAGE_START0)
    minus_depth = array("q", [0]) * (COVERAGE_END0 - COVERAGE_START0)
    running_plus = 0
    running_minus = 0
    for index in range(COVERAGE_END0 - COVERAGE_START0):
        running_plus += plus_difference[index]
        running_minus += minus_difference[index]
        plus_depth[index] = running_plus
        minus_depth[index] = running_minus

    coverage_path = out_dir / "locus_coverage.tsv.gz"
    retained_phage = counters["retained_phage_fragments"]
    with gzip.open(coverage_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "run_accession",
                "time_min",
                "culture",
                "position_1based",
                "plus_fragment_depth",
                "minus_fragment_depth",
                "plus_depth_per_million_retained_phage_fragments",
            ]
        )
        for offset, (plus, minus) in enumerate(zip(plus_depth, minus_depth)):
            writer.writerow(
                [
                    run,
                    time_min,
                    culture,
                    COVERAGE_START0 + offset + 1,
                    plus,
                    minus,
                    f"{(plus * 1_000_000 / retained_phage) if retained_phage else 0:.9f}",
                ]
            )

    window_path = out_dir / "window_counts.tsv"
    with window_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "run_accession",
            "sample_alias",
            "time_min",
            "culture",
            "window_id",
            "window_class",
            "start0",
            "end0",
            "length_nt",
            "sense_midpoint_fragments",
            "antisense_midpoint_fragments",
            "sense_midpoint_density_per_nt",
            "sense_span_bases",
            "antisense_span_bases",
            "sense_mean_fragment_depth",
            "antisense_mean_fragment_depth",
            "sense_coverage_breadth_ge1",
            "antisense_coverage_breadth_ge1",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for window in windows:
            name = str(window["window_id"])
            stats = window_stats[name]
            start = int(window["start0"])
            end = int(window["end0"])
            length = int(window["length_nt"])
            local_start = start - COVERAGE_START0
            local_end = end - COVERAGE_START0
            if local_start < 0 or local_end > len(plus_depth):
                raise RuntimeError(f"Coverage window outside stored interval: {name}")
            sense_vector = plus_depth[local_start:local_end] if window["strand"] == "+" else minus_depth[local_start:local_end]
            antisense_vector = minus_depth[local_start:local_end] if window["strand"] == "+" else plus_depth[local_start:local_end]
            writer.writerow(
                {
                    "run_accession": run,
                    "sample_alias": run_meta["sample_alias"],
                    "time_min": time_min,
                    "culture": culture,
                    "window_id": name,
                    "window_class": window["window_class"],
                    "start0": start,
                    "end0": end,
                    "length_nt": length,
                    **stats,
                    "sense_midpoint_density_per_nt": f"{stats['sense_midpoint_fragments'] / length:.12g}",
                    "sense_mean_fragment_depth": f"{stats['sense_span_bases'] / length:.12g}",
                    "antisense_mean_fragment_depth": f"{stats['antisense_span_bases'] / length:.12g}",
                    "sense_coverage_breadth_ge1": f"{sum(value >= 1 for value in sense_vector) / length:.12g}",
                    "antisense_coverage_breadth_ge1": f"{sum(value >= 1 for value in antisense_vector) / length:.12g}",
                }
            )

    array_stats = window_stats["array_whole"]
    array_rpk = array_stats["sense_midpoint_fragments"] / (1199 / 1000)
    feature_rpk = [count / (int(feature["length_nt"]) / 1000) for count, feature in zip(feature_sense, features)]
    tpm_denominator = array_rpk + sum(feature_rpk)
    feature_path = out_dir / "feature_counts.tsv"
    with feature_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "run_accession",
            "sample_alias",
            "time_min",
            "culture",
            "feature_id",
            "feature_type",
            "protein_id",
            "locus_tag",
            "product",
            "start0",
            "end0",
            "length_nt",
            "strand",
            "sense_midpoint_fragments",
            "antisense_midpoint_fragments",
            "sense_rpk",
            "sense_tpm_259_features",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {
                "run_accession": run,
                "sample_alias": run_meta["sample_alias"],
                "time_min": time_min,
                "culture": culture,
                "feature_id": "array_RNA",
                "feature_type": "array_RNA",
                "protein_id": "",
                "locus_tag": "",
                "product": "published whole non-coding region upstream of RT",
                "start0": 9447,
                "end0": 10646,
                "length_nt": 1199,
                "strand": "+",
                "sense_midpoint_fragments": array_stats["sense_midpoint_fragments"],
                "antisense_midpoint_fragments": array_stats["antisense_midpoint_fragments"],
                "sense_rpk": f"{array_rpk:.12g}",
                "sense_tpm_259_features": f"{(array_rpk * 1_000_000 / tpm_denominator) if tpm_denominator else 0:.12g}",
            }
        )
        for feature, sense, antisense, rpk in zip(features, feature_sense, feature_antisense, feature_rpk):
            writer.writerow(
                {
                    "run_accession": run,
                    "sample_alias": run_meta["sample_alias"],
                    "time_min": time_min,
                    "culture": culture,
                    "feature_id": feature["feature_id"],
                    "feature_type": "CDS",
                    "protein_id": feature["protein_id"],
                    "locus_tag": feature["locus_tag"],
                    "product": feature["product"],
                    "start0": feature["start0"],
                    "end0": feature["end0"],
                    "length_nt": feature["length_nt"],
                    "strand": feature["strand"],
                    "sense_midpoint_fragments": sense,
                    "antisense_midpoint_fragments": antisense,
                    "sense_rpk": f"{rpk:.12g}",
                    "sense_tpm_259_features": f"{(rpk * 1_000_000 / tpm_denominator) if tpm_denominator else 0:.12g}",
                }
            )

    with fastp_json.open(encoding="utf-8") as handle:
        fastp_report = json.load(handle)
    bowtie_text = bowtie_log.read_text(encoding="utf-8", errors="replace")
    alignment_rate_match = re.search(r"([0-9.]+)% overall alignment rate", bowtie_text)
    elapsed_seconds = time.perf_counter() - started_clock
    status.update(
        {
            "status": "complete",
            "finished_utc": iso_now(),
            "elapsed_seconds": elapsed_seconds,
            "fastp_returncode": fastp_return,
            "bowtie2_returncode": bowtie_return,
            "findstr_returncode": 0,
            "sam_filter_mode": "inline exact RNAME match to MW218148.1",
            "parser_reached_eof": True,
            "bowtie2_overall_alignment_rate_percent": float(alignment_rate_match.group(1)) if alignment_rate_match else None,
            "fastp_summary": fastp_report.get("summary", {}),
            "parser_counters": dict(counters),
            "array_sense_midpoint_fragments": array_stats["sense_midpoint_fragments"],
            "array_share_of_retained_phage_fragments": (
                array_stats["sense_midpoint_fragments"] / retained_phage if retained_phage else 0
            ),
            "output_files": {
                path.name: {"bytes": path.stat().st_size, "sha256": file_digest(path, "sha256")}
                for path in (coverage_path, window_path, feature_path, fastp_json, fastp_html, fastp_log, bowtie_log)
            },
        }
    )
    write_status(status_path, status)
    print(
        json.dumps(
            {
                "run_accession": run,
                "status": "complete",
                "elapsed_seconds": elapsed_seconds,
                "fastp_after_reads": fastp_report.get("summary", {}).get("after_filtering", {}).get("total_reads"),
                "bowtie2_overall_alignment_rate_percent": status["bowtie2_overall_alignment_rate_percent"],
                "retained_phage_fragments": retained_phage,
                "array_sense_midpoint_fragments": array_stats["sense_midpoint_fragments"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
