#!/usr/bin/env python3
"""Small synthetic interoperability check; never used as biological evidence."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FASTP = ROOT / "tools" / "p6_rnaseq" / "fastp-1.3.3-windows-ucrt64" / "fastp.exe"
BOWTIE2 = ROOT / "tools" / "p6_rnaseq" / "bowtie2-2.5.5-mingw-x86_64" / "bowtie2-align-s.exe"
INDEX = ROOT / "data" / "processed" / "p6_rnaseq" / "MW218148.1_NZ_CP059679.1"
REFERENCE = ROOT / "data" / "processed" / "p6_rnaseq" / "MW218148.1_NZ_CP059679.1.fna"


def resolve_executable(explicit: str | None, environment_name: str, names: tuple[str, ...], fallback: Path) -> Path:
    candidates = [explicit, os.environ.get(environment_name)]
    candidates.extend(shutil.which(name) for name in names)
    candidates.append(str(fallback))
    for candidate in candidates:
        if candidate and Path(candidate).expanduser().is_file():
            return Path(candidate).expanduser().resolve()
    raise SystemExit(f"Cannot find {names[0]}; use its CLI option, set {environment_name}, or add it to PATH")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def load_first_fasta_sequence(reference: Path) -> str:
    pieces = []
    with reference.open(encoding="ascii") as handle:
        next(handle)
        for line in handle:
            if line.startswith(">"):
                break
            pieces.append(line.strip())
    return "".join(pieces)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastp")
    parser.add_argument("--bowtie2")
    parser.add_argument("--bowtie2-index", type=Path, default=INDEX)
    parser.add_argument("--reference", type=Path, default=REFERENCE)
    parser.add_argument("--tmp", type=Path, default=ROOT / "results/toolchain_smoke_tmp")
    parser.add_argument("--receipt", type=Path, default=ROOT / "results/toolchain_smoke_receipt.json")
    args = parser.parse_args()
    fastp_executable = resolve_executable(args.fastp, "ART_ARRAY_FASTP", ("fastp",), FASTP)
    bowtie2_executable = resolve_executable(
        args.bowtie2, "ART_ARRAY_BOWTIE2", ("bowtie2", "bowtie2-align-s"), BOWTIE2
    )
    index_prefix = args.bowtie2_index if args.bowtie2_index.is_absolute() else ROOT / args.bowtie2_index
    reference = args.reference if args.reference.is_absolute() else ROOT / args.reference
    tmp = args.tmp if args.tmp.is_absolute() else ROOT / args.tmp
    receipt_path = args.receipt if args.receipt.is_absolute() else ROOT / args.receipt

    tmp.mkdir(parents=True, exist_ok=True)
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    sequence = load_first_fasta_sequence(reference)
    read1_path = tmp / "smoke_R1.fastq.gz"
    read2_path = tmp / "smoke_R2.fastq.gz"
    with gzip.open(read1_path, "wt", encoding="ascii", newline="\n") as read1, gzip.open(
        read2_path, "wt", encoding="ascii", newline="\n"
    ) as read2:
        for index in range(20):
            fragment_start = 9500 + index * 3
            fragment_end = fragment_start + 300
            read2_sequence = sequence[fragment_start : fragment_start + 100]
            read1_sequence = reverse_complement(sequence[fragment_end - 100 : fragment_end])
            quality = "I" * 100
            read1.write(f"@smoke_{index}/1\n{read1_sequence}\n+\n{quality}\n")
            read2.write(f"@smoke_{index}/2\n{read2_sequence}\n+\n{quality}\n")

    fastp_command = [
        str(fastp_executable),
        "--in1",
        str(read1_path),
        "--in2",
        str(read2_path),
        "--stdout",
        "--length_required",
        "30",
        "--thread",
        "2",
        "--json",
        str(tmp / "fastp.json"),
        "--html",
        str(tmp / "fastp.html"),
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
        "2",
        "-x",
        str(index_prefix),
        "--interleaved",
        "-",
    ]
    fastp = subprocess.Popen(fastp_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert fastp.stdout is not None
    bowtie = subprocess.Popen(bowtie_command, stdin=fastp.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    fastp.stdout.close()
    assert bowtie.stdout is not None
    bowtie_stdout, bowtie_stderr = bowtie.communicate()
    fastp_stderr = fastp.stderr.read() if fastp.stderr else b""
    bowtie_return = bowtie.returncode
    fastp_return = fastp.wait()

    records = [
        fields
        for line in bowtie_stdout.splitlines()
        if line and not line.startswith(b"@")
        for fields in [line.split(b"\t")]
        if len(fields) >= 3 and fields[2] == b"MW218148.1"
    ]
    proper_read2_plus = sum(
        len(fields) >= 11
        and fields[2] == b"MW218148.1"
        and (int(fields[1]) & 0x2)
        and (int(fields[1]) & 0x80)
        and not (int(fields[1]) & 0x10)
        and int(fields[4]) >= 10
        for fields in records
    )
    receipt = {
        "purpose": "Synthetic tool interoperability only; not a biological or partial-library result.",
        "synthetic_pairs": 20,
        "fastp_returncode": fastp_return,
        "bowtie2_returncode": bowtie_return,
        "sam_filter_mode": "inline exact RNAME match to MW218148.1",
        "phage_sam_records_after_filter": len(records),
        "proper_primary_read2_plus_mapq_ge10": proper_read2_plus,
        "fastp_stderr_tail": fastp_stderr.decode("utf-8", errors="replace")[-1000:],
        "bowtie2_stderr_tail": bowtie_stderr.decode("utf-8", errors="replace")[-1000:],
        "status": "pass"
        if fastp_return == 0 and bowtie_return == 0 and proper_read2_plus == 20
        else "fail",
    }
    with receipt_path.open("w", encoding="utf-8") as handle:
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if receipt["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
