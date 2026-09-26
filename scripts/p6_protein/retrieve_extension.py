#!/usr/bin/env python3
"""Retrieve and freeze public type-II ART-family neighborhood candidates.

Candidate selection uses RT and direct-partner protein evidence only. Upstream
sequence is extracted after selection fields have been evaluated; no repeat or
array outcome is used for inclusion, ranking, or deduplication.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "tools" / "python312_packages"))

from Bio import SeqIO  # noqa: E402
from Bio.Seq import Seq  # noqa: E402


QUERIES = [
    ("Q1_ST85_VIRAL", 'Viruses[Organism] AND retron_St85_RT[Protein Name]'),
    ("Q1B_ST85_VIRAL_ALL_FIELDS", 'Viruses[Organism] AND retron_St85_RT[All Fields]'),
    (
        "Q2_STAPH_RNA_DNA_POL",
        'Viruses[Organism] AND Staphylococcus[All Fields] AND "RNA-directed DNA polymerase"[Protein Name]',
    ),
    (
        "Q3_STAPH_REVERSE_TRANSCRIPTASE",
        'Viruses[Organism] AND Staphylococcus[All Fields] AND "reverse transcriptase"[Protein Name]',
    ),
    ("Q4_STAPH_RETRON", 'Viruses[Organism] AND Staphylococcus[All Fields] AND retron[All Fields]'),
]

FIXED_RT_FASTA = PROJECT / "data" / "processed" / "p1_adjust_assignment" / "rt_sequences.faa"
FIXED_PARTNER_FASTA = (
    PROJECT / "data" / "processed" / "p1_adjust_assignment" / "direct_downstream_partner_sequences.faa"
)
FIXED_10 = PROJECT / "data" / "processed" / "p2_loci" / "fixed_10_loci.tsv"
P1_PROTEIN_DIR = PROJECT / "data" / "raw" / "p1" / "proteins"
RAW_ROOT = PROJECT / "data" / "raw" / "p6_protein"
OUT_ROOT = PROJECT / "data" / "processed" / "p6_extension"
META_ROOT = PROJECT / "data" / "source_metadata" / "p6_protein"
BLAST_BIN = (
    PROJECT
    / "scripts"
    / "p1_data"
    / "tools"
    / "ncbi-blast-2.17.0+"
    / "bin"
    / "blastp.exe"
)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "ART-P6-public-locus-extension/1.0"


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def request(url: str, attempts: int = 3) -> tuple[bytes, dict]:
    last = None
    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=60) as response:
                payload = response.read()
                meta = {
                    "url": url,
                    "http_status": response.status,
                    "bytes": len(payload),
                    "attempt": attempt,
                    "retrieved_utc": now_utc(),
                }
                time.sleep(0.36)
                return payload, meta
        except Exception as exc:  # bounded network retry, captured in receipt
            last = exc
            if attempt < attempts:
                time.sleep(1.5 * attempt)
    raise RuntimeError("request failed after {} attempts: {} :: {}".format(attempts, url, last))


def eutils_url(endpoint: str, params: dict) -> str:
    return EUTILS + "/" + endpoint + "?" + urllib.parse.urlencode(params)


def read_fasta(path: Path) -> dict[str, str]:
    rows = {}
    for record in SeqIO.parse(str(path), "fasta"):
        rows[record.id] = str(record.seq).upper()
    return rows


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for name, sequence in records:
            handle.write(">{}\n".format(name))
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if not rows and fields is None:
        raise ValueError("fields required for empty TSV: {}".format(path))
    fields = fields or list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def sanitize(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def root_accession(value: str) -> str:
    return value.split(".", 1)[0]


def feature_qual(feature, name: str) -> str:
    values = feature.qualifiers.get(name, [])
    return values[0] if values else ""


def coded_by(record) -> str:
    for feature in record.features:
        if feature.type == "CDS":
            value = feature_qual(feature, "coded_by")
            if value:
                return value
    return ""


def parse_coded_by(value: str) -> tuple[str, int, int, int] | None:
    if not value or "join(" in value or "order(" in value:
        return None
    strand = -1 if value.startswith("complement(") else 1
    match = re.search(r"([A-Z][A-Z0-9_]*\.\d+):<?(\d+)\.\.>?(\d+)", value)
    if not match:
        return None
    return match.group(1), int(match.group(2)), int(match.group(3)), strand


def protein_regions(record) -> str:
    output = []
    for feature in record.features:
        if feature.type != "Region":
            continue
        name = feature_qual(feature, "region_name") or "unnamed"
        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        dbxref = ";".join(feature.qualifiers.get("db_xref", []))
        output.append("{}:{}-{}:{}".format(name, start, end, dbxref))
    return "|".join(output)


def run_blast(query: Path, subject: Path, output: Path) -> list[dict]:
    fields = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "qlen",
        "slen",
        "evalue",
        "bitscore",
        "qcovhsp",
        "qstart",
        "qend",
        "sstart",
        "send",
    ]
    command = [
        str(BLAST_BIN),
        "-query",
        str(query),
        "-subject",
        str(subject),
        "-seg",
        "no",
        "-max_target_seqs",
        "20",
        "-outfmt",
        "6 " + " ".join(fields),
        "-out",
        str(output),
    ]
    completed = subprocess.run(command, cwd=str(PROJECT), capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError("BLASTP failed: {}".format(completed.stderr.strip()))
    rows = []
    if output.exists():
        with output.open("r", encoding="utf-8") as handle:
            for line in handle:
                values = line.rstrip("\n").split("\t")
                rows.append(dict(zip(fields, values)))
    return rows


def best_hits(rows: list[dict]) -> dict[str, dict]:
    output = {}
    for row in rows:
        row = dict(row)
        for field in ("pident", "evalue", "bitscore", "qcovhsp"):
            row[field] = float(row[field])
        for field in ("length", "qlen", "slen", "qstart", "qend", "sstart", "send"):
            row[field] = int(row[field])
        row["subject_coverage"] = row["length"] / float(row["slen"])
        current = output.get(row["qseqid"])
        if current is None or (row["bitscore"], -row["evalue"], row["sseqid"]) > (
            current["bitscore"],
            -current["evalue"],
            current["sseqid"],
        ):
            output[row["qseqid"]] = row
    return output


def strict_yxdd(sequence: str) -> list[tuple[int, str]]:
    threshold = max(0, len(sequence) - 350)
    return [
        (match.start() + 1, match.group(0))
        for match in re.finditer(r"Y.DD", sequence)
        if match.start() >= threshold
    ]


def feature_interval(feature) -> tuple[int, int, int]:
    return int(feature.location.start), int(feature.location.end), int(feature.location.strand or 0)


def find_rt_feature(record, protein_accession: str, protein_sequence: str):
    root = root_accession(protein_accession)
    exact = []
    translated = []
    for feature in record.features:
        if feature.type != "CDS":
            continue
        protein_id = feature_qual(feature, "protein_id")
        translation = feature_qual(feature, "translation").replace(" ", "").upper()
        if protein_id and root_accession(protein_id) == root:
            exact.append(feature)
        elif translation == protein_sequence:
            translated.append(feature)
    matches = exact or translated
    return matches[0] if len(matches) == 1 else None


def first_downstream(record, rt):
    rt_start, rt_end, rt_strand = feature_interval(rt)
    candidates = []
    for feature in record.features:
        if feature.type != "CDS" or feature is rt:
            continue
        start, end, strand = feature_interval(feature)
        if rt_strand == 1 and start >= rt_end:
            candidates.append((start - rt_end, start, end, strand, feature))
        elif rt_strand == -1 and end <= rt_start:
            candidates.append((rt_start - end, -end, start, strand, feature))
    return min(candidates, key=lambda item: (item[0], item[1], item[2])) if candidates else None


def oriented_upstream(record, rt, requested: int = 6000) -> tuple[str, int]:
    sequence = str(record.seq).upper()
    start, end, strand = feature_interval(rt)
    topology = str(record.annotations.get("topology", "")).lower()
    if strand == 1:
        if topology == "circular" and len(sequence) >= requested:
            doubled = sequence + sequence
            return doubled[start + len(sequence) - requested : start + len(sequence)], requested
        available = min(requested, start)
        return sequence[start - available : start], available
    if strand == -1:
        if topology == "circular" and len(sequence) >= requested:
            doubled = sequence + sequence
            return str(Seq(doubled[end : end + requested]).reverse_complement()), requested
        available = min(requested, len(sequence) - end)
        return str(Seq(sequence[end : end + available]).reverse_complement()), available
    return "", 0


def oriented_feature_nt(record, feature) -> str:
    start, end, strand = feature_interval(feature)
    sequence = str(record.seq[start:end]).upper()
    return str(Seq(sequence).reverse_complement()) if strand == -1 else sequence


def main() -> None:
    started = now_utc()
    for path in (RAW_ROOT, OUT_ROOT, META_ROOT, RAW_ROOT / "protein_batches", RAW_ROOT / "genomes"):
        path.mkdir(parents=True, exist_ok=True)

    request_receipts = []
    query_rows = []
    uid_to_queries = defaultdict(set)
    for query_id, term in QUERIES:
        url = eutils_url(
            "esearch.fcgi",
            {"db": "protein", "term": term, "retmax": 500, "retmode": "json"},
        )
        payload, meta = request(url)
        request_receipts.append({"query_id": query_id, "kind": "protein_esearch", **meta})
        raw_path = RAW_ROOT / "{}_esearch.json".format(query_id)
        raw_path.write_bytes(payload)
        data = json.loads(payload.decode("utf-8"))
        result = data["esearchresult"]
        ids = result.get("idlist", [])
        count = int(result.get("count", 0))
        for uid in ids:
            uid_to_queries[uid].add(query_id)
        query_rows.append(
            {
                "query_id": query_id,
                "query": term,
                "reported_count": count,
                "returned_count": len(ids),
                "retmax": 500,
                "truncated_at_retmax": str(count > len(ids)).lower(),
                "raw_path": raw_path.relative_to(PROJECT).as_posix(),
                "request_url": url,
                "http_status": meta["http_status"],
            }
        )
    write_tsv(META_ROOT / "ncbi_protein_queries.tsv", query_rows)

    protein_records = {}
    protein_sources = {}
    accession_queries = defaultdict(set)
    uids = sorted(uid_to_queries, key=lambda value: int(value))
    for batch_index, offset in enumerate(range(0, len(uids), 50), 1):
        batch = uids[offset : offset + 50]
        summary_url = eutils_url(
            "esummary.fcgi",
            {"db": "protein", "id": ",".join(batch), "retmode": "json"},
        )
        summary_payload, summary_meta = request(summary_url)
        request_receipts.append({"query_id": "UNION", "kind": "protein_esummary", **summary_meta})
        summary_path = RAW_ROOT / "protein_batches" / "protein_batch_{:03d}_esummary.json".format(batch_index)
        summary_path.write_bytes(summary_payload)
        summary = json.loads(summary_payload.decode("utf-8")).get("result", {})
        for uid in batch:
            accession = summary.get(uid, {}).get("accessionversion", "")
            if accession:
                accession_queries[accession].update(uid_to_queries[uid])
        url = eutils_url(
            "efetch.fcgi",
            {"db": "protein", "id": ",".join(batch), "rettype": "gp", "retmode": "text"},
        )
        payload, meta = request(url)
        request_receipts.append({"query_id": "UNION", "kind": "protein_efetch", **meta})
        path = RAW_ROOT / "protein_batches" / "protein_batch_{:03d}.gp".format(batch_index)
        path.write_bytes(payload)
        for record in SeqIO.parse(str(path), "genbank"):
            accession = record.id
            protein_records[accession] = record
            protein_sources[accession] = path.relative_to(PROJECT).as_posix()

    fixed_rt = read_fasta(FIXED_RT_FASTA)
    fixed_partner = read_fasta(FIXED_PARTNER_FASTA)
    fixed_roots = {root_accession(value) for value in fixed_rt}

    protein_meta_rows = []
    search_rt_records = []
    for accession, record in sorted(protein_records.items()):
        sequence = str(record.seq).upper()
        search_rt_records.append((accession, sequence))
        protein_meta_rows.append(
            {
                "protein_accession": accession,
                "description": record.description,
                "organism": record.annotations.get("organism", ""),
                "length_aa": len(sequence),
                "coded_by": coded_by(record),
                "cdd_regions": protein_regions(record),
                "raw_path": protein_sources[accession],
                "source_query_ids": ";".join(sorted(accession_queries.get(accession, {"UNION_UNMAPPED"}))),
                "source_url": "https://www.ncbi.nlm.nih.gov/protein/{}".format(accession),
            }
        )
    write_tsv(OUT_ROOT / "searched_proteins.tsv", protein_meta_rows)
    searched_rt_fasta = OUT_ROOT / "searched_proteins.faa"
    fixed_rt_copy = OUT_ROOT / "fixed_10_rt_references.faa"
    fixed_partner_copy = OUT_ROOT / "fixed_10_partner_references.faa"
    write_fasta(searched_rt_fasta, search_rt_records)
    write_fasta(fixed_rt_copy, sorted(fixed_rt.items()))
    write_fasta(fixed_partner_copy, sorted(fixed_partner.items()))

    rt_blast_path = OUT_ROOT / "searched_vs_fixed_rt_blastp.tsv"
    rt_blast_rows = run_blast(searched_rt_fasta, fixed_rt_copy, rt_blast_path)
    rt_best = best_hits(rt_blast_rows)

    prepass = []
    protein_meta_by_id = {row["protein_accession"]: row for row in protein_meta_rows}
    for accession, sequence in search_rt_records:
        best = rt_best.get(accession)
        motifs = strict_yxdd(sequence)
        mapping = parse_coded_by(protein_meta_by_id[accession]["coded_by"])
        pass_length = 450 <= len(sequence) <= 600
        pass_homology = bool(
            best
            and best["evalue"] <= 1e-10
            and best["qcovhsp"] / 100.0 >= 0.70
            and best["subject_coverage"] >= 0.70
        )
        pass_motif = len(motifs) == 1
        pass_mapping = mapping is not None
        prepass.append(
            {
                "protein_accession": accession,
                "length_aa": len(sequence),
                "best_fixed_rt": best["sseqid"] if best else "",
                "rt_evalue": best["evalue"] if best else "",
                "rt_bitscore": best["bitscore"] if best else "",
                "rt_percent_identity": best["pident"] if best else "",
                "rt_query_coverage": best["qcovhsp"] / 100.0 if best else "",
                "rt_subject_coverage": best["subject_coverage"] if best else "",
                "yxdd_count_c_terminal_350": len(motifs),
                "yxdd_positions_1based": ";".join(str(item[0]) for item in motifs),
                "yxdd_motifs": ";".join(item[1] for item in motifs),
                "coded_by": protein_meta_by_id[accession]["coded_by"],
                "nucleotide_accession": mapping[0] if mapping else "",
                "passes_length": str(pass_length).lower(),
                "passes_rt_homology": str(pass_homology).lower(),
                "passes_yxdd": str(pass_motif).lower(),
                "passes_simple_coded_by": str(pass_mapping).lower(),
                "passes_rt_prepass": str(pass_length and pass_homology and pass_motif and pass_mapping).lower(),
            }
        )
    write_tsv(OUT_ROOT / "rt_prepass.tsv", prepass)

    passing_accessions = [row["protein_accession"] for row in prepass if row["passes_rt_prepass"] == "true"]
    genome_accessions = sorted({row["nucleotide_accession"] for row in prepass if row["passes_rt_prepass"] == "true"})
    nucleotide_records = {}
    nucleotide_paths = {}
    for accession in genome_accessions:
        reuse = None
        for candidate in (
            PROJECT / "data" / "raw" / "ncbi" / (accession + ".gb"),
            PROJECT / "data" / "raw" / "p1" / "genomes" / (accession + ".gb"),
        ):
            if candidate.exists():
                reuse = candidate
                break
        if reuse is None:
            url = eutils_url(
                "efetch.fcgi",
                {"db": "nuccore", "id": accession, "rettype": "gbwithparts", "retmode": "text"},
            )
            payload, meta = request(url)
            request_receipts.append({"query_id": "RT_PREPASS", "kind": "nuccore_efetch", **meta})
            reuse = RAW_ROOT / "genomes" / (sanitize(accession) + ".gb")
            reuse.write_bytes(payload)
        try:
            record = SeqIO.read(str(reuse), "genbank")
        except Exception:
            continue
        nucleotide_records[accession] = record
        nucleotide_paths[accession] = reuse.relative_to(PROJECT).as_posix()

    partner_candidates = []
    locus_work = []
    prepass_by_id = {row["protein_accession"]: row for row in prepass}
    protein_seq = {accession: sequence for accession, sequence in search_rt_records}
    for protein_accession in passing_accessions:
        pre = prepass_by_id[protein_accession]
        nucleotide_accession = pre["nucleotide_accession"]
        record = nucleotide_records.get(nucleotide_accession)
        if record is None:
            locus_work.append({"protein_accession": protein_accession, "resolution_error": "nucleotide_record_unavailable"})
            continue
        rt = find_rt_feature(record, protein_accession, protein_seq[protein_accession])
        if rt is None:
            locus_work.append({"protein_accession": protein_accession, "resolution_error": "rt_feature_not_unique"})
            continue
        downstream = first_downstream(record, rt)
        if downstream is None:
            locus_work.append({"protein_accession": protein_accession, "resolution_error": "no_downstream_cds"})
            continue
        gap, _, _, partner_strand, partner = downstream
        rt_start, rt_end, rt_strand = feature_interval(rt)
        partner_start, partner_end, _ = feature_interval(partner)
        partner_sequence = feature_qual(partner, "translation").replace(" ", "").upper()
        partner_accession = feature_qual(partner, "protein_id") or (
            nucleotide_accession + ":{}-{}".format(partner_start + 1, partner_end)
        )
        upstream, upstream_available = oriented_upstream(record, rt, 6000)
        boundary = (
            "complete_record_assessable"
            if upstream_available >= 6000
            else "partial_record_assessable"
            if upstream_available >= 3000
            else "boundary_incomplete"
        )
        candidate_id = "{}__{}-{}{}__{}".format(
            nucleotide_accession.replace(".", "_"),
            rt_start + 1,
            rt_end,
            "+" if rt_strand == 1 else "-",
            protein_accession.replace(".", "_"),
        )
        partner_query_id = candidate_id + "|" + partner_accession.replace(" ", "_")
        if partner_sequence:
            partner_candidates.append((partner_query_id, partner_sequence))
        locus_work.append(
            {
                "candidate_id": candidate_id,
                "protein_accession": protein_accession,
                "nucleotide_accession": nucleotide_accession,
                "record": record,
                "rt": rt,
                "rt_sequence": protein_seq[protein_accession],
                "rt_start": rt_start,
                "rt_end": rt_end,
                "rt_strand": rt_strand,
                "partner": partner,
                "partner_query_id": partner_query_id,
                "partner_accession": partner_accession,
                "partner_sequence": partner_sequence,
                "partner_start": partner_start,
                "partner_end": partner_end,
                "partner_strand": partner_strand,
                "gap": gap,
                "upstream": upstream,
                "upstream_available": upstream_available,
                "boundary": boundary,
                "raw_path": nucleotide_paths[nucleotide_accession],
            }
        )

    partner_query_fasta = OUT_ROOT / "resolved_direct_partners.faa"
    write_fasta(partner_query_fasta, partner_candidates)
    partner_blast_path = OUT_ROOT / "resolved_vs_fixed_partner_blastp.tsv"
    partner_blast_rows = run_blast(partner_query_fasta, fixed_partner_copy, partner_blast_path) if partner_candidates else []
    partner_best = best_hits(partner_blast_rows)

    fixed_rows = list(csv.DictReader(FIXED_10.open("r", encoding="utf-8"), delimiter="\t"))
    existing_proteins = {row["rt_protein_id"] for row in fixed_rows}
    manifest = []
    fasta_upstream = []
    fasta_rt = []
    fasta_partner = []
    fasta_rt_cds = []
    fasta_partner_cds = []
    for item in locus_work:
        if "resolution_error" in item:
            pre = prepass_by_id[item["protein_accession"]]
            manifest.append(
                {
                    "candidate_id": "UNRESOLVED_" + item["protein_accession"].replace(".", "_"),
                    "rt_protein_accession": item["protein_accession"],
                    "nucleotide_accession": pre.get("nucleotide_accession", ""),
                    "inclusion_status": "excluded_resolution",
                    "inclusion_reason": item["resolution_error"],
                }
            )
            continue
        best = partner_best.get(item["partner_query_id"])
        partner_len = len(item["partner_sequence"])
        same_strand = item["partner_strand"] == item["rt_strand"]
        pass_partner_architecture = same_strand and 0 <= item["gap"] <= 50 and 150 <= partner_len <= 400
        pass_partner_homology = bool(
            best
            and best["evalue"] <= 1e-5
            and best["qcovhsp"] / 100.0 >= 0.50
            and best["subject_coverage"] >= 0.50
        )
        pass_boundary = item["upstream_available"] >= 3000
        is_existing = item["protein_accession"] in existing_proteins
        included = pass_partner_architecture and pass_partner_homology and pass_boundary and not is_existing
        status = "existing_fixed_10_reference" if is_existing else "included_extension" if included else "excluded_neighborhood_rule"
        failures = []
        if not pass_partner_architecture:
            failures.append("direct_partner_architecture")
        if not pass_partner_homology:
            failures.append("direct_partner_homology")
        if not pass_boundary:
            failures.append("upstream_boundary")
        if is_existing:
            failures.append("already_in_fixed_10")
        record = item["record"]
        pre = prepass_by_id[item["protein_accession"]]
        digest = hashlib.sha256(str(record.seq).upper().encode("ascii")).hexdigest()
        row = {
            "candidate_id": item["candidate_id"],
            "label": record.annotations.get("organism", ""),
            "nucleotide_accession": item["nucleotide_accession"],
            "record_definition": record.description,
            "record_length_nt": len(record.seq),
            "record_topology": record.annotations.get("topology", ""),
            "raw_path": item["raw_path"],
            "record_source_url": "https://www.ncbi.nlm.nih.gov/nuccore/{}".format(item["nucleotide_accession"]),
            "rt_protein_accession": item["protein_accession"],
            "rt_protein_source_url": protein_meta_by_id[item["protein_accession"]]["source_url"],
            "source_query_ids": protein_meta_by_id[item["protein_accession"]]["source_query_ids"],
            "rt_cdd_regions": protein_meta_by_id[item["protein_accession"]]["cdd_regions"],
            "rt_start_1based": item["rt_start"] + 1,
            "rt_end_1based": item["rt_end"],
            "rt_strand": "+" if item["rt_strand"] == 1 else "-",
            "rt_length_aa": len(item["rt_sequence"]),
            "rt_product": feature_qual(item["rt"], "product"),
            "best_fixed_rt": pre["best_fixed_rt"],
            "rt_blast_evalue": pre["rt_evalue"],
            "rt_blast_bitscore": pre["rt_bitscore"],
            "rt_percent_identity": pre["rt_percent_identity"],
            "rt_query_coverage": pre["rt_query_coverage"],
            "rt_subject_coverage": pre["rt_subject_coverage"],
            "yxdd_position_1based": pre["yxdd_positions_1based"],
            "yxdd_motif": pre["yxdd_motifs"],
            "partner_protein_accession": item["partner_accession"],
            "partner_start_1based": item["partner_start"] + 1,
            "partner_end_1based": item["partner_end"],
            "partner_strand": "+" if item["partner_strand"] == 1 else "-",
            "partner_length_aa": partner_len,
            "partner_product": feature_qual(item["partner"], "product"),
            "rt_partner_gap_nt": item["gap"],
            "rt_partner_same_strand": str(same_strand).lower(),
            "best_fixed_partner": best["sseqid"] if best else "",
            "partner_blast_evalue": best["evalue"] if best else "",
            "partner_blast_bitscore": best["bitscore"] if best else "",
            "partner_percent_identity": best["pident"] if best else "",
            "partner_query_coverage": best["qcovhsp"] / 100.0 if best else "",
            "partner_subject_coverage": best["subject_coverage"] if best else "",
            "upstream_available_nt": item["upstream_available"],
            "boundary_status": item["boundary"],
            "technical_record_sha256": digest,
            "selection_array_blind": "true",
            "assignment_status": "type_II_ART_family_neighborhood_candidate" if included else "not_new_extension",
            "inclusion_status": status,
            "inclusion_reason": "passes frozen RT+partner+boundary rules" if included else ";".join(failures),
        }
        manifest.append(row)
        if included:
            fasta_upstream.append((item["candidate_id"], item["upstream"]))
            fasta_rt.append((item["candidate_id"] + "|" + item["protein_accession"], item["rt_sequence"]))
            fasta_partner.append((item["candidate_id"] + "|" + item["partner_accession"], item["partner_sequence"]))
            fasta_rt_cds.append((item["candidate_id"], oriented_feature_nt(record, item["rt"])))
            fasta_partner_cds.append((item["candidate_id"], oriented_feature_nt(record, item["partner"])))

    # Exact complete-record plus RT-coordinate technical duplicates only.
    duplicate_groups = defaultdict(list)
    for row in manifest:
        if row.get("technical_record_sha256"):
            key = (
                row["technical_record_sha256"],
                row["rt_start_1based"],
                row["rt_end_1based"],
                row["rt_strand"],
            )
            duplicate_groups[key].append(row)
    for index, members in enumerate(sorted(duplicate_groups.values(), key=lambda rows: rows[0]["candidate_id"]), 1):
        group = "TECH_{:03d}".format(index)
        ordered = sorted(members, key=lambda row: (row["nucleotide_accession"], row["candidate_id"]))
        aliases = ";".join(row["nucleotide_accession"] for row in ordered)
        for rank, row in enumerate(ordered):
            row["technical_duplicate_group"] = group
            row["technical_duplicate_aliases"] = aliases
            row["is_technical_representative"] = str(rank == 0).lower()
            if rank > 0 and row["inclusion_status"] == "included_extension":
                row["inclusion_status"] = "retained_technical_duplicate"
                row["inclusion_reason"] = "exact complete-record and RT-feature duplicate of {}".format(ordered[0]["candidate_id"])

    manifest.sort(key=lambda row: (row.get("nucleotide_accession", ""), row.get("rt_start_1based", ""), row["candidate_id"]))
    manifest_fields = [
        "candidate_id", "label", "nucleotide_accession", "record_definition", "record_length_nt", "record_topology",
        "raw_path", "record_source_url", "rt_protein_accession", "rt_protein_source_url", "source_query_ids", "rt_cdd_regions",
        "rt_start_1based", "rt_end_1based", "rt_strand",
        "rt_length_aa", "rt_product", "best_fixed_rt", "rt_blast_evalue", "rt_blast_bitscore", "rt_percent_identity",
        "rt_query_coverage", "rt_subject_coverage", "yxdd_position_1based", "yxdd_motif", "partner_protein_accession",
        "partner_start_1based", "partner_end_1based", "partner_strand", "partner_length_aa", "partner_product",
        "rt_partner_gap_nt", "rt_partner_same_strand", "best_fixed_partner", "partner_blast_evalue",
        "partner_blast_bitscore", "partner_percent_identity", "partner_query_coverage", "partner_subject_coverage",
        "upstream_available_nt", "boundary_status", "technical_duplicate_group", "technical_duplicate_aliases",
        "is_technical_representative", "selection_array_blind", "assignment_status", "inclusion_status", "inclusion_reason",
    ]
    write_tsv(OUT_ROOT / "extension_manifest.tsv", manifest, manifest_fields)
    write_fasta(OUT_ROOT / "included_upstream_6000nt.fna", fasta_upstream)
    write_fasta(OUT_ROOT / "included_rt_proteins.faa", fasta_rt)
    write_fasta(OUT_ROOT / "included_partner_proteins.faa", fasta_partner)
    write_fasta(OUT_ROOT / "included_rt_cds.fna", fasta_rt_cds)
    write_fasta(OUT_ROOT / "included_partner_cds.fna", fasta_partner_cds)

    receipt = {
        "analysis": "P6 array-blind public type-II ART-family neighborhood extension",
        "started_utc": started,
        "finished_utc": now_utc(),
        "parameter_record": "scripts/p6_protein/ANALYSIS_PARAMETERS.md",
        "query_count": len(QUERIES),
        "query_reported_total_nonunique": sum(row["reported_count"] for row in query_rows),
        "unique_protein_records_retrieved": len(protein_records),
        "rt_prepass_pass": len(passing_accessions),
        "resolved_loci": sum("candidate_id" in row for row in locus_work),
        "new_included_extensions": sum(row.get("inclusion_status") == "included_extension" for row in manifest),
        "existing_fixed_10_recovered": sum(row.get("inclusion_status") == "existing_fixed_10_reference" for row in manifest),
        "excluded_after_prepass": sum(row.get("inclusion_status", "").startswith("excluded") for row in manifest),
        "request_receipts": request_receipts,
        "selection_array_blind": True,
        "blast_executable": str(BLAST_BIN.relative_to(PROJECT)).replace("\\", "/"),
        "outputs": [
            "data/processed/p6_extension/extension_manifest.tsv",
            "data/processed/p6_extension/included_upstream_6000nt.fna",
            "data/processed/p6_extension/included_rt_proteins.faa",
            "data/processed/p6_extension/included_partner_proteins.faa",
        ],
    }
    (META_ROOT / "extension_retrieval_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: receipt[key] for key in (
        "unique_protein_records_retrieved", "rt_prepass_pass", "resolved_loci",
        "new_included_extensions", "existing_fixed_10_recovered", "excluded_after_prepass"
    )}, indent=2))


if __name__ == "__main__":
    main()
