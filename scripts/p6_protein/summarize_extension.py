#!/usr/bin/env python3
"""Build a local, machine-readable ledger for the frozen P6 extension search.

No network requests are made.  The script summarizes the exact records already
retrieved by ``retrieve_extension.py`` and preserves overlapping exclusion
reasons rather than forcing each protein into a single biological class.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SEARCHED = ROOT / "data" / "processed" / "p6_extension" / "searched_proteins.tsv"
PREPASS = ROOT / "data" / "processed" / "p6_extension" / "rt_prepass.tsv"
MANIFEST = ROOT / "data" / "processed" / "p6_extension" / "extension_manifest.tsv"
FIXED = ROOT / "data" / "processed" / "p2_loci" / "fixed_10_loci.tsv"
QUERY_LEDGER = ROOT / "data" / "source_metadata" / "p6_protein" / "ncbi_protein_queries.tsv"
RECEIPT = ROOT / "data" / "source_metadata" / "p6_protein" / "extension_retrieval_receipt.json"
OUT_LEDGER = ROOT / "data" / "processed" / "p6_extension" / "extension_exclusion_ledger.tsv"
OUT_SUMMARY = ROOT / "data" / "processed" / "p6_extension" / "extension_machine_summary.json"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def false(value: str) -> bool:
    return value.strip().lower() != "true"


def main() -> None:
    searched = read_tsv(SEARCHED)
    prepass = read_tsv(PREPASS)
    manifest = read_tsv(MANIFEST)
    fixed = read_tsv(FIXED)
    queries = read_tsv(QUERY_LEDGER)
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))

    searched_by_id = {r["protein_accession"]: r for r in searched}
    prepass_by_id = {r["protein_accession"]: r for r in prepass}
    manifest_by_id = {r["rt_protein_accession"]: r for r in manifest}
    fixed_by_id = {r["rt_protein_id"]: r["label"] for r in fixed}

    if set(searched_by_id) != set(prepass_by_id):
        raise ValueError("searched_proteins and rt_prepass accessions differ")

    ledger = []
    exclusion_combinations = Counter()
    reason_counts = Counter()
    for protein_id in sorted(searched_by_id):
        source = searched_by_id[protein_id]
        pre = prepass_by_id[protein_id]
        reasons = []
        if false(pre["passes_length"]):
            reasons.append("length_outside_450_600_aa")
        if false(pre["passes_rt_homology"]):
            reasons.append("fixed10_RT_homology_threshold_failed")
        if false(pre["passes_yxdd"]):
            reasons.append("single_C_terminal_350aa_YxDD_rule_failed")
        if false(pre["passes_simple_coded_by"]):
            reasons.append("coded_by_not_simple_resolvable_accession_location")
        row = manifest_by_id.get(protein_id)
        if reasons:
            outcome = "excluded_rt_prepass"
            exclusion_combinations[";".join(reasons)] += 1
            reason_counts.update(reasons)
            neighborhood_status = "not_evaluated_after_prepass_failure"
        elif row:
            outcome = row["inclusion_status"]
            neighborhood_status = row["inclusion_reason"]
        else:
            outcome = "passed_prepass_missing_manifest"
            neighborhood_status = "internal_resolution_gap"
        ledger.append(
            {
                "protein_accession": protein_id,
                "organism": source["organism"],
                "length_aa": source["length_aa"],
                "source_query_ids": source["source_query_ids"],
                "best_fixed_rt": pre["best_fixed_rt"],
                "rt_evalue": pre["rt_evalue"],
                "rt_query_coverage": pre["rt_query_coverage"],
                "rt_subject_coverage": pre["rt_subject_coverage"],
                "yxdd_count_c_terminal_350": pre["yxdd_count_c_terminal_350"],
                "rt_prepass_pass": pre["passes_rt_prepass"],
                "final_outcome": outcome,
                "exclusion_or_neighborhood_reason": ";".join(reasons) if reasons else neighborhood_status,
                "nucleotide_accession": row["nucleotide_accession"] if row else pre["nucleotide_accession"],
                "technical_duplicate_group": row["technical_duplicate_group"] if row else "",
                "technical_duplicate_aliases": row["technical_duplicate_aliases"] if row else "",
                "is_fixed10_RT": str(protein_id in fixed_by_id).lower(),
                "fixed10_label": fixed_by_id.get(protein_id, ""),
            }
        )

    OUT_LEDGER.parent.mkdir(parents=True, exist_ok=True)
    with OUT_LEDGER.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ledger[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(ledger)

    fixed_recovered = sorted((fixed_by_id[x], x) for x in manifest_by_id if x in fixed_by_id)
    fixed_not_recovered = sorted((label, protein_id) for protein_id, label in fixed_by_id.items() if protein_id not in searched_by_id)
    alias_size_by_group = {}
    for row in manifest:
        group = row["technical_duplicate_group"]
        if group and row["technical_duplicate_aliases"]:
            alias_size_by_group[group] = len(row["technical_duplicate_aliases"].split(";"))
    alias_sizes = list(alias_size_by_group.values())
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "selection_array_blind": True,
        "source_scope": "five bounded NCBI Protein ESearch queries followed by local fixed10 RT/partner and neighborhood rules",
        "record_counts": {
            "query_reported_nonunique": sum(int(r["reported_count"]) for r in queries),
            "unique_protein_records_retrieved_raw": len(searched),
            "rt_prepass_pass": sum(r["passes_rt_prepass"].lower() == "true" for r in prepass),
            "rt_prepass_excluded": sum(r["passes_rt_prepass"].lower() != "true" for r in prepass),
            "existing_fixed10_recovered": sum(r["inclusion_status"] == "existing_fixed_10_reference" for r in manifest),
            "new_included_extensions": sum(r["inclusion_status"] == "included_extension" for r in manifest),
            "postprepass_other_exclusions": sum(r["inclusion_status"].startswith("excluded") for r in manifest),
        },
        "fixed10_recovered": [{"label": label, "rt_protein_accession": protein_id} for label, protein_id in fixed_recovered],
        "fixed10_not_recovered_by_index_queries": [
            {"label": label, "rt_protein_accession": protein_id} for label, protein_id in fixed_not_recovered
        ],
        "explicit_index_gap": "SA1/QPI16926.1 was not returned by the bounded NCBI Protein index queries; this is an index-coverage observation, not a negative homology result",
        "exclusion_reason_counts_overlapping": dict(sorted(reason_counts.items())),
        "exclusion_reason_combinations": dict(sorted(exclusion_combinations.items())),
        "technical_deduplication": {
            "rule": "collapse only exact complete nucleotide record plus the same RT coordinates and strand; RT90 is descriptive only",
            "groups_evaluated_after_prepass": len(alias_size_by_group),
            "maximum_alias_count_in_group": max(alias_sizes) if alias_sizes else 0,
            "technical_alias_collapses_in_this_result": sum(max(0, size - 1) for size in alias_sizes),
        },
        "query_counts": [
            {
                "query_id": r["query_id"],
                "reported_count": int(r["reported_count"]),
                "returned_count": int(r["returned_count"]),
                "truncated_at_retmax": r["truncated_at_retmax"].lower() == "true",
            }
            for r in queries
        ],
        "interpretation": "The bounded public index search produced no new assessable Type-II ART-family locus. It does not establish absence from nr, GenBank, or unindexed records.",
        "inputs": [str(x.relative_to(ROOT)).replace("\\", "/") for x in (SEARCHED, PREPASS, MANIFEST, FIXED, QUERY_LEDGER, RECEIPT)],
        "outputs": [str(OUT_LEDGER.relative_to(ROOT)).replace("\\", "/")],
    }
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary["record_counts"], indent=2))


if __name__ == "__main__":
    main()
