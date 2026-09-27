# Reproducibility guide

This repository separates a short, offline saved-table reproduction from the longer component analyses and from raw-read processing. The levels answer different questions and should not be treated as interchangeable evidence.

## 1. Lightweight offline reproduction

The recommended route uses Python 3.12 to recompute key counts, medians, Spearman correlations, culture-paired changes, and the descriptive Table 1 rule comparison from supplied TSV files, then renders the requested figures. It uses no network access and no randomness.

Create an isolated environment from the repository root:

```bash
python -m venv .venv
```

Activate it with `.venv\Scripts\Activate.ps1` in Windows PowerShell or `source .venv/bin/activate` in a POSIX shell, then install the lightweight requirements and run:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python scripts/reproduce.py --output results/reproduced
```

The command writes:

- `results/reproduced/key_statistics.json`
- `results/reproduced/tables/structure_identity_correlations.tsv`
- `results/reproduced/tables/primary_temporal_changes.tsv`
- `results/reproduced/tables/rule_comparison.tsv`
- `results/reproduced/tables/rule_comparison_pairs.tsv`
- `results/reproduced/tables/rule_comparison_receipt.json`
- `results/reproduced/figures/Figure_1.{pdf,svg,png}` through `Figure_9.{pdf,svg,png}`
- `results/reproduced/figures/Figure_S1.{pdf,svg,png}` and `Figure_S2.{pdf,svg,png}`
- `results/reproduced/environment_versions.json`
- `results/reproduced/reproduction_summary.json`

The 14.60-second, 35-file receipt associated with version 1.0.0 remains
historical evidence for that release only. Version 1.1.0 adds Figure 1 and the
Table 1 outputs; a fresh `reproduction_summary.json` is the execution evidence
for a new invocation. Timings are hardware-specific and exclude environment
creation, package installation, and the raw RNA-seq route.

Useful bounded variants are:

```bash
# Validate required paths without writing results.
python scripts/reproduce.py --check-inputs

# Recompute tables and render only current Figures 6 and 8.
python scripts/reproduce.py --output results/reproduced --figures key

# Recompute tables without rendering figures.
python scripts/reproduce.py --output results/reproduced --figures none
```

`reproduction_summary.json` records the required inputs, consistency checks, generated files, and analyses intentionally not run. `environment_versions.json` records the Python, platform, NumPy, Matplotlib, Biopython, and ViennaRNA versions visible to that invocation. The presence of supplied outputs alone is not evidence of a fresh successful run; use the newly written summary and environment files for that purpose.

The reference figure exports under `figures/reference/` remain separate from regenerated figures. The lightweight environment pins Matplotlib 3.10.6. The accepted Figure 1 export records Matplotlib 3.11.2, while the accepted Figure 2–9 and Figure S1–S2 data-figure exports retain their earlier Matplotlib 2.2.3 provenance. The lightweight command checks semantic outputs and file structure; it does not compare pixels or claim identical font metrics across these environments.

## 2. Component-level regeneration

Install the extended Python dependencies before running sequence or folding components:

```bash
python -m pip install -r requirements-full.txt
```

Some component commands overwrite the supplied `data/processed/` tree; output paths are listed below. Use a disposable clone or copy when running an in-place component whose saved outputs must remain untouched.

| Component | Command from repository root | Principal inputs | Principal outputs | Recorded execution evidence |
|---|---|---|---|---|
| Current Figure 1 workflow | `python scripts/p6_revision_20260926_r8/render_workflow.py --project-root . --output-dir results/workflow` | Supplied unit, correspondence, BPP, SA1-window, and temporal tables | `results/workflow/Figure_1.{pdf,svg,png}` | The renderer validates the expected 37 units, 252 cross-group combinations, 35 selected pairs, representative BPP rows, and 12 saved libraries. |
| Descriptive correspondence-rule comparison | `python scripts/p6_revision_20260926_r8/summarize_rule_comparison.py --project-root . --output-dir results/rule-comparison` | Saved P2 best-hit, unit, copy, and supported-pair tables | Rule summary, pair-level evidence, and receipt | Recomputes 61/40/48/35 nonredundant pairs and verifies the full-rule 35-pair set exactly; it does not benchmark truth or accuracy. |
| Ten-locus annotation and correspondence | `python scripts/p2_analysis/run_unified_locus_analysis.py --project-root .` | `data/processed/p2_loci/` | `data/processed/p2_comparison/` | Saved receipt reports 10 loci, 47 primary copies, 37 primary units, and four secondary weak edges; portable wall time was not recorded. |
| Block-shuffled correspondence calibration | `python scripts/p3_calibration/run_spacer_order_calibration.py` | Saved P2 loci, copies, units, and correspondence objects | `data/processed/p3_calibration/` | 64 panels, seed 20260924; saved elapsed time 127.773 s. |
| Ensemble-pairing analysis | `python scripts/p6_rna/run_rna_ensemble_comparison.py` | P2 fixed loci, sequences, copies, units, identities, and 35 supported pairs | `data/processed/p6_rna/` | Python 3.12.14, ViennaRNA 2.7.2, one process; saved elapsed time 686.610 s. |
| Ensemble-output validation | `python scripts/p6_rna/validate_p6_rna_outputs.py` | Regenerated `data/processed/p6_rna/` | `data/processed/p6_rna/validation_receipt.json` | Validates object, pair, panel, and dinucleotide-preservation counts. |
| Fixed-ten protein summaries | `python scripts/p6_protein/analyze_protein_features.py` | Bundled fixed GenPept records and locus tables | `data/processed/p6_protein/` | Saved receipt records Python 3.12.14, Biopython 1.88, Matplotlib 3.10.6, and 45 unordered pairs; portable wall time was not recorded. |
| Bounded-screen summary | `python scripts/p6_protein/summarize_extension.py` | Saved evaluated-source ledger and bounded-search tables | `data/processed/p6_extension/extension_exclusion_ledger.tsv` and summary receipt | Summarizes saved retrieval results; it does not repeat remote queries. |
| RNA-seq saved-table aggregation | `python scripts/p6_rnaseq/aggregate_complete_runs.py` | Twelve bundled per-run `run_status.json`, `window_counts.tsv`, and `feature_counts.tsv` files | `data/processed/p6_rnaseq/aggregate/` | Reaggregates retained summaries without FASTQ, fastp, Bowtie2, SAM/BAM, or full coverage vectors. |
| Saved-table diagnostics | `python -X utf8 scripts/p6_revision_20260925/diagnostics/run_science_diagnostics.py --local-coverage-input data/processed/p6_revision_20260925/s4_terminal_rt_local_coverage_100nt.tsv` | Saved comparison, ensemble, temporal, partition, and bounded 100-nt coverage tables | `data/processed/p6_revision_20260925/` | Local-coverage mode validates 12 runs and 200 unique positions per run; it does not reopen full coverage vectors. |

Four large ensemble intermediates are intentionally absent from the lightweight data snapshot: `unit_bpp_edges.tsv`, `unit_position_profiles.tsv`, `dinucleotide_null_pair_scores.tsv`, and `dinucleotide_shuffle_manifest.tsv`. The ensemble command regenerates them.

The editable equation definitions are independent of manuscript files. Regenerate the supplied MathML JSON with:

```bash
python scripts/p6_revision_20260926_r8/define_equations.py
```

`native_equations.py` is an optional Windows Word helper, not part of the offline route. It requires `python-docx`, `lxml`, and a local Microsoft Office `MML2OMML.XSL`; pass that transform through `xsl_path=` or `ART_ARRAY_MML2OMML_XSL`. The helper inserts a selected equation into a caller-supplied Word document and does not rebuild or distribute a manuscript.

## 3. Raw RNA-seq route and its boundary

The public BioProject metadata file, `data/source_metadata/p6_rnaseq/ena_PRJNA836150_read_run.json`, records twelve paired-end runs, 24 FASTQ files, expected MD5 values, and 16,858,907,171 compressed bytes in total. Raw and trimmed FASTQ files, alignment streams, Bowtie2 indexes, full coverage vectors, tool binaries, and execution logs are not distributed in this repository.

The saved run used native Windows fastp 1.3.3 and Bowtie2 2.5.5 with eight threads. The source workflow used fastp 1.3.6; the saved Windows run's version difference is an explicit implementation limitation. aria2 1.37.0 was used for downloading. External tools retain their own licenses.

The repository already supplies `data/processed/p6_rnaseq/MW218148.1_NZ_CP059679.1.fna`. Build the six-part Bowtie2 index and start the complete downloader/runner with:

```bash
python scripts/p6_rnaseq/build_bowtie2_index.py --threads 8
python scripts/p6_rnaseq/orchestrate_complete_runs.py --threads 8 --force
python scripts/p6_rnaseq/aggregate_complete_runs.py
```

The index builder resolves `bowtie2-build` from `--bowtie2-build`, `ART_ARRAY_BOWTIE2_BUILD`, `PATH`, or the historical bundled-tool location. The orchestrator resolves aria2c, fastp, and Bowtie2 from explicit options, `ART_ARRAY_ARIA2` / `ART_ARRAY_FASTP` / `ART_ARRAY_BOWTIE2`, `PATH`, or historical bundled-tool locations. The repository already includes twelve complete per-library statuses, so the default mode reuses them and reports `processed=0, skipped=12`; it does not download or reprocess them. Use `--force` only in a disposable clone or copy to download and process all twelve libraries sequentially. To stop cleanly after the first complete library for a user-controlled inspection, add `--stop-after-first-library`. A normal resume without `--force` processes only runs not marked complete. Because the other eleven bundled statuses are already complete, continuing a full rerun after a forced first-library check requires `--force` with `--stop-after-first-library` removed; this repeats the first library before running the other eleven.

The executable entry for one already-downloaded run is also available:

```bash
python scripts/p6_rnaseq/run_library.py SRR19152327 \
  --threads 8 \
  --fastp /path/to/fastp \
  --bowtie2 /path/to/bowtie2 \
  --bowtie2-index data/processed/p6_rnaseq/MW218148.1_NZ_CP059679.1 \
  --force
```

`run_library.py` verifies the expected FASTQ byte counts and MD5 values, streams fastp and Bowtie2 output, and writes one set of retained per-library summaries. It expects the two downloaded FASTQ files beneath `data/raw/p6_rnaseq/`, the frozen feature and window tables, and a Bowtie2 index with the supplied prefix. Because the repository includes a complete status for this run, keep `--force` when using the example in a disposable clone or copy; omit it only to reuse the bundled output. The public runner is a portable adaptation: it performs an inline exact RNAME match where the saved Windows workflow used `findstr`. Bundled per-library receipts retain the real `findstr` return code from that original run; new inline-run receipts record the filter mode explicitly and leave `findstr_returncode` null. The inline adaptation has not been exercised in a fresh full raw-data run, so strict equivalence is not claimed.

`scripts/p6_rnaseq/prepare_reference.py` rebuilds the joint FASTA and fixed SA1 tables after exact GenBank files have been placed at `data/raw/ncbi/MW218148.1.gb` and `data/raw/p6_rnaseq/NZ_CP059679.1.gb`. The index command above then consumes that combined FASTA.

Saved per-library receipts report 342.090–551.647 s for processing each library and 4,969.866 s when those twelve elapsed values are summed. These values exclude download time and are not a portable performance benchmark. Peak memory, peak disk use, full download time, and a fresh raw-to-figure run for this public layout were not measured. The 16.86 GB compressed-download total is therefore a lower bound on storage planning, not a peak-disk estimate.

## 4. Other provenance limits

- The fixed-locus sequence analysis begins from supplied oriented sequence objects and coordinate maps. Public-record acquisition and construction of those fixed objects are documented by source metadata but are not exposed as a single raw-to-P2 command here.
- Remote NCBI retrieval scripts are retained for provenance. Public databases can change, and a future query need not reproduce the saved returned records.
- The unfiltered NCBI protein request stopped at a server CPU limit before homology screening completed. Its empty hit list is not biological no-hit evidence.
- The bounded terminal/RT coverage export contains aggregate fragment-span depth. It cannot identify boundary-crossing molecules, a shared transcript, or a mature RNA end.
- The two 64-panel calibrations are conditional analyses with dependent components. The lightweight command reads their saved panel summaries; it does not rerun their randomization or folding stages.

## 5. Inputs, figures, and tables

[`figure_map.tsv`](figure_map.tsv) gives the exact saved inputs, executable command, regenerated output, and fixed reference export for Figure 1–9 and Figure S1–S2. It also maps Table 1 and each key supplementary table to their machine-readable inputs and outputs. Parameter records beside each component script define thresholds, coordinate conventions, random seeds, and interpretation limits.
