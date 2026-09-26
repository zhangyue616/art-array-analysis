# P6 science-diagnostic reproduction

Run commands from the repository or lightweight-package root.

The full workspace can extract the fixed terminal/RT windows from the saved per-library coverage vectors:

```powershell
python -X utf8 scripts/p6_revision_20260925/diagnostics/run_science_diagnostics.py
```

The lightweight package omits those full vectors. Reuse its bounded 2,400-row export explicitly:

```powershell
python -X utf8 scripts/p6_revision_20260925/diagnostics/run_science_diagnostics.py --local-coverage-input data/processed/p6_revision_20260925/s4_terminal_rt_local_coverage_100nt.tsv
```

`--local-coverage-input` does not claim to reopen the original vectors. It requires the 12 run accessions in `whole_array_partitions.tsv`, exactly 200 unique positions per run (terminal 1-based 10547–10646 and RT 10647–10746), matching run metadata, boundary/relative coordinates and region labels, numeric nonnegative depth fields, and a nonempty `source_coverage_file` provenance field. It then recomputes `s4_terminal_rt_continuity_summary.tsv` with the same aggregation used by the full-vector mode.

Both modes use only saved tables. Neither reruns read processing, alignment, RNA folding, or correspondence selection.
