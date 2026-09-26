# RNA ensemble comparison parameters

Frozen on 2026-09-24 before the main structure run. This analysis is a new comparison of the accepted ten-locus annotation and does not alter the earlier copy, unit, or correspondence calls.

## Research question and fixed objects

The primary question is whether the 35 accepted cross-group unit correspondences, especially the distal U1 and U2 correspondences, retain ensemble base-pair-probability (BPP) similarity outside the shared repeat sequence. The objects are operational sequence intervals, not claims about mature RNA endpoints.

- Use the 47 periodic-chain copies and 37 intervening chain units from `data/processed/p2_comparison/`.
- Exclude the four low-confidence edge units from all primary and background calculations.
- Keep the saved 35 cross-group supported pairs fixed. Do not rerun unit discovery or select pairs using structure.
- Use all 217 other cross-group unit combinations within the same 21 locus pairs as the sequence-defined misaligned-unit comparison.
- Treat pairwise units and locus pairs as dependent observations. Summaries are descriptive; no pair-level p-values or covariation tests are used.

## Sequence cuts and masks

Coordinates are oriented toward the RT, with the RT 5-prime nucleotide at zero. The saved calibrated copy start is a group-selected seed register and is not assumed to be the source 14-nt-core start in both groups. The saved `MarsHill14_best_offset_from_copy_start` field gives offset 0 for all 35 reference-group chain copies and offset -2 for all 12 PB50-related chain copies, consistent with the PB50-group `ATGAATACGT` seed starting two nucleotides inside the source `ATATGAATACGTAT` register. Before extracting any structure object, define `core14_anchor_start = calibrated_copy_start + phase_offset`, with phase offset 0 for the reference group and -2 for the PB50-related group. Save this audit explicitly. Following the source report, define the 22-nt operational repeat analog as `[core14_anchor_start-4, core14_anchor_start+18)`.

For every primary unit with left and right phase-corrected 14-mer anchor starts:

1. **Two-repeat context (primary):** `[left_anchor-4, right_anchor+18)`, comprising the left 22-nt repeat, the intervening spacer `[left_anchor+18, right_anchor-4)`, and the right 22-nt repeat.
2. **One-repeat context (boundary sensitivity):** `[left_anchor-4, right_anchor-4)`, comprising the left 22-nt repeat and the same spacer but excluding the right repeat.

These are the only boundary variants. The same spacer is used in both contexts. Extract sequences from `oriented_anchor_to_partner_locus.fna` using the saved oriented-coordinate map. Save RT-relative, oriented-FASTA, and genomic interval coordinates for each object. RNA1–RNA3 from the source study were coverage-defined intervals and are not equated with U1–U4.

## Folding model and stored quantities

- Use the ViennaRNA 2.7.2 Python bindings under project-local Python 3.12.14 at 37 degrees C with default model details. Convert DNA T to RNA U only for folding.
- For each sequence, run MFE followed by partition-function folding with MFE rescaling, then obtain the full BPP matrix, ensemble free energy, centroid structure, and ensemble diversity.
- Store BPP edges at probability at least `1e-6` and per-position pairedness. Classify every pair as repeat-repeat, spacer-spacer, or mixed according to the fixed masks.
- MFE and ensemble-energy values are secondary descriptors. A low MFE or an individual drawing is not treated as functional evidence.

## Structural correspondence metric

Globally align each unit pair with the accepted affine nucleotide scores: match 2, mismatch -1, gap-open -5, and gap-extension -1. Map every nucleotide to its alignment column. For each pair class separately, compare BPP vectors in this shared column-pair coordinate system.

The primary metric is the spacer-spacer BPP overlap in the two-repeat context:

`overlap = 2 * sum_k min(pA_k, pB_k) / (sum_k pA_k + sum_k pB_k)`.

Unaligned sequence positions retain distinct alignment-column keys and therefore contribute pairing mass to the denominator without creating artificial overlap. Secondary metrics are repeat-repeat overlap, mixed overlap, all-pair overlap, and the same metrics in the one-repeat context.

## Two fixed comparisons

1. **Misaligned-unit comparison:** compare the 35 supported pairs with the 217 non-supported cross-group unit combinations. For each of the 18 locus pairs with at least one supported unit pair, also compute the difference between its median supported and median misaligned spacer-spacer overlap.
2. **Dinucleotide-preserving spacer background:** use seed `20260924` and 64 whole-panel replicates. Within each replicate, independently permute the source-defined spacer of every one of the 37 units by a randomized Eulerian trail that exactly preserves sequence length and every adjacent dinucleotide count. Keep both 22-nt repeat sequences unchanged. Reuse one shuffled spacer per unit in both contexts and across all 35 supported comparisons in that panel.

The main null statistic is the median two-repeat spacer-spacer BPP overlap across the fixed 35 supported pairs. Also summarize U1 (15 pairs) and U2 (18 pairs) separately, plus repeat-repeat and mixed overlap. Report the observed value, null median/range/nearest-rank p95, and the number of 64 panels meeting or exceeding the observation. The panel is the replicate; the 35 pairs are not independent tests.

The choice of 64 panels was locked after a single 212-nt ViennaRNA probe required 0.0666 seconds for MFE, partition function, and BPP extraction. The projected 4,810 folds (37 units, two contexts, observed plus 64 panels) fit a single-process run of minutes rather than hours. No replicate count or metric will be changed in response to the result.

## Sequence relationship and positional interpretation

- Relate BPP overlap descriptively to the saved trimmed-spacer global identity using Spearman rho without an independence p-value.
- Report U1, U2, and the two supported U3 pairs separately. Do not infer a general proximal trend from two U3 observations.
- Localize representative agreement or disagreement with BPP/pairedness plots only after the panel summaries are fixed. Do not use a displayed example to select or redefine the main metric.
- Do not interpret predicted structure as proof of transcript boundaries, processing, biochemical function, selection, or ancestral gain/loss.

## Outputs and validation

Write sequence objects, per-unit ensemble summaries, sparse BPP edges, per-position profiles, all 252 cross-group pair scores in both contexts, locus-pair summaries, all null pair scores and panel summaries, and a structured summary under `data/processed/p6_rna/`. Validate exact dinucleotide preservation, unchanged repeats, shared within-panel shuffled objects, fixed object counts, and the 64 completed panels. Produce two or three figures that directly display the resulting effect and its limits.
