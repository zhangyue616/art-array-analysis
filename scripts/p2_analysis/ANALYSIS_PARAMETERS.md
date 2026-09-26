# P2 unified array annotation and locus comparison parameters

Frozen on 2026-09-24 before running the P2 annotation on the ten-locus input
package.  P1 calls and tables remain historical inputs; all P2 annotations are
written separately.  The fixed main set is the seven P0 Staphylococcus phage
loci plus AH12, Machias, and PB50.  LPJP1 is outside this comparison set.

## Coordinate and locus objects

- Orient every locus from the reliable distal coding anchor toward the RT and
  retain the complete anchor-to-RT segment.  Coordinates relative to the RT
  5-prime base use zero at the first RT base and negative values upstream.
- Use the data-lane anchor, RT, direct partner, genomic position map, and
  assignment confidence without re-inferring source or protein assignment.
  The repeat search interval is the physical interval between the proximal end
  of the distal anchor and the RT; coding overlaps inside that interval are
  reported rather than silently removed.
- Analyze the seven-locus seed group and the PB50/AH12/Machias group
  independently first.  Cross-group copy or unit correspondence is reported
  only when sequence alignment supports it; a shared RT-relative ordinal is
  not evidence by itself.

## Main repeat-chain annotation

1. Enumerate every exact 10-nt word occurring at least three times in the
   anchor-to-RT search interval.  Seeds must contain at least three distinct
   bases and no homopolymer run of six or more.  For each seed, collect all
   non-overlapping occurrences with at most one mismatch; when candidate starts
   are less than 10 nt apart, retain the lower-mismatch start and then the lower
   coordinate.
2. Build chains with observed start-to-start gaps of 60--600 nt.  A chain has a
   base period in 120--350 nt, each gap must fit one period or one skipped copy
   (two periods) within 30% of the base period, and at most one skipped slot is
   allowed in a preliminary chain.  At least three observed occurrences are
   required.  The period is refined as the median of gap/multiplier values.
3. Score a chain as `(observed copies - 1) * information content`, where
   information content is the 26-column seed-anchored copy profile relative to
   base composition of that locus search interval.  Rank physical chains by
   this score, then observed-copy count, lower normalized spacing deviation,
   proximity of the proximal copy to the RT, and seed word.  These choices are
   an operationalization of the source report's chain/information-content
   method without reconstructing its shuffle gate.
4. Merge seed views into one physical-chain hypothesis when their copy starts
   have a one-to-one match within 10 nt after a constant shift, or when at
   least 80% of starts match within 10 nt and their refined periods differ by
   at most 10%.  Keep the best-ranked representative and the member seed/view
   audit.  Overlapping low-complexity occurrences do not become separate
   copies.  Each locus receives one top physical-chain annotation; a distinct
   runner-up region remains an explicit boundary alternative rather than a
   second array.
   Selection is then made at group level so that the same physical chain is
   represented in every locus.  A main-chain topology must be supported, after
   medoid-coordinate transfer, by at least `ceil(group_size/2)` loci.  An extra
   edge occurrence below that support remains a locus-specific boundary
   alternative.  Within the winning physical-chain cluster, use the seed phase
   supported by the largest number of group loci, then mean chain score,
   10-mer entropy, and lexical order.  This prevents a chance RT-proximal word
   or a 2--10 nt shifted view of the same repeat from changing copy number or
   unit coordinates in only one locus.
5. Do not discard a start already in the selected periodic seed chain because
   its 26-nt flank profile is weak.  Report seed mismatches, 26-nt consensus
   score, local-alignment score, coding overlap, and group-synteny support for
   every selected start.  This rule specifically prevents a weak edge copy
   from disappearing before gap assessment.
6. Estimate usual spacing from the selected chain after dividing any accepted
   two-period gap by two.  Inspect every internal skipped slot and one expected
   slot beyond each chain edge, once, within plus/minus 40 nt.  A new
   degenerate copy must be non-overlapping, score at least 18/26 to the locus
   profile or contain a seed match within two mismatches, and map within 10 nt
   of the same expected slot in at least two loci of its group.  Choose highest
   profile score, then smallest expected-position offset, then lower
   coordinate.  Added copies do not trigger another search round.

The 18/26 value above is only a floor for a one-pass degenerate candidate at an
already expected, group-supported slot.  It is not the complete array
definition.  Exact/one-mismatch periodic-chain copies remain annotated even
below that flank score.

## Repeat boundary, units, and alignments

- Extract 61 nt around each selected seed start (`[-20,+41)`).  Within each
  group, infer the repeat block as the maximal seed-containing interval whose
  columns have at least 80% consensus, allowing one internal lapse.  Report
  both this data-driven block and the fixed 26-nt `[-6,+20)` comparison window.
- Before stacking these windows, calibrate the seed register by mapping every
  member view to the group-selected seed phase through the whole-locus medoid
  alignment.  Preserve the raw seed word, raw start, and applied integer shift
  in the audit table.  Repeat boundaries and units use the calibrated physical
  start, so seed-word phase differences are not reported as biological length
  differences.
- A core-to-core unit is `[copy_start,next_copy_start)`.  The interval between
  the fixed 26-nt windows is `[left+20,right-6)`.  The P0-comparable trimmed
  spacer is `[left+40,right-26)`, length `period-66`.  These sequences and
  lengths are stored separately.
- Use Biopython `Bio.Align.PairwiseAligner` from the project-local dependency
  directory for reproducible affine global alignments: match 2, mismatch -1,
  gap-open -5, gap-extension -1.  Select the medoid anchor-to-RT locus within
  each group by highest mean pairwise global identity.  Save every medoid pair
  alignment, alignment-coordinate blocks, and mapped anchor/repeat/RT
  positions as the whole-locus coordinate evidence.
- Compare repeat blocks and 20-nt flanks in their seed-anchored orientation;
  compare trimmed spacers and full core-to-core units with the same affine
  global nucleotide alignment.  Free best-hit unit mapping is computed within
  each group before any ordinal label is examined.  Cross-group best hits are
  reported with identity, runner-up margin, and reciprocity, but named as a
  correspondence only if the best is unique, reciprocal, and supported by the
  whole-locus coordinate evidence.

## One alternative/background check

Perform one leave-one-locus-out group check: rebuild the group repeat profile
and medoid-coordinate transfer without the queried locus, then ask whether its
main copy starts are recovered within 10 nt.  Report disagreements as concrete
boundary uncertainty.  Do not run further threshold grids, new shuffles,
phylogenies, or detector comparisons.

## Interpretation and figures

- Copy-count or span differences may be described, but without an ancestral
  state they are not called gain or loss.  Madawaska coding overlap remains an
  annotation conflict, not proof for or against an RNA array.
- Produce two or three 175-mm research figures from P2 tables: the ten-locus
  unified annotation, group-wise whole-locus/copy alignment evidence, and a
  repeat/flank/unit comparison only if it adds a distinct result.  Export PDF,
  SVG, and PNG from the same Matplotlib source.  Missing or unresolved states
  receive explicit symbols; dependent pairwise values do not receive invented
  error bars or p-values.
