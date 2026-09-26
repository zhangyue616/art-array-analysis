# P3 composition-and-geometry calibration parameters

Frozen on 2026-09-24 before the P3 calibration run.  P2 annotations and
parameters remain unchanged.  This calibration conditions on the P2 chain
annotation and asks whether shuffled spacer order often survives the complete
cross-group selection and coordinate-support workflow.

## Fixed objects

- Use the same ten loci and RT-oriented distal-anchor-through-RT sequences as
  P2.
- Retain only the 47 `participates_in_seed_chain=true` copies and the 37 units
  whose two boundary copies both belong to that chain.  The four weak proximal
  edge units are excluded.  Repeat windows, flanks, distal coding anchor, RT,
  all coordinates, and every base outside the 37 trimmed intervals remain
  fixed.  Existing coding annotations remain attached to their coordinates;
  if a retained Madawaska conflict CDS overlaps a trimmed interval, those DNA
  bases are shuffled with that interval rather than protected by a new coding
  mask.
- For each retained unit, the only shuffled interval is its saved P0-comparable
  trimmed spacer `[left copy start + 40, right copy start - 26)`.  Its saved
  sequence and length must agree exactly with that slice of the locus object.

## Null generation

- Random seed: `20260924`; generator: Python `random.Random` used as one
  continuous deterministic stream.
- Replicates: exactly 64.  Loci are visited in P2 `locus_order`, units in
  distal-to-proximal ordinal order, and blocks from left to right.
- Split each trimmed spacer into consecutive 50-nt blocks; the final shorter
  block is kept as its own block.  Independently permute bases within each
  block.  Thus every block preserves length and exact A/C/G/T counts.
- Generate one complete ten-locus replacement panel per replicate and reuse
  that same panel for all 21 P0-by-PB50-group locus pairs.  Do not shuffle a
  selected pair separately and do not condition on the 35 observed pairs.

## Recomputed selection workflow

- For every replicate, rerun all 21 cross-group distal-anchor-through-RT global
  affine alignments and all 28-by-9=252 cross-group core-to-core and trimmed-
  spacer unit alignments with the P2 settings: match 2, mismatch -1, gap-open
  -5, gap-extension -1, identity denominator including gap columns.
- For each query unit and each locus in the opposite group, rank every target
  unit by core-to-core global identity.  A best hit is unique only if no other
  target is within `1e-12` of its score.  Require the reverse query to select it
  uniquely, then require both copy boundaries from the recomputed whole-locus
  alignment to fall within 10 nt of the selected target boundaries.
- Deduplicate reciprocal directions by the sorted pair of unit IDs.  Ordinal
  labels do not enter ranking or support.

## Prespecified summaries

- Primary statistic: number of deduplicated cross-group supported unit pairs in
  the full replicate.  Observed value is reused from the saved P2 table after
  restricting it to the 37 retained units; no observed alignment is rerun.
- Secondary statistics: among the 21 cross-group locus pairs, coverage by a
  supported same-position `P0 U1 <-> PB50-group U1` pair and by a `U2 <-> U2`
  pair; and the median trimmed-spacer global identity among all deduplicated
  supported pairs in that replicate.  Missing identity is reported if a
  replicate has no supported pair.
- Report observed value, null median, minimum, maximum, and nearest-rank 95th
  percentile.  A panel-level Monte Carlo exceedance fraction may be shown as
  `(1 + number of null replicates >= observed) / 65`; it is descriptive and is
  not a test over 21 locus pairs or individual unit pairs.

This null is conditional on the selected copy geometry, fixed boundaries, and
the two sampled groups.  It calibrates sequence correspondence after local
mononucleotide order is broken; it does not validate repeat discovery, copy
number, or biological array endpoints.  No threshold scan, alternative null,
weak-edge sensitivity rerun, or database expansion is performed.
