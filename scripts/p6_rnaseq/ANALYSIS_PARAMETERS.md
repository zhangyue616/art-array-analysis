# P6 conditional RNA-seq analysis parameters

Frozen before inspecting any per-unit RNA-seq result on 2026-09-24.

## Biological question and closest prior work

The published analysis already shows sense-strand coverage across the SA1 ART locus at 5, 15, and 55 min and reports whole-array, RT, and partner TPM. This analysis asks a narrower question that was not reported there: whether the relative abundance and coverage of the four sequence-defined SA1 core-to-core units (U1-U4) are stable across infection or show a culture-reproducible temporal shift. The 0-min libraries are a background control. Unit-resolved coverage is descriptive evidence about the sequenced RNA population; it does not establish mature RNA boundaries, processing, or function.

## Libraries and culture pairing

The complete paired-end libraries in ENA BioProject PRJNA836150 are used. Culture is the repeated-measures unit.

| Time | Culture 1 | Culture 2 | Culture 3 |
|---|---|---|---|
| 0 min | SRR19152335 | SRR19152334 | SRR19152331 |
| 5 min | SRR19152330 | SRR19152329 | SRR19152328 |
| 15 min | SRR19152327 | SRR19152326 | SRR19152325 |
| 55 min | SRR19152324 | SRR19152333 | SRR19152332 |

The `sample_alias` and `library_name` fields are retained from ENA. The suffix `-1`, `-2`, or `-3` defines culture identity. Reads, bases, genomic positions, and the four units are not treated as biological replicates.

## Reference, preprocessing, and fragment filter

- Reference: joint SA1 phage MW218148.1 and *Staphylococcus lentus* NZ_CP059679.1.
- Preprocessing baseline from the source report: fastp 1.3.6, paired-end mode, minimum retained read length 30 nt. The actual native run uses fastp 1.3.3 with the same explicit minimum-length setting because upstream fastp publishes no Windows 1.3.6 binary and the registered WSL Ubuntu instance failed to mount with `E_ACCESSDENIED`. The Windows build is source-identified and hashed in the tool receipt. This is a declared implementation deviation; no equivalence to 1.3.6 is assumed.
- Alignment: Bowtie2 2.5.5, `--very-sensitive -X 1000 --no-unal`, eight threads.
- The stream is aggregated directly; no complete SAM/BAM or trimmed FASTQ is written.
- A fragment is counted once from its primary read-2 SAM record when it is properly paired, both mates are mapped, MAPQ is at least 10, nonzero absolute TLEN is at most 1,500 nt, and the alignment is neither secondary nor supplementary.
- Fragment strand is the genomic strand of read 2, matching the published statement that read 2 is the sense read. SA1 array, RT, and partner are on the plus strand.
- For same-reference proper pairs, the outer template span is reconstructed from SAM as `start0=min(POS-1,PNEXT-1)` and `end0=start0+abs(TLEN)`. This remains correct when one mate alignment is contained in the other; the current read's CIGAR end is not used to infer the opposite template edge. Midpoint is the floor of this 0-based half-open template-span midpoint. Coverage is the template span intersected with the reference. These implement, respectively, the published midpoint feature assignment and fragment-span coverage concepts.

The published Bowtie2 insert ceiling is 1,000 nt while the downstream retained-template ceiling is 1,500 nt. The latter is retained verbatim as a post-alignment audit, though alignments admitted by the stated Bowtie2 command should already be at most 1,000 nt.

## Frozen SA1 intervals

All internal intervals below are 0-based half-open. The RT 5-prime coordinate is relative position 0 (genomic position 10,647, 1-based). Primary core starts are relative positions -1,049, -837, -666, -465, and -289. For a core start `c`, the repeat mask is `[c-4,c+18)`.

- Published whole array: genomic `[9447,10646)` (positions 9,448-10,646, 1-based inclusive; 1,199 nt).
- RT CDS: genomic `[10646,12152)`.
- Direct partner CDS: genomic `[12153,12951)`.
- Mutually exclusive core-to-core units: U1 `[9597,9809)`, U2 `[9809,9980)`, U3 `[9980,10181)`, U4 `[10181,10357)`.
- Repeat-excluded spacers, defined as `[left_core+18,right_core-4)`: U1 `[9615,9805)`, U2 `[9827,9976)`, U3 `[9998,10177)`, U4 `[10199,10353)`.
- Whole-array count partition: distal residual `[9447,9597)`, U1-U4 as above, and terminal residual `[10357,10646)`. Residual regions are not named as additional units or mature RNAs.

The interval table used by the program is emitted as `data/processed/p6_rnaseq/sa1_windows.tsv`. When `data/processed/p6_rna/unit_sequence_objects.tsv` becomes available, its SA1 seed-phase intervals are compared with these frozen coordinates before the final run. A mismatch stops unit-level aggregation rather than silently changing the windows.

## End points

Primary abundance endpoint per library and unit:

1. Count sense fragments whose midpoint lies in the repeat-excluded spacer.
2. Divide by spacer length to obtain midpoint density.
3. Divide each density by the sum of the four unit densities to obtain the four-part within-chain density share.

Primary temporal contrast is the paired 5-to-55-min change in density share. A unit is called a supported, material temporal change only when all three culture-paired changes have the same nonzero direction, the median absolute share change is at least 0.05, and both libraries in every pair contain at least 100 sense midpoint fragments across the four spacers. The 5-to-15 and 15-to-55 contrasts are secondary and use the same rule. With three cultures, effect sizes and culture-level values are primary; low-powered p-values are not used to turn ambiguous patterns into claims.

Confirmatory unit endpoint: mean sense fragment-span depth in each repeat-excluded spacer and its four-part depth share. Core-to-core midpoint density/share and depth/share are sensitivity summaries because they include repeated sequence. Agreement or disagreement between spacer-only and core-to-core results is reported.

Whole-array context is kept separate:

- midpoint fragments in the distal residual, each core-to-core unit, and the terminal residual as a fraction of all sense midpoint fragments in the 1,199-nt array;
- integrated sense fragment-span coverage for the same partition;
- the fraction of whole-array signal outside U1-U4;
- the published-style array share among all fragments aligned to SA1.

The U1-U4 denominator is never described as the whole array.

## Feature and infection-stage context

Sense midpoint counts and length-normalized TPM-like values are retained for the 258 annotated SA1 CDS features plus the published array interval. The RT and direct partner are reported explicitly. A small set of early/replication and late/structural CDS controls may be selected from existing GenBank product annotations after counting, with the annotation-based selection rule shown. Four time points are not used as independent observations for gene-gene correlation or functional inference.

## Background, estimability, and interpretation

- A post-infection library is unit-estimable when it has at least 100 sense midpoint fragments across the four repeat-excluded spacers. All counts remain in the machine-readable output if this threshold is not met.
- The nominal 0-min libraries quantify baseline/background mapping and are not subtracted from post-infection counts. Unit shares at 0 min are descriptive only when they satisfy the same estimability threshold.
- A stable conclusion requires that no unit meets the predeclared material-change rule and that the largest time-point mean range for every unit is below 0.05; otherwise the result is called indeterminate rather than stable.
- Coverage peaks do not establish RNA-processing boundaries. Relative composition changes do not establish ART function. Similar unit profiles can arise from a longer precursor or shared transcription, and mapping uncertainty over repeated sequence is assessed through MAPQ filtering and spacer-only sensitivity analysis.

## Run completeness and provenance

Every result row records ENA run accession, expected FASTQ bytes and MD5, observed file sizes, fastp before/after read counts, Bowtie2 alignment summary, and fragment-filter counts. A run is complete only when both full FASTQ files pass ENA MD5, fastp and Bowtie2 exit successfully, and the parser reaches EOF. Partial read prefixes are forbidden as inferential results.
