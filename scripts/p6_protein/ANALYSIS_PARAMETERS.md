# P6 protein and public-locus extension parameters

Frozen before inspecting any newly retrieved upstream array outcome. The ten accepted type-II loci remain the reference set; P6 extension candidates are selected from RT and coding-neighborhood evidence only.

## Questions

1. Do additional public *Staphylococcus* phage records contain the type-II ART-family RT plus its directly downstream partner in the same orientation, with enough upstream sequence for later array evaluation?
2. Across the fixed ten loci and any frozen extensions, is the catalytic RT core conserved while the long N-terminal domain (NTD) and direct partner carry localized lineage-specific changes that could plausibly alter the interaction context?
3. Does the available variation support a concrete protein feature linked to the array comparison, beyond the partner classes, NTD divergence, and six type-II co-folds already reported by Yoon et al.?

## Public-source search and inclusion

The bounded official NCBI Protein queries are recorded verbatim before retrieval:

- `Viruses[Organism] AND retron_St85_RT[Protein Name]`
- `Viruses[Organism] AND retron_St85_RT[All Fields]` (the broader conserved-domain index used as the bounded P1 fallback)
- `Viruses[Organism] AND Staphylococcus[All Fields] AND "RNA-directed DNA polymerase"[Protein Name]`
- `Viruses[Organism] AND Staphylococcus[All Fields] AND "reverse transcriptase"[Protein Name]`
- `Viruses[Organism] AND Staphylococcus[All Fields] AND retron[All Fields]`

The existing P1/P2 local results from the first query and the closed MarsHill remote-nr search are retained as searched-source evidence, including their negative boundary. New ESearch results are unioned by versioned protein accession. Only records resolvable through official NCBI Protein and Nucleotide links are evaluated further.

### Frozen bounded lineage amendment

The five initial Protein index queries (including the broader St85 conserved-domain field query) recovered no new passing locus beyond the fixed set. Before inspecting any new upstream sequence, the search was therefore extended to complete public nucleotide records from the three viral genera already represented by the ten fixed type-II loci: `Lentusvirus[Organism]`, `Madawaskavirus[Organism]`, and `Machiasvirus[Organism]`, each combined with `complete genome[Title]`. Every annotated CDS in those bounded records is screened against the same fixed RT and partner references and the unchanged numeric rules below. This lineage amendment expands annotation coverage without lowering thresholds or selecting on arrays.

Those three genus queries and the encompassing `Wallmarkvirinae[Organism] AND complete genome[Title]` query resolved only the already fixed public genomes. A final array-blind homology source is therefore frozen: one NCBI BLASTP remote search of nr restricted to viruses, using SA1, MarsHill, and PB50 RTs as three lineage-spanning queries (`E-value <=1e-5`, at most 200 reported targets per query). Returned proteins still must pass the unchanged multi-reference RT, direct-partner, synteny, motif, and boundary rules; the BLAST result alone does not confer ART status. The submitted job is retrieved with bounded waiting and is not resubmitted to obtain a preferred result.

### Recoverable Common URL API source completion

After the first BLAST+ remote CLI process ended without exposing an RID or returning XML, the controller authorized one distinct, recoverable Common URL API submission; this is not recovery of the earlier job. The current official Common URL API table (checked 2026-09-24 at `https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html`) states that only listed parameters are supported. `ENTREZ_QUERY` is not listed, so the new request does **not** claim a reliable server-side virus restriction. It submits the same SA1, MarsHill, and PB50 multi-FASTA once to `nr` with `PROGRAM=blastp`, `EXPECT=1e-5`, `HITLIST_SIZE=200`, `WORD_SIZE=3`, `GAPCOSTS=11 2`, `MATRIX=BLOSUM62`, composition-based statistics `2`, `FILTER=F`, and requested `XML2` output. Any returned record must be source-resolved and post-filtered before applying the unchanged RT, direct-partner, synteny, motif, and boundary rules; no array result enters this filtering. The submitter saves the raw Put response, RID, RTOE, Get URLs, and next allowed check time, then stops without polling. NCBI's current developer guidance (`https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html`) sets a minimum one-minute interval between polls of a single RID.

A public locus is a **type-II ART-family neighborhood candidate** only if all of the following hold before any array scan:

1. A complete RT CDS is present and encodes 450–600 aa.
2. Local BLASTP against the fixed ten RTs finds at least one hit with E-value `<=1e-10`, query coverage `>=0.70`, and subject coverage `>=0.70`. This is an operational P6 bounded-database rule, not the original report's HMM definition.
3. The strict catalytic motif `YxDD` occurs once in the C-terminal 350 aa. Absence is retained as a failure or truncation state, not silently rescued with a looser motif.
4. The first coding sequence downstream in the RT transcription direction is on the same strand, begins within 0–50 nt after the RT, and encodes 150–400 aa.
5. That direct partner has a local BLASTP hit to at least one of the ten fixed type-II partners with E-value `<=1e-5`, query coverage `>=0.50`, and subject coverage `>=0.50`.
6. The nucleotide record supplies the full RT and partner CDSs plus at least 3,000 nt upstream of the RT. Records with a full 6,000-nt RT-oriented upstream window are `complete_record_assessable`; 3,000–5,999 nt are `partial_record_assessable`; shorter or unresolved records are retained in the exclusion ledger. The 6-kb boundary matches the accepted locus workflow and can contain arrays longer than a 3-kb observation window.

Passing these rules supports a homologous type-II ART-family neighborhood. It does not by itself prove an upstream array, biochemical activity, or source-report ART assignment. Array status is evaluated only after the extension manifest is frozen.

## Sampling, redundancy, and ordering

- Exact duplicates require the same complete nucleotide record sequence and the same RT feature coordinates/strand. One stable accession is the representative and all aliases remain mapped.
- Distinct genomes or distinct RT contexts are retained even when their RT proteins are identical. Exact 3-kb upstream contexts are assigned a context group but are not deleted.
- RT clusters at 90% amino-acid identity and 80% length coverage are descriptive lineage strata, not a deduplication rule.
- If more than 30 new assessable neighborhoods pass, representatives are allocated round-robin across RT90 strata, then sorted by boundary status, descending minimum RT/partner coverage, and accession. No array, repeat, motif-count, or RNA-structure result enters this ordering.
- Existing ten loci are never counted as new extensions. Report-named LPJP1 remains a record-level ART reference unless a specific feature also passes the frozen protein/neighborhood rules.

## Protein measurements

- Record full RT and partner lengths, exact YxDD sequence and 1-based position, RT–partner gap, direction, and source annotations.
- Transfer the MarsHill RT-core boundary (MarsHill residue 211) through a saved global affine protein alignment. Residues before the mapped boundary are the operational NTD; the remainder is the operational RT core. Report alignment coverage and do not map across an unresolved gap.
- Protein alignments use a fixed global affine score (match `2`, mismatch `-1`, gap open `-5`, gap extension `-1`). Compute identity as identical residues divided by all alignment columns, separately for full RT, operational NTD, operational core, and direct partner. Pairwise values are dependent summaries, not independent samples.
- Report conserved blocks and group-specific insertions only when they are reproducible in at least two distinct nucleotide-context groups. Single-locus changes remain observations.
- Existing report evidence is background: all 95 ART RTs retain YxDD when complete; type-II partners are approximately 270-aa all-helical proteins without clear database homologs; six type-II RT–partner pairs were already co-folded. P6 does not relabel these as new findings.
- Domain or structure-database results must be accession- and source-traceable. A database absence or low-confidence model is not evidence of no domain or no structure.
- A protein–array association is attempted only if at least three non-identical context groups contain variation in both the protein feature and the post-freeze array outcome. Otherwise the protein result remains descriptive background.

## Outputs and stopping rules

The frozen extension manifest records inclusion and exclusion, source query, versioned accessions, RT/partner coordinates and sequences, homology statistics, upstream boundary status, raw path, and source URL. Each official ESearch query is capped at 500 returned protein records and any larger reported count is preserved as a truncation boundary. Oriented FASTA files contain the RT-direction 6-kb upstream window when available, RT CDS/protein, and direct partner CDS/protein. Retrieval retries are bounded; an unavailable or empty official source is reported, not treated as biological absence. No whole-nr download, exhaustive phylogeny, generic Pfam/Foldseek sweep, or new large structure-prediction batch is run.
