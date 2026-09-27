# Data

This directory contains the public source records, source metadata, and saved derived tables used by the ART-array analysis.

## Directory layout

- `raw/`: versioned third-party NCBI nucleotide and protein records plus the preserved NCBI BLAST URL API terminal response. No RNA-seq FASTQ, SRA, SAM, BAM, Bowtie2 index, executable, or full per-base coverage file is included.
- `source_metadata/`: public archive metadata, retrieval receipts, tool-version records, and source-processing notes. Machine-specific commands, absolute paths, process IDs, local orchestration fields, and omitted-log hashes were removed.
- `processed/`: saved project-derived tables and sequence objects used by the lightweight reproduction entry point and figure renderers.
- `sources.tsv`: source/accession mapping for the ten fixed loci, the RNA-seq host reference and project, and the two source publications.
- `rnaseq_libraries.tsv`: all 12 public RNA-seq runs with culture/time mapping, ENA and SRA pages, direct paired FASTQ HTTPS URLs, expected MD5 values, and byte counts.
- `MANIFEST.tsv`: file inventory, size, category, rights scope, and simple record counts. The manifest intentionally excludes itself.

The figure and key-table input/output map is `../docs/figure_map.tsv`.

## Reproduction boundary

The saved derived tables in `processed/` are sufficient for the lightweight reproduction entry point, the Table 1 component comparison, and all 11 figure renderers. The `p6_revision_20260926_r8/` directory adds the descriptive rule comparison and its public receipt; all three files are repository-derived data under CC BY 4.0. A full raw-data rerun additionally requires downloading the 12 paired FASTQ libraries listed in `rnaseq_libraries.tsv` and installing the documented external tools. Raw reads were not redownloaded or committed to this repository.

The saved joint RNA-seq reference contains the SA1 phage and *Staphylococcus lentus* reference sequences. The TPM-like denominator is the RPK sum over the array feature plus 258 SA1 CDS features (259 features); host features are not part of that denominator.

The unfiltered NCBI `nr` BLASTP URL API request for three RT queries did **not** complete homology screening. NCBI returned a server CPU-limit failure (`SIGXCPU 24`) for all three queries with `db-num=0` and `db-len=0`. Empty hit lists in this response are not database no-hit results, and this route is excluded from evaluated screening totals. The raw terminal response and parsed interpretation are retained under `raw/p6_protein/` and `source_metadata/p6_protein/`.

## Source studies

- Yoon PH et al. *Autonomous AI agents discover reverse transcriptases with tandem repeat arrays*. Anthropic technical preprint, 23 September 2026. [Official release](https://www.anthropic.com/news/claude-discovers-novel-enzyme-system); [public PDF](https://www-cdn.anthropic.com/22573675ada52a8ca8a97a1a4b4326b2f208a071.pdf). Accessed 24 September 2026. No paper copy is bundled.
- Zhang B, Xu J, He X, Tong Y, Ren H. Interactions between Jumbo Phage SA1 and *Staphylococcus*: A Global Transcriptomic Analysis. *Microorganisms*. 2022;10(8):1590. doi:10.3390/microorganisms10081590. [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9414953/).

The 2022 publication describes the original cultures and sequencing experiment. The processing, paired contrasts, estimability, and ART-array analyses in this repository are the later reanalysis and are not attributed to that paper.

## Rights and attribution

Versioned NCBI/ENA records and their source metadata remain third-party materials. They are provided with accession-level attribution and are not relicensed by this repository. Source submitters or other parties may retain rights; this repository does not characterize all NCBI records as public domain. See the [NCBI Molecular Data Usage policy](https://www.ncbi.nlm.nih.gov/home/about/policies/) and [EMBL-EBI terms of use](https://www.ebi.ac.uk/about/terms-of-use/).

Repository-originated derived tables follow the repository's data license. Consult the top-level license and third-party notices for the exact scope. Publication PDFs and raw RNA-seq files are not redistributed here.
