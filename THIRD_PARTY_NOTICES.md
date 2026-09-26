# Source data, attribution and license scope

The repository's [MIT License](LICENSE) applies to original analysis/plotting
software and accompanying documentation. Original analytical tables, summary
matrices and figures use [CC BY 4.0](LICENSE-DATA.md). Neither grant changes the
rights in source records or external dependencies.

| Material | Repository location / source | Terms and attribution |
|---|---|---|
| Public genome and protein records, source annotations, and sequence extracts | `data/raw/`, FASTA/GenBank-derived files under `data/processed/`; individual accessions and source links in [data/README.md](data/README.md) and its manifests | NCBI's [Molecular Data Usage policy](https://www.ncbi.nlm.nih.gov/home/about/policies/) allows use and redistribution of molecular data while noting that submitters may retain rights. Inclusion here is not a claim that all source records are public domain. Retain accession and original study attribution. |
| RNA-seq reads and archive metadata | [PRJNA836150](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA836150), ENA/SRA run accessions and download manifest under `data/`; raw FASTQ reads are retrieved from their public archive rather than bundled | ENA is subject to [EMBL-EBI terms of use](https://www.ebi.ac.uk/about/terms-of-use/). EMBL-EBI does not add restrictions beyond source rights; contributed data may retain third-party rights and should be attributed to their source. The source archive records govern use. |
| Derived counts, comparisons and summary matrices | Analytical TSV/JSON outputs under `data/processed/` | [CC BY 4.0](LICENSE-DATA.md) covers this project's original analysis contributions. Source sequence, annotation and metadata fields retain their source rights and accessions. |
| Research figures and workflow | `figures/reference/`, `docs/analysis_workflow.svg`, `docs/analysis_workflow.png` | Original plotted results and workflow artwork: [CC BY 4.0](LICENSE-DATA.md). No third-party article pages, logos, gallery images or manuscript files are included. |
| Python and command-line dependencies | Declared requirements and the software named in [reproduction documentation](docs/reproducibility.md) | Dependencies are installed separately, not redistributed as vendored packages or executable binaries. Their own licenses apply. Cite the relevant scientific methods and software as described in the documentation. |

Source accessions, acquisition dates where recorded, and download URLs are
provided with the data manifests. The research citations in the README and
reproduction documentation identify prior discoveries and experimental data;
this repository does not claim ownership of those discoveries or their datasets.
