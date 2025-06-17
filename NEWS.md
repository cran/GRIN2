
# GRIN2 2.0.0

## Major updates
- Significant improvements to the `prob.hits()` function, greatly enhancing the performance and speed of probability convolution calculations for measuring statistical significance of lesion frequencies.
- Removed the `row.bern.conv()` function and introduced `pbc()` and `rpbc()`, which compute the probability that a series of independent Bernoulli trials yields *x* or more successes.
- Improved compatibility with the Human GRCh38 (hg38) genome assembly.
- Some packages that handle Human GRCh37 (hg19) assembly annotations are either outdated or have platform compatibility issues. Therefore, the support for hg19 has been discontinued in `get.ensembl.annotation()` and `lsn.transcripts.plot()`. These functions now exclusively support hg38, and users are encouraged to convert their lesion data coordinates to hg38 before running GRIN2 analyses.

## New features
- Enhanced error handling with informative messages for missing annotations or data inputs.
- Added a new vignette, `GRIN2`, which demonstrates the package’s preprocessing, analysis, and plotting capabilities.

## Data and Annotation
- Included bundled datasets for GRCh38: `lesion_data`, `expr_data`, `hg38_gene_annotation`, `hg38_chrom_size`, and `hg38_cytoband`.
- Ensembl and regulatory annotations are retrieved directly from Ensembl BioMart v110 with graceful fallback mechanisms.

## Bug fixes and improvements
- Improved Roxygen documentation for better CRAN compliance.
- `get.ensembl.annotation()` and `get.chrom.length()` now handle database connection issues with informative error messages.
- Fixed minor bugs in `genomewide.lsn.plot()`, specifically regarding color assignment for lesion groups when not automatically specified by `default.grin.colors()`.
