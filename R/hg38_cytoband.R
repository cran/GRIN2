
#' GRCh38 Chromosome Cytobands
#'
#' The dataset has the start and end positions in base pairs of all 22 autosomes in addition to X and Y chromosome cytobands for Human-GRCh38 (hg38) genome assembly.
#'
#' @format ## `hg38_cytoband`
#' A data frame with 1,549 rows and 5 columns:
#' \describe{
#'   \item{chrom}{The chromosome number.}
#'   \item{chromStart}{The cytoband start position on the chromosome in base pairs.}
#'   \item{chromEnd}{The cytoband end position on the chromosome in base pairs.}
#'   \item{name}{The cytoband name.}
#'   \item{gieStain}{The coloring scheme of the cytobands.}
#' }
#' @source The Chromosome cytobands data file was downloaded from the UCSC genome browser for GRCh38 genome assembly <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/>.
"hg38_cytoband"
