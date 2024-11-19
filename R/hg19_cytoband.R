
#' GRCh37 Chromosome Cytobands
#'
#' The dataset has the start and end positions in base pairs of all 22 autosomes in addition to X and Y chromosome cytobands for Human-GRCh37 (hg19) genome assembly.
#'
#' @format ## `hg19_cytoband`
#' A data frame with 862 rows and 5 columns:
#' \describe{
#'   \item{chrom}{The chromosome number.}
#'   \item{chromStart}{The cytoband start position on the chromosome in base pairs.}
#'   \item{chromEnd}{The cytoband end position on the chromosome in base pairs.}
#'   \item{name}{The cytoband name.}
#'   \item{gieStain}{The coloring scheme of the cytobands.}
#' }
#' @source The Chromosome cytobands data file was downloaded from the UCSC genome browser for GRCh37 genome assembly <https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/>.
"hg19_cytoband"
