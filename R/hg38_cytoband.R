
#' GRCh38 Chromosome Cytobands
#'
#' This dataset contains the start and end positions (in base pairs) of cytogenetic bands (cytobands) for all 22 autosomes, as well as the X and Y chromosomes, based on the Human GRCh38 (hg38) genome assembly.
#'
#' @format ## `hg38_cytoband`
#' A data frame with 1,549 rows and 5 columns:
#' \describe{
#'   \item{chrom}{Chromosome name (1:22, X, Y).}
#'   \item{chromStart}{Start position of the cytoband in base pairs.}
#'   \item{chromEnd}{End position of the cytoband in base pairs.}
#'   \item{name}{Name of the cytoband (e.g., p11.1, q22.3).}
#'   \item{gieStain}{Staining pattern (e.g., gpos, gneg, acen) used for cytogenetic visualization.}
#' }
#'
#' @source Retrieved from the UCSC Genome Browser (GRCh38 assembly):
#' <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/>
"hg38_cytoband"
