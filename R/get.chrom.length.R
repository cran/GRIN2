
#' Get Chromosome Length
#'
#' @description
#' Retrieves chromosome size data for the human GRCh38 (hg38) genome assembly using UCSC `chromInfo` data, accessed via the `circlize` package.
#'
#' @param genome.assembly Character string specifying the genome assembly. Currently, only `"Human_GRCh38"` is supported.
#'
#' @details
#' The function fetches chromosome size information from the UCSC genome browser via the `circlize::read.chromInfo()` function. It returns data for all 22 autosomes in addition to X and Y chromosomes in the human GRCh38 (hg38) assembly. Chromosome names are formatted without the `"chr"` prefix.
#'
#' @return
#' A data frame with the following columns:
#' \describe{
#'   \item{chrom}{Chromosome identifier (e.g., 1, 2, ..., X, Y).}
#'   \item{size}{Chromosome size in base pairs.}
#' }
#'
#' @export
#'
#' @importFrom circlize read.chromInfo
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link[circlize]{read.chromInfo}}
#'
#' @examples
#' # Retrieve chromosome size data for the GRCh38 genome assembly
#' hg38.chrom.size <- get.chrom.length("Human_GRCh38")

get.chrom.length <- function(genome.assembly) {
  # Check for required package
  if (!requireNamespace("circlize", quietly = TRUE)) {
    message("Package 'circlize' is not installed. Please install it using install.packages('circlize')")
    return(NULL)
  }

  # Check supported genome assembly
  if (genome.assembly != "Human_GRCh38") {
    stop("Unsupported genome assembly. Only 'Human_GRCh38' is currently supported.")
  }

  # Retrieve chromosome sizes using circlize
  chr.size.hg38 <- circlize::read.chromInfo(species = "hg38")
  chr.size <- data.frame(
    chrom = gsub("chr", "", as.character(chr.size.hg38$chromosome)),
    size = chr.size.hg38$chr.len
  )

  return(chr.size)
}
