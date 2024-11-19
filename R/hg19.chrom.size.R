
#' Chromosome Length Data File
#'
#' The file has the size of 22 autosomes in addition to X and Y chromosomes in base pairs directly retrieved from chr.info txt files available on the UCSC genome browser using get.chrom.length function and "Human-GRCh37" as a genome assembly option (hg19).
#'
#' @format ## `hg19.chrom.size`
#' A data frame with 24 rows and 2 columns:
#' \describe{
#'   \item{chrom}{The chromosome number.}
#'   \item{size}{The chromosome length in base pairs.}
#' }
#' @source Chromosome size data directly retrieved from chr.info txt files available on the UCSC genome browser using get.chrom.length function and "Human-GRCh37" as a genome assembly option (hg19).
"hg19.chrom.size"
