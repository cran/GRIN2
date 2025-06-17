
#' Chromosome Length Data (hg38)
#'
#' This dataset contains the lengths (in base pairs) of the 22 autosomes in addition to the X and Y chromosomes based on the GRCh38 human genome assembly. The data was retrieved from the UCSC Genome Browser using the `get.chrom.length` function with "Human-GRCh38" as the selected genome.
#'
#' @format ## `hg38_chrom_size`
#' A data frame with 24 rows and 2 columns:
#' \describe{
#'   \item{chrom}{Chromosome name (1,2, 3, ..., X, Y).}
#'   \item{size}{Chromosome length in base pairs.}
#' }
#'
#' @source Retrieved from UCSC Genome Browser `chr.info` text files using the `get.chrom.length` function with the GRCh38 genome assembly.
"hg38_chrom_size"
