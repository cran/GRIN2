
#' Get Chromosome Length
#'
#' @description
#' Retrieve chromosome size data from chr.info txt files available on the UCSC genome browser based on the user specified genome assembly.
#'
#' @param genome.assembly User can specify one of four supported genome assemblies that include "Human_GRCh38", "Human_GRCh37", "Mouse_HGCm39" and "Mouse_HGCm38".
#'
#' @details
#' Based on the genome assembly specified by the user, the function will directly retrieve chromosome size data from chr.info txt file available on the UCSC genome browser.
#'
#' @return
#' A data table with the following two columns:
#' \item{chrom}{column has the chromosome number denoted as 1, 2, X, Y, etc..}
#' \item{size}{column has the chromosome size in base pairs.}
#'
#' @export
#'
#' @importFrom circlize read.chromInfo
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author {Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [circlize::read.chromInfo()]
#'
#' @examples
#' # To retreive chromosome size data for hg19 genome assembly:
#' hg19.chrom.size=get.chrom.length("Human_GRCh37")
#' # "Human_GRCh38" can be used to retreive chromosome size data for hg38 genome assembly.

get.chrom.length=function(genome.assembly)   # function support four genome assemblies that include "Human_GRCh38", "Human_GRCh37", "Mouse_HGCm39" and "Mouse_HGCm38"
{
  # retrieve chromosome size data for GRCh38 (hg38) genome build
  if (genome.assembly=="Human_GRCh38")
  {
    chr.size.hg38= circlize::read.chromInfo(species = "hg38")
    chr.size.hg38=as.data.frame(chr.size.hg38)
    chr.size=cbind.data.frame(chrom=chr.size.hg38$chromosome,
                              size=chr.size.hg38$chr.len)
    chr.size$chrom<-gsub("chr","",as.character(chr.size$chrom))

    return(chr.size)
  }

  # retrieve chromosome size data for GRCh37 (hg19) genome build
  if (genome.assembly=="Human_GRCh37")
  {
    chr.size.hg19= circlize::read.chromInfo(species = "hg19")
    chr.size.hg19=as.data.frame(chr.size.hg19)
    chr.size=cbind.data.frame(chrom=chr.size.hg19$chromosome,
                              size=chr.size.hg19$chr.len)
    chr.size$chrom<-gsub("chr","",as.character(chr.size$chrom))

    return(chr.size)
  }

  # retrieve chromosome size data for Mouse_HGCm39 (mm39) genome build
  if (genome.assembly=="Mouse_HGCm39")
  {
    chr.size.mm39= circlize::read.chromInfo(species = "mm39")
    chr.size.mm39=as.data.frame(chr.size.mm39)
    chr.size=cbind.data.frame(chrom=chr.size.mm39$chromosome,
                              size=chr.size.mm39$chr.len)
    chr.size$chrom<-gsub("chr","",as.character(chr.size$chrom))

    return(chr.size)
  }

  # retrieve chromosome size data for Mouse_HGCm38 (mm10) genome build
  if (genome.assembly=="Mouse_HGCm38")
  {
    chr.size.mm38= circlize::read.chromInfo(species = "mm10")
    chr.size.mm38=as.data.frame(chr.size.mm38)
    chr.size=cbind.data.frame(chrom=chr.size.mm38$chromosome,
                              size=chr.size.mm38$chr.len)
    chr.size$chrom<-gsub("chr","",as.character(chr.size$chrom))

    return(chr.size)
  }
}
