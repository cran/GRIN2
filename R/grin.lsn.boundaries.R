
#' GRIN Evaluate Lesion Boundaries
#'
#' @description
#' The function evaluates Copy number variations that include gain and deletions as boundaries based on unique lesion start and end positions. This analysis is lesion type specific and covers the entire genome.
#'
#' @param lsn.data Lesion data file that should be limited to include either gain or deletions. If gains are splitted to gain and amplifications based on the log2Ratio value of the CNV segmentation file, the two categories can be included in the same data table, same for homozygous and heterozygous deletions.
#' @param chrom.size Chromosize size table that should include two columns "chrom" with the chromosome number and "size" with the chromosome size in base pairs.
#'
#' @details
#' The function evaluates Copy number variations that include gain and deletions as boundaries and return a table of ordered boundaries based on the unique start and end positions of different lesions on each chromosome. If gains are splitted to gain and amplifications based on the log2Ratio value of the CNV segmentation file, the two categories can be included in the same analysis, same for homozygous and heterozygous deletions. Boundary will be the region between each unique start and end positions where large size lesions will be splitted into multiple boundaries based on other smaller size lesions that affect the same region in other patients if any. This analysis is meant to cover the entire genome, so regions without any annotated genes or regulatory features will be incuded will be assesed in the analysis. The first boundary for each chromosome will start from the first nucleotide base on the chromosome till the start position of the first lesion that affect the chromosome. Similarly, the last boundary will start from the end position of the last lesion that affect the chromosome till the last base on the chromosome.
#'
#' @return
#' Function return a data.frame with five columns:
#' \item{gene}{Ordered boundaries by unique start and end positions of different lesions on each chromosome.}
#' \item{chrom}{Chromosome on which the bounday is located.}
#' \item{loc.start}{Boundary start position.}
#' \item{loc.end}{Boundary end position.}
#' \item{diff}{Boundary size in base pairs.}
#'
#' @export
#'
#' @author {Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @examples
#' data(lesion.data)
#' data(hg19.chrom.size)
#'
#' # This analysis is lesion type specific. So, user should first data extract data for a specific
#' # lesion group of interest for example gains from the lesion data file:
#' gain=lesion.data[lesion.data$lsn.type=="gain",]
#' # Return lesion boundaries for gains:
#' lsn.bound.gain=grin.lsn.boundaries(gain, hg19.chrom.size)
#' # Run GRIN analysis Using Lesion Boundaries markers Instead of the gene annotation file:
#' GRIN.results.gain.bound=grin.stats(gain, lsn.bound.gain, hg19.chrom.size)
#'
#' # same analysis can be done for mutations, deletions and structural rearrangments.

grin.lsn.boundaries=function(lsn.data, # Lesion data file that should be limited to include gain OR deletions (If gains are splitted to gain and amplification based on the log2Ratio value of the CNV segmentation file, the two categories can be included, same for homozygous and heterozygous deletions)
                             chrom.size) # Chromosize size table
{
  lsn.data=lsn.data
  chr.size=chrom.size
  chroms=chr.size$chrom
  lsn.bound = data.frame()
  {
    for (i in chroms) {
      chr=i
      chr.lsns=lsn.data[lsn.data$chrom==chr,]
      ord=order(chr.lsns$loc.start,
                chr.lsns$loc.end)

      # specify all unique start and end positions for included CNVs
      uniq.start=unique(chr.lsns$loc.start)
      uniq.end=unique(chr.lsns$loc.end)
      chrom.size=chr.size$size[chr.size$chrom==chr]

      # boundary will be the region between each unique start and end positions
      # large size lesions will be splitted into multiple boundaries based on other smaller size lesions that affect the same region in other patients if any
      all.start=unique(c(1,uniq.start,uniq.end+1))
      all.end=unique(c(uniq.start-1,uniq.end,chrom.size))
      all.end=all.end[!all.end ==0]
      all.start=sort(all.start)
      all.end=sort(all.end)
      chr.lsn.loci=cbind.data.frame(gene=paste0("chr",chr,"_",all.start,"_",all.end),
                                    chrom=chr,
                                    loc.start=all.start,
                                    loc.end=all.end)
      chr.lsn.bound = chr.lsn.loci
      lsn.bound=rbind(lsn.bound, chr.lsn.loci)

    }
    lsn.bound$diff=(lsn.bound$loc.end - lsn.bound$loc.start)
    lsn.bound=lsn.bound[!lsn.bound$diff<0,]

    return(lsn.bound)

  }
}
