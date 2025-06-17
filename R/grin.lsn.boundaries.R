
#' GRIN Evaluation of Lesion Boundaries
#'
#' @description
#' This function evaluates copy number variations (CNVs), specifically gains and deletions, using unique lesion start and end positions to define genomic boundaries. The analysis is lesion-type specific and spans the entire genome.
#'
#' @param lsn.data A lesion data table containing only gain or deletion events. If gains are subdivided (e.g., into gains and amplifications based on log2Ratio values), both subtypes can be included. The same applies for deletions (e.g., homozygous and heterozygous).
#' @param chrom.size A chromosome size table with two required columns: `"chrom"` (chromosome identifier) and `"size"` (chromosome size in base pairs).
#'
#' @details
#' The function identifies unique CNV boundaries by evaluating the start and end positions of lesions on each chromosome. Large lesions may be split into smaller boundaries based on overlapping smaller lesions in other samples. This boundary-based approach ensures comprehensive genome-wide coverage, including intergenic regions and areas without annotated features.
#'
#' The first boundary on each chromosome spans from the start of the chromosome to the start of the first lesion affecting any of the included patient samples. Similarly, the last boundary extends from the end of the last lesion to the end of the chromosome.
#'
#' This method is particularly useful when analyzing CNV data at a finer resolution than gene-level annotation allows.
#'
#' @return
#' A `data.frame` with five columns:
#' \item{gene}{Boundary identifier, based on unique start and end positions.}
#' \item{chrom}{Chromosome on which the boundary resides.}
#' \item{loc.start}{Start position of the boundary.}
#' \item{loc.end}{End position of the boundary.}
#' \item{diff}{Length of the boundary in base pairs.}
#'
#' @export
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link{grin.stats}}
#'
#' @examples
#' data(lesion_data)
#' data(hg38_chrom_size)
#'
#' # This analysis is lesion-type specific. For example, extract gains:
#' gain <- lesion_data[lesion_data$lsn.type == "gain", ]
#'
#' # Generate lesion boundaries for gains:
#' lsn.bound.gain <- grin.lsn.boundaries(gain, hg38_chrom_size)
#'
#' # Run GRIN using lesion boundaries as markers instead of gene annotations:
#' GRIN.results.gain.bound <- grin.stats(gain, lsn.bound.gain, hg38_chrom_size)
#'
#' # The same analysis can be applied to deletions, mutations, or structural rearrangements.

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
