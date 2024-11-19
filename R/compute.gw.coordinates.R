
#' Compute Genome-wide Coordinates
#'
#' @description
#' The function assign plotting coordinates necessary for the genome-wide lesion plot.
#'
#' @param grin.res GRIN results (output of the grin.stats function).
#' @param scl length of chromosome units in base pairs. Default is 1,000,000 which means that each chromosome will be divided into multiple pieces each is 1 million base pair in length.
#'
#' @details
#' The function divides each chromosome into multiple units based on the specified scl value. In addition, it orders and adds two columns x.start and x.end to the chromosme size file (x.start for chr2 is equal to x.end of chr1). Function also adds x.start and x.end columns to lesion and gene annotation data files (x.start is the start position of the lesion or the gene divided by scl and x.end is the end position of the lesion or the gene divided by scl taking into consideration that the start position of the chromosomes is added consecutively based on the chromosomes length).
#'
#' @return
#' Function return a list of GRIN results with the following changes to allow adding genome-wide plotting coordinates:
#' \item{gene.hits}{No changes, a data table of GRIN results that includes gene annotation, number of subjects and number of hits affecting each locus, p and FDR adjusted q-values showing the probability of each locus to be affected by one or a constellation of multiple types of lesions.}
#' \item{gene.lsn.data}{No changes, each row represent a gene overlapped by a certain lesion. Column "gene" shows the overlapped gene ensembl ID, and ID column has the patient ID}
#' \item{lsn.data}{input lesion data with two additional columns (x.start and x.end). x.start is the start position of the lesion divided by scl and x.end is the end position of the lesion divided by scl taking into consideration that the start position of the chromosomes is added consecutively based on the chromosomes length.}
#' \item{gene.data}{input gene annotation data with two additional columns (x.start and x.end). x.start is the start position of the gene divided by scl and x.end is the end position of the gene divided by scl taking into consideration that the start position of the chromosomes is added consecutively based on the chromosomes length.}
#' \item{chr.size}{data table showing the size of the 22 autosomes, in addition to X and Y chromosomes in base pairs with two additional columns (x.start and x.end). x.start is the start position of the chromosome divided by scl and x.end is the end position of the chromosome divided by scl taking into consideration that the start position of the chromosomes is added consecutively based on the chromosomes length.}
#' \item{gene.index}{data.frame with overlapped gene-lesion data rows that belong to each chromosome in the gene.lsn.data table.}
#' \item{lsn.index}{data.frame that shows the overlapped gene-lesion data rows taht belong to each lesion in the gene.lsn.data table.}
#'
#' @export
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author {Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [grin.stats()]
#'
#' @examples
#' data(lesion.data)
#' data(hg19.gene.annotation)
#' data(hg19.chrom.size)
#'
#' # Run GRIN model using grin.stats function
#' grin.results=grin.stats(lesion.data,
#'                         hg19.gene.annotation,
#'                         hg19.chrom.size)
#' # assign genomewide coordinates and prepare the results for the genomewide.lsn.plot function
#' genome.coord=compute.gw.coordinates(grin.results)

compute.gw.coordinates=function(grin.res,    # GRIN results (output of the grin.stats function)
                                scl=1000000) # length of chromosome units in base pairs

{
  # Compute new coordinates for chromosomes
  cum.size=cumsum(grin.res$chr.size$size/scl)
  n.chr=nrow(grin.res$chr.size)
  grin.res$chr.size$x.start=c(0,cum.size[-n.chr])
  grin.res$chr.size$x.end=cum.size

  grin.res$gene.data$x.start=NA
  grin.res$gene.data$x.end=NA
  grin.res$gene.hits$x.start=NA
  grin.res$gene.hits$x.end=NA
  grin.res$lsn.data$x.start=NA
  grin.res$lsn.data$x.end=NA

  gene.data.ord=order(grin.res$gene.data$gene.row)
  grin.res$gene.data=grin.res$gene.data[gene.data.ord,]

  gene.hits.ord=order(grin.res$gene.hits$gene.row)
  grin.res$gene.hits=grin.res$gene.hits[gene.hits.ord,]

  lsn.ord=order(grin.res$lsn.data$lsn.row)
  grin.res$lsn.data=grin.res$lsn.data[lsn.ord,]


  ngi=nrow(grin.res$gene.index)
  for (i in 1:ngi)
  {
    chr.mtch=which(grin.res$gene.index$chrom[i]==grin.res$chr.size$chrom)
    chr.start=grin.res$chr.size$x.start[chr.mtch]

    gene.rows=(grin.res$gene.index$row.start[i]:grin.res$gene.index$row.end[i])
    grin.res$gene.data$x.start[gene.rows]=chr.start+grin.res$gene.data$loc.start[gene.rows]/scl
    grin.res$gene.data$x.end[gene.rows]=chr.start+grin.res$gene.data$loc.end[gene.rows]/scl

    grin.res$gene.hits$x.start[gene.rows]=chr.start+grin.res$gene.hits$loc.start[gene.rows]/scl
    grin.res$gene.hits$x.end[gene.rows]=chr.start+grin.res$gene.hits$loc.end[gene.rows]/scl
  }

  nli=nrow(grin.res$lsn.index)
  for (i in 1:nli)
  {
    chr.mtch=which(grin.res$lsn.index$chrom[i]==grin.res$chr.size$chrom)
    chr.start=grin.res$chr.size$x.start[chr.mtch]

    lsn.rows=(grin.res$lsn.index$row.start[i]:grin.res$lsn.index$row.end[i])
    grin.res$lsn.data$x.start[lsn.rows]=chr.start+grin.res$lsn.data$loc.start[lsn.rows]/scl
    grin.res$lsn.data$x.end[lsn.rows]=chr.start+grin.res$lsn.data$loc.end[lsn.rows]/scl
  }

  return(grin.res) # modified GRIN results to allow adding genome-wide coordinates
}
