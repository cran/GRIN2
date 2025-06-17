
#' Compute Genome-wide Plotting Coordinates
#'
#' @description
#' Computes and assigns genome-wide plotting coordinates to lesion, gene, and chromosome data for use in genome-wide lesion plots.
#'
#' @param grin.res GRIN results, typically the output of the `grin.stats` function.
#' @param scl Chromosome unit length in base pairs. Default is 1,000,000, meaning each chromosome is divided into segments of 1 million base pairs for plotting.
#'
#' @details
#' This function processes the GRIN results to add genome-wide x-axis coordinates necessary for plotting lesions and genes across all chromosomes. It divides each chromosome into segments based on the specified `scl` value and computes cumulative start and end positions across chromosomes to ensure a continuous x-axis. Specifically:
#' \itemize{
#'   \item Chromosome sizes are updated to include `x.start` and `x.end` columns, where each chromosome starts where the previous one ends.
#'   \item Gene and lesion data are similarly updated with `x.start` and `x.end` coordinates, scaled by `scl`, and adjusted for cumulative chromosome positions.
#' }
#'
#' @return
#' A list identical in structure to the original `grin.res` object, with the following additions:
#' \describe{
#'   \item{gene.hits}{Unchanged. GRIN gene-level summary statistics, including hit counts and p/q-values.}
#'   \item{gene.lsn.data}{Unchanged. Gene-lesion overlaps showing which lesion affects which gene for each patient.}
#'   \item{lsn.data}{Input lesion data with added `x.start` and `x.end` columns for genome-wide coordinates.}
#'   \item{gene.data}{Input gene annotation data with added `x.start` and `x.end` columns for genome-wide coordinates.}
#'   \item{chr.size}{Chromosome size table (22 autosomes + X and Y) with added `x.start` and `x.end` columns for plotting.}
#'   \item{gene.index}{Mapping of `gene.lsn.data` rows to their corresponding chromosomes.}
#'   \item{lsn.index}{Mapping of `gene.lsn.data` rows to their corresponding lesions.}
#' }
#'
#' @export
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link{grin.stats}}
#'
#' @examples
#' data(lesion_data)
#' data(hg38_gene_annotation)
#' data(hg38_chrom_size)
#'
#' # Run GRIN model
#' grin.results <- grin.stats(lesion_data,
#'                            hg38_gene_annotation,
#'                            hg38_chrom_size)
#'
#' # Assign genome-wide coordinates for plotting
#' genome.coord <- compute.gw.coordinates(grin.results)

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
