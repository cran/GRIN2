
#' Genome-wide Lesion Plot
#'
#' @description
#' Generates a genome-wide lesion plot displaying all lesion types affecting different chromosomes.
#'
#' @param grin.res GRIN results (output from the `grin.stats` function).
#' @param ordered Logical; if `TRUE`, patient IDs will be reordered according to the `pt.order` data frame. If `FALSE` (default), patient IDs are ordered alphabetically.
#' @param pt.order A data frame with two columns: `"ID"` (patient identifiers matching those in the lesion data) and `"pts.order"` (numeric vector specifying the new patient order from 1 to n patients). Only required if `ordered = TRUE`.
#' @param lsn.colors A named vector of colors assigned to lesion types. If not provided, colors will be automatically assigned using the `default.grin.colors` function.
#' @param max.log10q Numeric; maximum value for -log10(q-value) used in the plot. Any value greater than `max.log10q` will be capped at this value in the left panel of the plot.
#'
#' @details
#' This function uses genome-wide coordinates (from `compute.gw.coordinates`) to generate a three-panel plot. The middle panel shows lesions by chromosome across patients. The left panel displays the -log10(q-values) from the GRIN results for each gene, and the right panel shows the number of patients affected at each locus, color-coded by lesion type.
#'
#' @return
#' A genome-wide lesion plot consisting of three aligned panels:
#' \itemize{
#'   \item Middle panel: genome-wide lesion map across all chromosomes and patients.
#'   \item Left panel: -log10(q-values) of each locus from GRIN results showing Statistical Significance of Lesion Frequencies.
#'   \item Right panel: number of affected patients at each locus, colored by lesion category.
#' }
#'
#' @export
#'
#' @importFrom graphics rect legend text segments
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link{compute.gw.coordinates}}
#'
#' @examples
#' data(lesion_data)
#' data(hg38_gene_annotation)
#' data(hg38_chrom_size)
#'
#' # Run GRIN analysis
#' grin.results <- grin.stats(lesion_data,
#'                            hg38_gene_annotation,
#'                            hg38_chrom_size)
#'
#' # Generate genome-wide lesion plot with alphabetical patient ordering
#' genomewide.plot <- genomewide.lsn.plot(grin.results, max.log10q = 50)

genomewide.lsn.plot=function(grin.res,        # GRIN results (output of the grin.stats function)
                             ordered=FALSE,   # user can specify a certain patients order in the genowide lesion plot by specifying ordered = TRUE and pass pts.order data.frame with the new patients order , otherwise patients will be ordered alphabetically
                             pt.order=NULL,   # data.frame with two columns "ID" that has the patient ID and "pts.order" tha has the patient order in numbers
                             lsn.colors=NULL, # Lesion colors (If not provided by the user, colors will be automatically assigned using default.grin.colors function).
                             max.log10q=NULL) # Maximum log10 q value for genes in the GRIN results table to be added to the plot. If max.log10q=100, all -log10q values<100, will be adjusted to 100 in the plot

{
  if (!is.element("x.start",colnames(grin.res$lsn.data)))
    grin.res=compute.gw.coordinates(grin.res)

  if (ordered == FALSE) {
    grin.res$lsn.data$pts.order=as.numeric(as.factor(grin.res$lsn.data$ID))
    n=max(grin.res$lsn.data$pts.order)
    n.chr=nrow(grin.res$chr.size)
  } else {
    ordered.pts=pt.order
    n=max(ordered.pts$pts.order)
    lesions=grin.res$lsn.data
    merged.data=merge(ordered.pts,lesions,by="ID", all.y=TRUE)
    grin.res$lsn.data=merged.data

    n.chr=nrow(grin.res$chr.size)
  }

  # set up plotting region
  plot(c(-0.2,1.2)*n,c(0,-1.1*grin.res$chr.size$x.end[n.chr]),
       type="n",axes=F,xlab="",ylab="")

  # background colors for chromosomes
  graphics::rect(0,-grin.res$chr.size$x.start,
                 n,-grin.res$chr.size$x.end,
                 col=c("lightgray","gray"),
                 border=NA)

  graphics::rect(-0.075*n,-grin.res$chr.size$x.start,
                 -0.2*n,-grin.res$chr.size$x.end,
                 col=c("lightgray","gray"),
                 border=NA)

  graphics::rect(1.075*n,-grin.res$chr.size$x.start,
                 1.2*n,-grin.res$chr.size$x.end,
                 col=c("lightgray","gray"),
                 border=NA)

  lsn.types=sort(unique(grin.res$lsn.index$lsn.type))

  if (is.null(lsn.colors))
  {
    lsn.colors=default.grin.colors(lsn.types)
  }
  grin.res$lsn.data$lsn.colors=lsn.colors[grin.res$lsn.data$lsn.type]

  grin.res$lsn.data$lsn.size=grin.res$lsn.data$x.end-grin.res$lsn.data$x.start
  ord=order(grin.res$lsn.data$lsn.size,decreasing=T)

  graphics::rect(grin.res$lsn.data$pts.order-1,
                 -grin.res$lsn.data$x.start,
                 grin.res$lsn.data$pts.order,
                 -grin.res$lsn.data$x.end,
                 col=grin.res$lsn.data$lsn.colors,
                 border=grin.res$lsn.data$lsn.colors)

  graphics::text(c(n,0)[(1:n.chr)%%2+1],
                 pos=c(4,2)[(1:n.chr)%%2+1],
                 -(grin.res$chr.size$x.start+grin.res$chr.size$x.end)/2,
                 grin.res$chr.size$chrom,
                 cex=0.5)

  graphics::legend(n/2,-1.025*grin.res$chr.size$x.end[n.chr],
                   fill=lsn.colors,cex=0.75,
                   legend=names(lsn.colors),
                   xjust=0.5,
                   ncol=length(lsn.colors),
                   border=NA,bty="n")

  nsubj.mtx=unlist(grin.res$gene.hits[,paste0("nsubj.",lsn.types)])
  qval.mtx=unlist(grin.res$gene.hits[,paste0("q.nsubj.",lsn.types)])
  nsubj.data=cbind.data.frame(gene=grin.res$gene.hits$gene,
                              x.start=grin.res$gene.hits$x.start,
                              x.end=grin.res$gene.hits$x.end,
                              nsubj=nsubj.mtx,
                              log10q=-log10(qval.mtx),
                              lsn.type=rep(lsn.types,each=nrow(grin.res$gene.hits)))
  nsubj.data=nsubj.data[nsubj.data$nsubj>0,]
  nsubj.data$lsn.colors=lsn.colors[nsubj.data$lsn.type]

  ord=order(nsubj.data$nsubj,decreasing=T)
  nsubj.data=nsubj.data[ord,]

  graphics::segments(1.075*n,
                     -(nsubj.data$x.start+nsubj.data$x.end)/2,
                     1.075*n+0.125*nsubj.data$nsubj/max(nsubj.data$nsubj)*n,
                     col=nsubj.data$lsn.colors)

  nsubj.data$log10q[nsubj.data$log10q>max.log10q]=max.log10q

  ord=order(nsubj.data$log10q,decreasing=T)
  nsubj.data=nsubj.data[ord,]

  graphics::segments(-0.075*n,
                     -(nsubj.data$x.start+nsubj.data$x.end)/2,
                     -0.075*n-0.125*nsubj.data$log10q/max(nsubj.data$log10q)*n,
                     col=nsubj.data$lsn.colors)

  graphics::text(-(0.075+0.20)*n/2,0,
                 "-log10(q)",cex=0.75,
                 pos=3)

  graphics::text((1.075+1.2)*n/2,0,
                 "Subjects",cex=0.75,
                 pos=3)

  graphics::text(c(-0.075,1.075)*n,
                 -grin.res$chr.size$x.end[n.chr],
                 0,cex=0.75,pos=1)

  graphics::text(c(-0.2,1.2)*n,
                 -grin.res$chr.size$x.end[n.chr],
                 c(max(nsubj.data$log10q),
                   max(nsubj.data$nsubj)),
                 cex=0.75,pos=1)

}
