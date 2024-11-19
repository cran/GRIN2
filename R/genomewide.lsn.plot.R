
#' Genome-wide Lesion Plot
#'
#' @description
#' Function return a genomewide lesion plot for all lesion types affecting different chromosomes.
#'
#' @param grin.res GRIN results (output of the grin.stats function).
#' @param ordered By default the function will order the patient IDs alphabetically. However, users can specify a certain patient's order in the genomewide lesion plot by specifying ordered=TRUE and pass a data frame with new patient's order to the pt.order argument.
#' @param pt.order data.frame of two columns "ID" that has patient IDs matching the unique IDs in the lesion data file and "pts.order" that has the new patient's order listed as numbers that range from 1:n.patients (Should be only specified if ordered=TRUE).
#' @param lsn.colors a vector of lesion colors (If not provided by the user, colors will be automatically assigned using default.grin.colors function).
#' @param max.log10q Maximum log10 q value for genes in the GRIN results table to be added to the plot. If max.log10q=100 for example, all -log10q values>100, will be adjusted to 100 in the plot.
#'
#' @details
#' The function use the genome-wide plotting coordinates obtained from the compute.gw.coordinates function and plot the whole set of lesions affecting subjects included in the dataset in the middle panel of the figure. Two additional side panels show the number of affected subjects and -log10 q value of each locus to be affected by all different types of lesions.
#'
#' @return
#' The function return a genome-wide lesion plot (all chromosomes) in the middle panel. For each locus, Panel on the left shows -log10 q value and the Panel on the right show the number of subjects affected by all different types of lesions color coded by lesion category.
#'
#' @export
#'
#' @importFrom graphics rect legend text segments
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author {Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [compute.gw.coordinates()]
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
#'
#' # prepare the genomewide lesion plot using genomewide.lsn.plot function with patient IDs ordered
#' # alphabetically:
#' genomewide.plot=genomewide.lsn.plot(grin.results, max.log10q=50)
#'
#' # To pass certain patients order to the genomewide.lsn.plot function, the user should specify
#' # a certain patients order using the pt.order argument.

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
    lesions=grin.res$lsn.data
    merged.data=merge(ordered.pts,lesions,by="ID", all.y=TRUE)
    grin.res$lsn.data=merged.data
    n=max(grin.res$lsn.data$pts.order)
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
