
#' Genomewide log10q Plot
#'
#' @description
#' The function return a genome-wide plot based on -log(10) q-value of each of the evaluated annotated genes or lesion boundaries on each chromosome. The plot is lesion type specific (gain, loss, mutation, etc...).
#'
#' @param grin.res GRIN results evaluating annotated genes or lesion boundaries (output of the grin.stats function using either lesion boundaries or annotated genes as a marker input file).
#' @param lsn.grps Selected lesion groups to be added to the plot.
#' @param lsn.colors Colors assigned to each lesion group (if NULL, lesion colors will be assigned automatically using default.grin.colors function).
#' @param max.log10q Maximum log10 q value to be added to the plot. All boundaries or genes with -log10 q smaller than this value will be set automatically to max.log10q.
#'
#' @return
#' The function return a genome-wide plot based on -log(10) q-value of each of the evaluated annotated genes or lesion boundaries to be affected by a certain type of lesions (gain, loss, mutation, etc...).
#'
#' @export
#'
#' @importFrom graphics rect legend text segments
#'
#' @author {Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [grin.lsn.boundaries()]
#'
#' @examples
#' data(lesion.data)
#' data(hg19.gene.annotation)
#' data(hg19.chrom.size)
#'
#' # This analysis is lesion type specific. So, user should first data extract data for a specific
#' # lesion group of interest for example gains from the lesion data file:
#' gain=lesion.data[lesion.data$lsn.type=="gain",]
#' # Return lesion boundaries for gains:
#' lsn.bound.gain=grin.lsn.boundaries(gain, hg19.chrom.size)
#' # Run GRIN analysis Using Lesion Boundaries as Markers Instead of the Gene Annotation File:
#' GRIN.results.gain.bound=grin.stats(gain, lsn.bound.gain, hg19.chrom.size)
#' # Return genomewide -log10q plot for association between lesion boundaries and gain:
#' genomewide.log10q.plot(GRIN.results.gain.bound, lsn.grps=c("gain"),
#'                        lsn.colors=c("gain" = "red"), max.log10q = 10)
#'
#' # instead of lesion boundaries, users can also plot -log10q values for annotated genes using
#' # genes annotation data as a marker data file:
#' grin.results=grin.stats(lesion.data,
#'                         hg19.gene.annotation,
#'                         hg19.chrom.size)
#'
#' genomewide.log10q.plot(grin.results, lsn.grps=c("gain"), lsn.colors=c("gain" = "red"),
#'                        max.log10q = 10)
#'
#' # User can run this same analysis for other lesion types such as mutations and deletions.

genomewide.log10q.plot=function(grin.res,        # GRIN results (output of the grin.stats function)
                                lsn.grps,        # selected lesion groups to be added to the plot
                                lsn.colors=NULL, # Lesion colors
                                max.log10q=NULL) # Maximum log10 q value to be added to the plot

{
  if (!is.element("x.start",colnames(grin.res$lsn.data)))
    grin.res=compute.gw.coordinates(grin.res)

  grin.res$lsn.data$x.ID=as.numeric(as.factor(grin.res$lsn.data$ID))
  n=max(grin.res$lsn.data$x.ID)
  n.chr=nrow(grin.res$chr.size)

  # set up plotting region
  plot(c(-0.05,1.2)*n,c(0,-1.1*grin.res$chr.size$x.end[n.chr]),
       type="n",axes=F,xlab="",ylab="")

  # background colors for chromosomes
  graphics::rect(0*n,-grin.res$chr.size$x.start,
       1*n,-grin.res$chr.size$x.end,
       col=c("lightgray","gray"),
       border=NA)

  lsn.types=lsn.grps
  selected.clms = list()

  for (i in 1:length(lsn.types)) {
    temp=grin.res$gene.hits[grepl(lsn.types[i],colnames(grin.res$gene.hits))]
    selected.clms[[i]] <- temp
  }

  final.data = do.call(cbind, selected.clms)

  if (is.null(lsn.colors))
  {
    lsn.colors=default.grin.colors(lsn.types)
  }
  grin.res$lsn.data$lsn.colors=lsn.colors[grin.res$lsn.data$lsn.type]

  grin.res$lsn.data$lsn.size=grin.res$lsn.data$x.end-grin.res$lsn.data$x.start
  ord=order(grin.res$lsn.data$lsn.size,decreasing=T)

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

  nsubj.data$log10q[nsubj.data$log10q>max.log10q]=max.log10q

  ord=order(nsubj.data$log10q,decreasing=T)
  nsubj.data=nsubj.data[ord,]

  graphics::segments(0*n,
           -(nsubj.data$x.start+nsubj.data$x.end)/2,
           0*n+1*nsubj.data$log10q/max(nsubj.data$log10q)*n,
           col=nsubj.data$lsn.colors)

  graphics::text(c(n,0)[(1:n.chr)%%2+1],
       pos=c(4,2)[(1:n.chr)%%2+1],
       -(grin.res$chr.size$x.start+grin.res$chr.size$x.end)/2,
       grin.res$chr.size$chrom,
       cex=0.75)

  graphics::legend(n/2,-1.05*grin.res$chr.size$x.end[n.chr],
         fill=lsn.colors,cex=0.9,
         legend=names(lsn.colors),
         xjust=0.5,
         ncol=length(lsn.colors),
         border=NA,bty="n")

  graphics::text(1*n/2,0,
       "-log10(q)",cex=0.95,
       pos=3)

  graphics::text(c(0,0.25, 0.5, 0.75, 1)*n,
       -grin.res$chr.size$x.end[n.chr],
       c(0,max.log10q/4, max.log10q/2, round(max.log10q/1.3333333333, 1), max.log10q),
       cex=0.75,pos=1)
}

