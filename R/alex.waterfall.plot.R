
#' Prepare Waterfall Plot of Lesion and Expression Data
#'
#' @description
#' Function return a waterfall plot for expression data by lesion groups of a selected gene.
#'
#' @param waterfall.prep Output of the alex.waterfall.prep function. It's a list of three data tables that include "gene.lsn.exp" that has patient ID, lesion type that affect this gene if any and expression level of the selected gene, "lsns" which is a data table with all lesions affecting the gene of interest in a GRIN compatible format and "stats" which is one row with the Kruskal-Wallis test result (output of the KW.hit.express function).
#' @param lsn.data Lesion data in a GRIN compatible format.
#' @param lsn.clrs Assigned colors per lesion types. If not specified, colors will be automatically assigned using default.grin.colors function.
#' @param delta Spacing argument for the waterfall plot.
#'
#' @details
#' Function return a waterfall plot for expression data by lesion groups of a selected gene. This plot offers a side by side graphical representation of lesion and expression data for each patient where lesion groups are ordered alphabetically. For each lesion category, expression data is ordered from the lowest to the highest with patient with the median expression of the gene in the middle of the panel.
#'
#' @return Function return a waterfall plot for expression data by lesion groups of a selected gene.
#'
#' @export
#'
#' @importFrom stats median aggregate
#' @importFrom graphics rect legend text segments
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author {Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [alex.prep.lsn.expr()], [KW.hit.express()], [alex.waterfall.prep()]
#'
#' @examples
#' data(expr.data)
#' data(lesion.data)
#' data(hg19.gene.annotation)
#'
#' # prepare expression, lesion data and return the set of genes with both types of data available
#' # ordered by gene IDs in rows and patient IDs in columns:
#' alex.data=alex.prep.lsn.expr(expr.data, lesion.data,
#'                              hg19.gene.annotation, min.expr=1, min.pts.lsn=5)
#'
#' # run KW test for association between lesion groups and expression level of the same gene:
#' alex.kw.results=KW.hit.express(alex.data, hg19.gene.annotation, min.grp.size=5)
#'
#' # To prepare lesion and expression data for a waterfall plot (WT1 gene):
#' WT1.waterfall.prep=alex.waterfall.prep(alex.data, alex.kw.results, "WT1", lesion.data)
#'
#' # waterfall plot of WT1 gene:
#' WT1.waterfall.plot=alex.waterfall.plot(WT1.waterfall.prep, lesion.data)

alex.waterfall.plot=function(waterfall.prep,   # Output of the alex.waterfall.prep function
                             lsn.data,         # Lesion data in a GRIN compatible format
                             lsn.clrs=NULL,    # Colors assigned for each lesion group. If NULL, the default.grin.colors function will be used to assign lesion colors automatically
                             delta=0.5)        # spacing argument for the waterfall plot

{
  unique.grps=sort(unique(lsn.data$lsn.type))

  # assign colors for multiple and none lesion groups (colors to be assigned automatically for other lesion groups)
  if (is.null(lsn.clrs))
  {
    lsn.grps.clr=default.grin.colors(unique.grps)
    common.grps.clr=c(none="gray", multiple="violet")
    lsn.clrs=c(lsn.grps.clr, common.grps.clr)
  }

  gene.ID=waterfall.prep$gene.ID
  gene.lsn.exp=waterfall.prep$gene.lsn.exp

  lsn.clm=paste0(gene.ID,".lsn")
  rna.clm=paste0(gene.ID,".RNA")

  gene.lsn.exp.ord=order(gene.lsn.exp[,lsn.clm],
                         gene.lsn.exp[,rna.clm])

  gene.lsn.exp=gene.lsn.exp[gene.lsn.exp.ord,]

  pt.num=1:nrow(gene.lsn.exp)
  names(pt.num)=gene.lsn.exp$ID


  ####################################
  # Set up plotting region
  plot(c(-1.1,+1.5),
       c(0.1,-1.1)*nrow(gene.lsn.exp),
       type="n",axes=F,
       xlab="",ylab="")

  ###################################
  # DNA lesion plot

  # gene locus
  loc.rng=c(waterfall.prep$stats$loc.start,
            waterfall.prep$stats$loc.end)
  loc.lng=diff(loc.rng)

  pos.rng=loc.rng+c(-1,1)*delta*loc.lng
  waterfall.prep$lsns$x.start=(waterfall.prep$lsns$loc.start-pos.rng[1])/(diff(pos.rng))-1.1
  waterfall.prep$lsns$x.end=(waterfall.prep$lsns$loc.end-pos.rng[1])/(diff(pos.rng))-1.1

  x.locus=(loc.rng-pos.rng[1])/diff(pos.rng)-1.1

  # background for lesion plot
  graphics::rect(-1.1,-nrow(gene.lsn.exp),
       -0.1,0,
       col=lsn.clrs["none"],
       border=lsn.clrs["none"])

  waterfall.prep$lsns$size=waterfall.prep$lsns$loc.end-
    waterfall.prep$lsns$loc.start+1
  ord=rev(order(waterfall.prep$lsns$size))
  waterfall.prep$lsns=waterfall.prep$lsns[ord,]
  graphics::rect(pmax(waterfall.prep$lsns$x.start,-1.1),
       -pt.num[waterfall.prep$lsns$ID],
       pmin(waterfall.prep$lsns$x.end,-0.1),
       -pt.num[waterfall.prep$lsns$ID]+1,
       col=lsn.clrs[waterfall.prep$lsns$lsn.type],
       border=lsn.clrs[waterfall.prep$lsns$lsn.type])

  graphics::segments(x.locus,-nrow(gene.lsn.exp),
           x.locus,0,col="white",lty=3)

  graphics::text(x.locus,-1.05*nrow(gene.lsn.exp),
       loc.rng,cex=0.75)
  graphics::text(-0.55,+0.05*nrow(gene.lsn.exp),
       paste0(gene.ID," DNA Lesions"))

  ##############################################
  # RNA expression plot

  rna.lbls=pretty(gene.lsn.exp[,rna.clm])
  ok.lbls=(rna.lbls>min(gene.lsn.exp[,rna.clm],na.rm=T))&
    (rna.lbls<max(gene.lsn.exp[,rna.clm],na.rm=T))
  rna.lbls=rna.lbls[ok.lbls]

  gene.lsn.exp$x.rna=(gene.lsn.exp[,rna.clm]-min(gene.lsn.exp[,rna.clm],na.rm=T))/diff(range(gene.lsn.exp[,rna.clm],na.rm=T))
  gene.lsn.exp$x.rna=gene.lsn.exp$x.rna+0.1

  x.rna.lbls=(rna.lbls-min(gene.lsn.exp[,rna.clm],na.rm=T))/diff(range(gene.lsn.exp[,rna.clm],na.rm=T))+0.1

  n.mdn=function(x)
  {
    mdn=stats::median(x,na.rm=T)
    n=sum(!is.na(x))
    res=unlist(c(n=n,mdn=mdn))
    return(res)
  }

  lsn.mdn=stats::aggregate(x=gene.lsn.exp$x.rna,
                    by=list(lsn=gene.lsn.exp[,lsn.clm]),
                    FUN=n.mdn)
  rownames(lsn.mdn$x)=lsn.mdn$lsn

  x.mdn=lsn.mdn$x[gene.lsn.exp[,lsn.clm],"mdn"]

  graphics::text(x.rna.lbls,-1.05*nrow(gene.lsn.exp),
       rna.lbls,cex=0.75)
  graphics::segments(x.rna.lbls,0,
           x.rna.lbls,-nrow(gene.lsn.exp),
           lty=3)

  graphics::rect(pmin(gene.lsn.exp$x.rna,x.mdn),-(1:nrow(gene.lsn.exp)),
       pmax(gene.lsn.exp$x.rna,x.mdn),-(1:nrow(gene.lsn.exp))+1,
       col=lsn.clrs[gene.lsn.exp[,lsn.clm]],
       border=NA)

  graphics::segments(x.mdn,-(1:nrow(gene.lsn.exp)),
           x.mdn,-(1:nrow(gene.lsn.exp))+1,
           col=lsn.clrs[gene.lsn.exp[,lsn.clm]])

  graphics::text(0.55,0.05*nrow(gene.lsn.exp),
       paste0(gene.ID," RNA Expression"))

  #############################
  # Add legend

  lsn.inc=names(lsn.clrs)%in%gene.lsn.exp[,lsn.clm]
  graphics::legend(1.15,-0.25*nrow(gene.lsn.exp),
         fill=unlist(lsn.clrs[lsn.inc]),
         legend=names(lsn.clrs[lsn.inc]),
         cex=0.75,border=NA,bty="n")

}
