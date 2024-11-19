
#' Associate Lesions with Expression Data on the Pathway Level
#'
#' @description
#' Function compute the distance between subjects in the dataset based on the lesions that affect different genes assigned to the pathway of interest and return two panels of lesion and expression data of ordered subjects based on the computed distances.
#'
#' @param alex.data output of the alex.prep.lsn.expr function. It's a list of three data tables that include "row.mtch", "alex.expr" with expression data, "alex.lsn" with lesion data. Rows of alex.expr, and "alex.lsn" matrices are ordered by gene ensembl IDs and columns are ordered by patient ID.
#' @param lsn.data Lesion data in a GRIN compatible format. data.frame should has five columns that include "ID" with patient ID, "chrom" which is the chromosome on which the lesion is located, "loc.start" which is the lesion start position, "loc.end" the lesion end position and "lsn.type" which is the lesion type for example gain, loss, mutation, fusion, etc...
#' @param pathways data.frame with three columns "gene.name" that has gene symbols, "ensembl.id" with gene ensembl ID and "pathway" that has the pathway name.
#' @param selected.pathway The pathway of interest.
#'
#' @details
#' Function compute the distance between subjects in th dataset based on lesions affecting different genes assigned to the pathway of interest and return two panels of lesion and expression data of ordered subjects based on the computed distances. Function also return a data.frame with lesion and expression data of the pathway genes ordered based on the hierarchical clustering analysis (same order of the subjects in the lesion and expression panels of the figure).
#'
#' @return
#' Function will return two panels figure of lesion and expression data of ordered subjects based on the computed distances of lesions in all genes assigned to the pathway of interest. The function will also return:
#' \item{ordered.path.data}{data.frame with lesion and expression data of the pathway genes ordered based on the hiearchial clustering analysis (same order of the subjects in the lesion and expression panels of the figure).}
#'
#' @export
#'
#' @importFrom tibble rownames_to_column
#' @importFrom stats hclust as.dist sd
#' @importFrom grDevices rgb
#' @importFrom graphics par rect legend text
#'
#' @author {Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [alex.prep.lsn.expr()], [stats::hclust()]
#'
#' @examples
#' data(expr.data)
#' data(lesion.data)
#' data(hg19.gene.annotation)
#' data(pathways)
#'
#' # prepare expression, lesion data and return the set of genes with both types of data available
#' # ordered by gene IDs in rows and patient IDs in columns:
#' alex.data=alex.prep.lsn.expr(expr.data, lesion.data,
#'                              hg19.gene.annotation, min.expr=5,
#'                              min.pts.lsn=5)
#'
#' # use lesions in all genes assigned to the jak_pathway as an example pathway:
#' alex.path=alex.pathway(alex.data, lesion.data, pathways, "Jak_Pathway")
#' # extract expression and lesion data (same subjects order in the figure)
#' alex.path

alex.pathway=function(alex.data,          # output of the alex.prep.lsn.expr function (list of three data table "alex.expr" with expression data ready for KW test, "alex.lsn" with overlapped gene lesion data and row.mtch)
                      lsn.data,           # lesion data in a GRIN compatible format
                      pathways,           # data.table with three columns "gene.name" with gene symbols, "ensembl.id" with gene ensembl ID and "pathway" that has the pathway name
                      selected.pathway)   # pathway of interest

{
  expr=alex.data$alex.expr
  lsn=alex.data$alex.lsn
  selected.genes=pathways[pathways$pathway==selected.pathway,]
  pathway.ensembl=selected.genes$ensembl.id

  # extract lesion and expression data for genes assigned to the pathway of interest
  path.expr=expr[(rownames(expr)%in%pathway.ensembl),]
  path.lsn=lsn[(rownames(lsn)%in%pathway.ensembl),]
  path.expr=as.matrix(path.expr)
  path.lsn=as.matrix(path.lsn)

  # assign colors to lesion groups
  unique.grps=sort(unique(lsn.data$lsn.type))
  lsn.grps.clr=default.grin.colors(unique.grps)
  common.grps.clr=c(none="gray", multiple="violet")
  lsn.clrs=c(lsn.grps.clr, common.grps.clr)

  lsn.grps=names(lsn.clrs)
  clrs=as.character(lsn.clrs)
  path.lsn.clr=path.lsn

  # to replace lesion groups with colors in the lsn matrix
  path.lsn.clr[path.lsn.clr %in% lsn.grps] <- clrs[match(path.lsn.clr, lsn.grps, nomatch = 0)]
  path.genes=selected.genes$ensembl.id
  path.gene.names=selected.genes$gene.name
  names(path.gene.names)=path.genes

  # to replace ensembl IDs with gene name
  rownames(path.expr)=path.gene.names[rownames(path.expr)]
  rownames(path.lsn)=path.gene.names[rownames(path.lsn)]

  # compute the distance between each two genes based on the lesion data using dist.lsn function
  dist.lsn.genes=dist.lsn(t(path.lsn))
  hcl.lsn.genes=stats::hclust(stats::as.dist(dist.lsn.genes),"complete")

  # compute distance between subjects based on lesions affecting pathway genes using dist.lsn function
  path.lsn.dist=dist.lsn(path.lsn)
  path.lsn.dist=as.matrix(path.lsn.dist)

  path.ls.hcl=stats::hclust(stats::as.dist(path.lsn.dist),method="ward.D2")
  sum.none=rowSums(path.lsn.clr=="gray")
  subj.ord=path.ls.hcl$order
  subj.labels=path.ls.hcl$labels

  #####################################
  # side-by-side heatmap
  opar<-graphics::par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  on.exit(par(opar), add = TRUE)
  plot(c(0,+1.1*ncol(path.lsn.clr)),
       c(0,-(nrow(path.lsn.clr)+1+nrow(path.expr))),
       type="n",xlab="",ylab="",
       axes=F)

  for (i in 1:nrow(path.lsn.clr))
    graphics::rect(0:(ncol(path.lsn.clr)-1),-(i-1),
         1:ncol(path.lsn.clr),-i,
         col=path.lsn.clr[hcl.lsn.genes$order[i],subj.ord],
         border=NA)
  graphics::text(ncol(path.lsn),-(1:nrow(path.lsn))+0.5,
       rownames(path.expr)[hcl.lsn.genes$order],
       pos=4,cex=0.75)

  # Add legend
  lsn.inc=names(lsn.clrs)%in%path.lsn
  graphics::legend("topright", inset=c(-0.25,0.01),
         fill=unlist(lsn.clrs[lsn.inc]),
         legend=names(lsn.clrs[lsn.inc]),
         cex=0.72,border=NA,bty="n")

  for (i in 1:nrow(path.expr))
  {
    y=path.expr[hcl.lsn.genes$order[i],]
    z=(y-mean(y))/stats::sd(y)
    clr=grDevices::rgb((z>0),0,(z<0),alpha=sqrt(1-exp(-abs(z))))
    graphics::rect(0:(ncol(path.lsn.clr)-1),-nrow(path.lsn.clr)-3-(i-1),
         1:ncol(path.lsn.clr),-nrow(path.lsn.clr)-3-i,
         col=clr[subj.ord],
         border=NA)
  }

  graphics::legend("bottomright", inset=c(-0.25,0.2),
         legend=c("Z.expr<0",0,"Z.expr>0"),
         fill = c("blue", "white", "red"),
         cex=0.72,border=NA,bty="n")

  graphics::text(ncol(path.expr),-nrow(path.lsn.clr)-3-1:nrow(path.expr)+0.5,
       rownames(path.expr)[hcl.lsn.genes$order],cex=0.75,pos=4)

  # Extract ordered lesion and expression data for the pathway genes
  pts.labels=as.data.frame(subj.labels)
  pts.labels<-tibble::rownames_to_column(pts.labels, "subj.ord")
  pts.order=as.data.frame(subj.ord)
  pts.order$index=1:nrow(pts.order)
  ordered.subj.final=merge(pts.labels,pts.order,by="subj.ord", all.y=TRUE)
  ordered.subj.final=ordered.subj.final[order(ordered.subj.final$index),]
  path.lsn.df=as.data.frame(path.lsn)
  ordered.lsn.data=path.lsn.df[ordered.subj.final$subj.labels]
  rownames(ordered.lsn.data) = paste(rownames(ordered.lsn.data),"_lsn")
  path.expr.df=as.data.frame(path.expr)
  ordered.expr.data=path.expr.df[ordered.subj.final$subj.labels]
  rownames(ordered.expr.data) = paste(rownames(ordered.expr.data),"_expr")
  ordered.path.data=rbind(ordered.lsn.data, ordered.expr.data)

  return(ordered.path.data)

}
