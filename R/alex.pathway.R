
#' Visualize Lesion and Expression Data by Pathway
#'
#' @description
#' Computes pairwise distances between subjects based on lesion profiles in genes associated with a specified pathway, and returns a figure with two panels: one showing lesion data and another showing expression data, both ordered based on the computed distances (useful for hierarchical clustering). It also returns the corresponding ordered data.
#'
#' @param alex.data Output of the \code{alex.prep.lsn.expr} function. A list of three data tables: \code{"row.mtch"}, \code{"alex.expr"} (expression matrix), and \code{"alex.lsn"} (lesion matrix). Rows in both matrices are ordered by Ensembl gene IDs; columns represent patient IDs.
#' @param lsn.clrs Optional. A named vector of colors for lesion types. If not specified, default colors will be assigned using \code{default.grin.colors()}.
#' @param lsn.data A data frame of lesion data in GRIN-compatible format with the following columns:
#' \code{"ID"} (patient ID), \code{"chrom"} (chromosome), \code{"loc.start"} (lesion start), \code{"loc.end"} (lesion end), and \code{"lsn.type"} (e.g., gain, loss, mutation, fusion).
#' @param pathways A data frame with three columns: \code{"gene.name"} (gene symbol), \code{"ensembl.id"} (Ensembl gene ID), and \code{"pathway"} (pathway name).
#' @param selected.pathway A character string indicating the pathway of interest.
#'
#' @details
#' This function identifies all genes associated with the specified pathway, extracts lesion and expression data for those genes, and computes pairwise distances between subjects based on their lesion profiles. It uses hierarchical clustering to order the subjects and visualizes lesion and expression matrices in two aligned panels. It also returns a data frame containing the ordered expression and lesion data for all genes in the pathway.
#'
#' @return
#' A list with the following element:
#' \item{ordered.path.data}{A data frame with lesion and expression data for pathway genes, ordered according to hierarchical clustering (matching the order used in the plot).}
#'
#' A figure with two panels will also be generated showing:
#' \enumerate{
#'   \item Lesion data of pathway genes across subjects
#'   \item Expression data of the same genes across subjects
#' }
#' Both panels are ordered by subject similarity based on lesion profiles.
#'
#' @export
#'
#' @importFrom tibble rownames_to_column
#' @importFrom stats hclust as.dist sd
#' @importFrom grDevices rgb
#' @importFrom graphics par rect legend text
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link{alex.prep.lsn.expr}}, \code{\link[stats]{hclust}}
#'
#' @examples
#' data(expr_data)
#' data(lesion_data)
#' data(hg38_gene_annotation)
#' data(pathways)
#'
#' # Prepare matched expression and lesion data
#' alex.data <- alex.prep.lsn.expr(expr_data, lesion_data,
#'                                 hg38_gene_annotation, min.expr = 5, min.pts.lsn = 5)
#'
#' # Visualize pathway-level association using JAK pathway as an example
#' alex.path=alex.pathway(alex.data,
#'                        lsn.data = lesion_data,
#'                        pathways = pathways,
#'                        selected.pathway = "Jak_Pathway")
#'
#' # Access the ordered data matrix used in the plot
#' head(alex.path$ordered.path.data)

alex.pathway=function(alex.data,          # output of the alex.prep.lsn.expr function (list of three data table "alex.expr" with expression data ready for KW test, "alex.lsn" with overlapped gene lesion data and row.mtch)
                      lsn.clrs=NULL,      # Specified colors per lesion types (gene plots when gene name is specified). If not specified, colors will be automatically assigned using default.grin.colors function
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
  lesions=lsn.data
  unique.grps=sort(unique(lesions$lsn.type))
  if (is.null(lsn.clrs))
  {
    lsn.grps.clr=default.grin.colors(unique.grps)
    common.grps.clr=c(none="gray", multiple="violet")
    lsn.clrs=c(lsn.grps.clr, common.grps.clr)
  }

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
