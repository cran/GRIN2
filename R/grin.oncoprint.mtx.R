
#' GRIN OncoPrint Matrix
#'
#' @description
#' Function use GRIN results table and prepare the lesion matrix that the user can pass to the oncoprint function from ComplexHeatmap package to geneate an OncoPrint for a selcted list of genes.
#'
#' @param grin.res GRIN results (output of the grin.stats function).
#' @param oncoprint.genes Vector of ensembl IDs for the selected list of genes to be added to the OncoPrint.
#'
#' @details
#' Function will use the input list of ensembl IDs to prepare a data table of lesions that affect these genes (each row is a gene and each column is a patient ID). This lesion matrix is compatible and can be passed to oncoprint function in ComplexHeatmap library to prepare an OncoPrint for lesions in the selected list of genes.
#'
#' @return
#' Function uses the output results of grin.stats function and return data table of lesions that affect a group of selected genes (each row is a gene and each column is a patient ID).
#'
#' @export
#'
#' @importFrom tibble rownames_to_column
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author {Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [grin.stats()]
#'
#' @examples
#' data(lesion.data)
#' data(hg19.gene.annotation)
#' data(hg19.chrom.size)
#'
#' # Run GRIN analysis using grin.stats function:
#' grin.results=grin.stats(lesion.data,
#'                         hg19.gene.annotation,
#'                         hg19.chrom.size)
#'
#' # specify a list of genes to be included in the oncoprint (driver genes):
#' oncoprint.genes=as.vector(c("ENSG00000148400", "ENSG00000171862", "ENSG00000171843",
#'                             "ENSG00000156531", "ENSG00000162367", "ENSG00000096968",
#'                             "ENSG00000105639", "ENSG00000118513","ENSG00000102974",
#'                             "ENSG00000133703"))
#'
#' # prepare the oncoprint lesion matrix:
#' oncoprint.mtx=grin.oncoprint.mtx(grin.results,
#'                                  oncoprint.genes)
#'
#' # user can also specify a list of top significant genes in the GRIN constellation test:
#' # for example: select genes affected by two types of lesion with q2.nsubj<0.01:
#' genes.const = grin.results$gene.hits[grin.results$gene.hits$q2.nsubj < 0.01, ]
#' # get ensembl.ids for this list of genes
#' selected.genes=as.vector(genes.const$gene)
#' oncoprint.mtx.const=grin.oncoprint.mtx(grin.results,
#'                                        selected.genes)

grin.oncoprint.mtx=function(grin.res, # GRIN results (output of the grin.stats function)
                            oncoprint.genes) # vector of ensembl IDs for the selected list of genes

{
  selected=unlist(oncoprint.genes)
  selected=as.vector(selected)
  selected.genes= grin.res$gene.lsn.data[grin.res$gene.lsn.data$gene %in% selected,]
  selected.genes=selected.genes[,c(2,7,11)]  # extract patient IDs and lsn type for each gene in the selected genes list
  row.data=paste(selected.genes[,1],
                 selected.genes[,2],
                 selected.genes[,3],
                 sep="_")
  dup.data=duplicated(row.data)
  select.genes=selected.genes[!dup.data,]

  ord=order(select.genes$gene,
            select.genes$ID,
            select.genes$lsn.type)
  select.genes=select.genes[ord,]

  uniq.genes=unique(select.genes$gene)
  uniq.subj=unique(select.genes$ID)
  n.genes=length(uniq.genes)
  n.subj=length(uniq.subj)
  mtx=matrix("",n.genes,n.subj)   # create a matrix with each gene as a row
  colnames(mtx)=uniq.subj
  rownames(mtx)=uniq.genes

  k=nrow(select.genes)
  for (i in 1:k)
  {
    subj.id=select.genes[i,"ID"]
    gene.id=select.genes[i,"gene"]
    mtx[gene.id,subj.id]=paste0(mtx[gene.id,subj.id],
                                select.genes[i,"lsn.type"],";")
  }
  mtx=as.data.frame(mtx)
  mtx<-tibble::rownames_to_column(mtx, "ensembl.ID")

  gene.annotation= grin.res$gene.data
  ensembl.annotation=cbind(gene.annotation$gene, gene.annotation$gene.name)
  colnames(ensembl.annotation)=c("ensembl.ID", "gene.name")

  mtx.final=merge(ensembl.annotation,mtx,by="ensembl.ID", all.y=TRUE)  # add gene name
  mtx.final=mtx.final[,-1]
  rownames(mtx.final)=mtx.final[,1]
  mtx.final=mtx.final[,-1]

  return(mtx.final)

}
