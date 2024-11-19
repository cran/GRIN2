
#' Prepare Lesion and Expression Data for Waterfall Plots
#'
#' @description
#' Function prepares lesion and expression data of a selected gene for the alex.waterfall.plot function.
#'
#' @param alex.data output of the alex.prep.lsn.expr function. It's a list of three data tables that include "row.mtch", "alex.expr" with expression data, "alex.lsn" with lesion data. Rows of alex.expr, and "alex.lsn" matrices are ordered by the gene ensembl IDs and columns are ordered by patient IDs.
#' @param alex.kw.results ALEX Kruskal-Wallis test results (output of the KW.hit.express function).
#' @param gene Gene name or ensembl ID of the gene of interest.
#' @param lsn.data Lesion data in a GRIN compatible format. Object should has five columns that include "ID" with patient ID, "chrom" which is the chromosome on which the lesion is located, "loc.start" which is the lesion start position, "loc.end" the lesion end position and "lsn.type" which is the lesion category for example gain, loss, mutation, fusion, etc...
#'
#' @details
#' Function prepares lesion and expression data of a selected gene for the alex.waterfall.plot function. It return a data table with patient ID, lesion types that affect each patient if any and expression level of the gene of interest. It also extract the kruskal-wallis test result and all lesions that affect the gene of interest.
#'
#' @return
#' A list of four components:
#' \item{gene.lsn.exp}{Data table with three columns("ID" with patient ID, "gene.name_lsn" has the type of lesion affecting the patient which can be none, gain, mutation, multiple, etc.. and "gene.name_expr" which has the expression level of the gene of interst in this particular patient.}
#' \item{lsns}{Data table with all lesions affecting the gene of interst in a GRIN compatible format extracted from the lesion data file.}
#' \item{stats}{One row with the Kruskal-Wallis test result for the gene of interst (output of the KW.hit.express function).}
#' \item{gene.ID}{Gene name of the gene of interst.}
#'
#' @export
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author {Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [alex.prep.lsn.expr()], [KW.hit.express()]
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

alex.waterfall.prep=function(alex.data,             # output of the alex.prep.lsn.expr function (list of three data table "alex.expr" with expression data, "alex.lsn" with lesion data and row.mtch)
                             alex.kw.results,       # ALEX Kruskal-Wallis test results (output of the KW.hit.express function)
                             gene,                  # gene name or ensembl.ID
                             lsn.data)              # lesion data in a GRIN compatible format


{
  # Find the row of the KW test with results for the gene
  int.row=which((alex.kw.results$gene==gene)|
                  (alex.kw.results$gene.name==gene))

  if (length(int.row)!=1) stop("gene must match exactly ONE row of alex.kw.results")

  ens.ID=alex.kw.results$gene[int.row]
  gene.ID=alex.kw.results$gene.name[int.row]

  gene.lsn.exp=as.data.frame(colnames(alex.data$alex.lsn))
  colnames(gene.lsn.exp)="ID"
  # Add the lesion data for this gene to the gene.lsn.exp data
  rownames(gene.lsn.exp)=gene.lsn.exp$ID
  gene.lsn=paste0(gene.ID,".lsn")
  lsn.row=which(rownames(alex.data$alex.lsn)==ens.ID)
  gene.lsn.exp[,gene.lsn]=""
  lsn.ids=intersect(gene.lsn.exp$ID,colnames(alex.data$alex.lsn))
  gene.lsn.exp[lsn.ids,gene.lsn]=unlist(alex.data$alex.lsn[lsn.row,lsn.ids])

  # Add the expression data for this gene to the gene.lsn.exp data
  gene.rna=paste0(gene.ID,".RNA")
  gene.lsn.exp[,gene.rna]=NA
  rna.ids=intersect(gene.lsn.exp$ID,colnames(alex.data$alex.expr))
  rna.row=which(rownames(alex.data$alex.expr)==ens.ID)
  gene.lsn.exp[rna.ids,gene.rna]=unlist(alex.data$alex.expr[rna.row,rna.ids])

  # Find lesions that overlap this gene
  lsn.mtch=(lsn.data$chrom==alex.kw.results$chrom[int.row])&
    (lsn.data$loc.end>=alex.kw.results$loc.start[int.row])&
    (lsn.data$loc.start<=alex.kw.results$loc.end[int.row])

  lsns=lsn.data[lsn.mtch,]

  res=list(gene.lsn.exp=gene.lsn.exp,       # data table with three columns that has patient ID, type of lesions that affect this gene if any, expression level of the gene
           lsns=lsns,                       # all lesions that affect this gene of interest
           stats=alex.kw.results[int.row,], # KW test results for the gene of interest
           gene.ID=gene.ID)                 # gene name
}
