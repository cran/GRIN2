
#' Prepare Lesion and Expression Data for Waterfall Plots
#'
#' @description
#' Prepares matched lesion and expression data for a selected gene to be used with the \code{alex.waterfall.plot} function.
#'
#' @param alex.data Output from \code{alex.prep.lsn.expr}. A list of three data tables: \code{"row.mtch"} (gene matching info), \code{"alex.expr"} (expression matrix), and \code{"alex.lsn"} (lesion matrix). Rows are ordered by Ensembl gene IDs, and columns are ordered by patient IDs.
#' @param alex.kw.results Kruskal Wallis test results for gene expression by lesion group, as returned by \code{KW.hit.express}.
#' @param gene Gene of interest, specified by either its gene symbol or Ensembl ID.
#' @param lsn.data Lesion data in GRIN-compatible format. A data frame with five required columns: \code{"ID"} (patient ID), \code{"chrom"} (chromosome), \code{"loc.start"} (lesion start), \code{"loc.end"} (lesion end), and \code{"lsn.type"} (lesion type; e.g., gain, loss, mutation, fusion, etc...).
#'
#' @details
#' This function extracts and combines lesion and expression data for a specified gene across patients. It returns a data table showing each patient's lesion status and expression level for the gene. It also extracts the corresponding Kruskal Wallis test result and all lesions that affect the gene from the lesion data.
#'
#' @return
#' A list with the following components:
#' \item{gene.lsn.exp}{A data table with three columns: \code{"ID"} (patient ID), \code{"<gene_name>_lsn"} (lesion status: e.g., none, gain, mutation, multiple), and \code{"<gene_name>_expr"} (expression level of the gene in that patient).}
#' \item{lsns}{A data table of all lesions affecting the gene of interest, extracted from the input lesion data (GRIN-compatible format).}
#' \item{stats}{A one-row data frame containing the Kruskal Wallis test result for the gene, from \code{KW.hit.express}.}
#' \item{gene.ID}{The gene name (symbol or Ensembl ID) provided as input.}
#'
#' @export
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link{alex.prep.lsn.expr}}, \code{\link{KW.hit.express}}
#'
#' @examples
#' data(expr_data)
#' data(lesion_data)
#' data(hg38_gene_annotation)
#'
#' # Prepare matched expression and lesion data
#' alex.data <- alex.prep.lsn.expr(expr_data, lesion_data,
#'                                 hg38_gene_annotation, min.expr = 1, min.pts.lsn = 5)
#'
#' # Run Kruskal Wallis test
#' alex.kw.results <- KW.hit.express(alex.data, hg38_gene_annotation, min.grp.size = 5)
#'
#' # Prepare lesion and expression data for waterfall plot of WT1
#' WT1.waterfall.prep <- alex.waterfall.prep(alex.data, alex.kw.results, "WT1", lesion_data)

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
