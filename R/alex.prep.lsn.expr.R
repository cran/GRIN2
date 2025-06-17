
#' Prepare Lesion and Expression Data for Kruskal Wallis Test
#'
#' @description
#' Prepares matched lesion and expression data matrices for use with the \code{KW.hit.express} function, which performs Kruskal Wallis tests to assess associations between lesion groups and gene expression levels.
#'
#' @param expr.mtx A data frame or matrix of normalized, log2-transformed gene expression values with genes in rows and subjects in columns. The first column must be named \code{"ensembl.ID"} and contain Ensembl gene IDs.
#' @param lsn.data A data frame of lesion data in GRIN-compatible format. Must contain five columns:
#' \code{"ID"} (patient ID), \code{"chrom"} (chromosome), \code{"loc.start"} (lesion start position), \code{"loc.end"} (lesion end position), and \code{"lsn.type"} (lesion type; e.g., gain, loss, mutation, fusion, etc..).
#' @param gene.annotation A gene annotation data frame, either user-provided or retrieved via \code{get.ensembl.annotation()} from the GRIN2.0 package. Must contain four columns: \code{"gene"} (Ensembl gene ID), \code{"chrom"} (chromosome), \code{"loc.start"} (gene start), and \code{"loc.end"} (gene end).
#' @param min.expr Minimum total expression level required for a gene to be retained (i.e., sum of expression values across all subjects). Useful to filter out genes with very low expression.
#' @param min.pts.lsn Minimum number of subjects required to have a lesion in a gene for that gene to be retained for the KW test.
#'
#' @details
#' The function uses \code{prep.lsn.type.matrix()} internally to create a lesion matrix where each gene is represented by one row and all lesion types are included. It filters genes to retain only those with both sufficient expression and lesion data. The final expression and lesion matrices are matched by gene and patient IDs, with rows ordered by Ensembl gene ID and columns by patient ID.
#'
#' @return
#' A list with the following components:
#' \item{alex.expr}{A matrix of gene expression data with Ensembl gene IDs as row names and patient IDs as column names.}
#' \item{alex.lsn}{A matrix of lesion data for the same genes and patients as in \code{alex.expr}, similarly ordered.}
#' \item{alex.row.mtch}{A data frame with two columns showing the matched Ensembl gene IDs from the expression and lesion matrices.}
#'
#' @export
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link{KW.hit.express}}
#'
#' @examples
#' data(expr_data)
#' data(lesion_data)
#' data(hg38_gene_annotation)
#'
#' # Prepare matched lesion and expression data
#' alex.data <- alex.prep.lsn.expr(expr_data, lesion_data,
#'                                 hg38_gene_annotation, min.expr = 1,
#'                                 min.pts.lsn = 5)

alex.prep.lsn.expr=function(expr.mtx,          # Normalized log2 transformed expression data with genes in rows and subjects in columns (first column "ensembl.ID" should be gene ensembl IDs)
                            lsn.data,          # lesion data in GRIN compatible format. Data frame should has five columns that include "ID" with patient ID, "chrom" which is the chromosome on which the lesion is located, "loc.start" which is the lesion start position, "loc.end" the lesion end position and "lsn.type" which is the lesion type for example gain, loss, mutation, fusion, etc...
                            gene.annotation,   # gene annotation data either provided by the user or directly retreived from ensembl biomaRT database using get.ensembl.annotation function included in the GRIN2.0 library if the genome.version is specified. Object should has four columns "gene" which is the ensembl ID of annotated genes, "chrom" which is the chromosome on which the gene is located, "loc.start" which is the gene start position, and "loc.end" the gene end position.
                            min.expr=NULL,     # minimum allowed expression level of the gene (the sum of expression level of the gene in all patients; helpful to exclude genes with very low expression)
                            min.pts.lsn=NULL)  # minimum number of patients with any type of lesions in a certain gene

{
  gene.lsn=prep.gene.lsn.data(lsn.data, gene.annotation)    # Prepare gene and lesion data for later computations
  gene.lsn.overlap= find.gene.lsn.overlaps(gene.lsn)        # Use the results of prep.gene.lsn.data to find lesion-gene overlaps

  message(paste0("Preparing gene-lesion matrix: ",date()))
  lsn.type.matrix=prep.lsn.type.matrix(gene.lsn.overlap)    # prepare lesion type matrix (just one row for each gene with all lesion types included)
  lsn.grp.mtx=as.data.frame(lsn.type.matrix)
  id.hits = colnames(lsn.grp.mtx)

  # prepare expression data for ALEX analysis
  expr.mtx=expr.mtx
  rownames(expr.mtx)=expr.mtx[,1]
  expr.mtx=expr.mtx[,-1]

  id.expr=colnames(expr.mtx)
  lsn.grp.mtx=lsn.grp.mtx[ ,(names(lsn.grp.mtx) %in% id.expr)]

  id.lsn.final=colnames(lsn.grp.mtx)
  expr.mtx=expr.mtx[ ,(names(expr.mtx) %in% id.lsn.final)]

  if (is.numeric(min.expr))
  {
    # To keep only genes with total expression level > value specified by the user in the min.expr argument
    total.expr=rowSums(expr.mtx)
    expr.grp.keep=total.expr>=min.expr
    expr.mtx=expr.mtx[expr.grp.keep,]
  }
  ############################
  # To exclude any gene that has no lesions in at least the number of patients specified in min.pts.lsn
  if (is.numeric(min.pts.lsn))
  {
    n.none=rowSums(lsn.grp.mtx== "none")
    lsn.grp.keep=(ncol(lsn.grp.mtx)-n.none)>=min.pts.lsn
    lsn.grp.mtx=lsn.grp.mtx[lsn.grp.keep,]
  }
  # To keep only shared set of genes between expression and lesion matrices
  lsn.genes=rownames(lsn.grp.mtx)
  expr.mtx=expr.mtx[rownames(expr.mtx) %in% lsn.genes,]

  expr.genes=rownames(expr.mtx)
  lsn.grp.mtx=lsn.grp.mtx[rownames(lsn.grp.mtx) %in% expr.genes,]

  ## order lsn matrix by gene ID first then colnames for patient IDs alphabetically
  lsn.grp.mtx=lsn.grp.mtx[order(rownames(lsn.grp.mtx)),]
  lsn.grp.mtx=lsn.grp.mtx[,order(colnames(lsn.grp.mtx))]

  ## order expr matrix by gene ID first then colnames for patient IDs alphabetically
  expr.mtx=expr.mtx[order(rownames(expr.mtx)),]
  expr.mtx=expr.mtx[,order(colnames(expr.mtx))]

  # To make sure that both lesion and expression matrices have same patients and genes order
  check.rownames=all(rownames(lsn.grp.mtx)==rownames(expr.mtx))
  if (check.rownames==FALSE)
    stop("Gene-lesion matrix rownames must match expression matrix rownames.")

  check.colnames=all(colnames(lsn.grp.mtx)==colnames(expr.mtx))
  if (check.colnames==FALSE)
    stop("Gene-lesion matrix patient IDs must match expression matrix patient IDs.")

  # To generate row.mtch file
  lsn.ensembl.id=rownames(lsn.grp.mtx)
  expr.ensembl.id=rownames(expr.mtx)

  row.mtch=cbind.data.frame(expr.row=expr.ensembl.id,
                            hit.row=lsn.ensembl.id)


  res=list(alex.expr=expr.mtx,      # expression data (data table with gene ensembl IDs as row names and each column is a patient)
           alex.lsn=lsn.grp.mtx,    # overlapped gene lesion data table with gene ensembl IDs as row names and each column is a patient
           alex.row.mtch=row.mtch)  # ensembl ID for one row of the alex.lsn table to be paired with one row of the expression data for the Kruskal Wallis test

  return(res)

}
