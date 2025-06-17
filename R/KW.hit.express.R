
#' Associate Lesion Groups with Gene Expression
#'
#' @description
#' Performs the Kruskal Wallis (KW) test to evaluate the association between lesion groups and the expression level of the corresponding gene.
#'
#' @param alex.data Output from the \code{\link{alex.prep.lsn.expr}} function. A list of three data frames:
#' \itemize{
#'   \item \code{row.mtch}: Table of matched lesion expression entries, including gene IDs.
#'   \item \code{alex.expr}: Gene expression matrix (rows = genes by Ensembl ID, columns = patient IDs).
#'   \item \code{alex.lsn}: Lesion status matrix with the same dimensions/order as \code{alex.expr}.
#' }
#'
#' @param gene.annotation Gene annotation table. Must be a \code{data.frame} with the following columns:
#' \code{gene} (Ensembl gene ID), \code{chrom} (chromosome), \code{loc.start} (start position), and \code{loc.end} (end position).
#' Can be obtained manually or via the \code{\link{get.ensembl.annotation}} function.
#'
#' @param min.grp.size Minimum number of patients required in a lesion group to be included in the test. For a gene to be tested, there must be at least two groups with more than \code{min.grp.size} patients each.
#'
#' @details
#' For each row in the \code{row.mtch} table, the function performs a Kruskal Wallis test comparing expression values of the gene across lesion groups.
#' The lesion groups are defined in the \code{alex.lsn} matrix. Patients with multiple types of lesions in a gene are assigned the label \code{"multiple"},
#' and those with no lesion are labeled \code{"none"}. The expression values are obtained from the corresponding row in \code{alex.expr}.
#'
#' The function tests whether expression levels significantly differ across lesion groups for the same gene.
#'
#' @return
#' A data frame where each row corresponds to a gene tested. Columns include:
#' \item{gene}{Ensembl gene ID.}
#' \item{gene.name}{HGNC gene symbol.}
#' \item{p.KW}{Kruskal Wallis test p-value.}
#' \item{q.KW}{FDR-adjusted q-value for multiple testing correction.}
#' \item{_n.subjects}{Number of subjects in each lesion group, including "none" and "multiple".}
#' \item{_mean}{Mean expression per lesion group.}
#' \item{_median}{Median expression per lesion group.}
#' \item{_sd}{Standard deviation of expression per lesion group.}
#'
#' @export
#'
#' @importFrom stats kruskal.test p.adjust
#'
#' @references
#' Hollander, M., & Wolfe, D. A. (1973). Nonparametric Statistical Methods. New York: Wiley. pp. 115-120.
#'
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link{alex.prep.lsn.expr}}
#'
#' @examples
#' data(expr_data)
#' data(lesion_data)
#' data(hg38_gene_annotation)
#'
#' # Prepare matched lesion-expression data
#' alex.data <- alex.prep.lsn.expr(expr_data, lesion_data,
#'                                 hg38_gene_annotation,
#'                                 min.expr = 1, min.pts.lsn = 5)
#'
#' # Perform Kruskal Wallis test between lesion groups and expression levels:
#' alex.kw.results <- KW.hit.express(alex.data,
#'                                   hg38_gene_annotation,
#'                                   min.grp.size = 5)

KW.hit.express=function(alex.data,          # output of the alex.prep.lsn.expr function (list of three data tables "alex.expr" with expression data ready for KW test, "alex.lsn" with lesion data and row.mtch)
                        gene.annotation,    # gene annotation data. First column "gene" should has ensembl IDs
                        min.grp.size=NULL)  # minimum group size to perform the test (there should be at least two groups with number of patients > min.grp.size)

{
  expr.matrix=as.matrix(alex.data$alex.expr)
  lsn.matrix=as.matrix(alex.data$alex.lsn)
  row.mtch=alex.data$alex.row.mtch

  message(paste0("Computing KW P-value: ",date()))
  p=apply(row.mtch,1,one.KW.pvalue,
          expr.mtx=expr.matrix,
          hit.grps=lsn.matrix,
          min.grp.size=min.grp.size)

  message(paste0("Preparing results table: ",date()))
  kw.res=cbind(row.mtch,p.KW=p)
  # To prepare and print KW test results

  colnames(kw.res)=c("expr.row", "gene", "p.KW")
  kw.res.annotated=merge(gene.annotation,kw.res,by="gene", all.y=TRUE)

  # Compute FDR adjusted q values
  kw.res.annotated$q.KW=stats::p.adjust(kw.res.annotated$p.KW,method="fdr")

  # To compute number of subjects, mean, median and stdev by lesion type
  unique.grps=sort(unique(lsn.matrix))
  unique.grps=sort(unique(unique.grps))
  count.by.lesion <- sapply(unique.grps,function(x)rowSums(lsn.matrix==x))
  colnames(count.by.lesion) = paste(colnames(count.by.lesion),"n.subjects",sep="_")
  count.by.lesion=as.matrix(count.by.lesion)

  mean.by.lsn=row.stats.by.group(expr.matrix, lsn.matrix, mean)
  colnames(mean.by.lsn) = paste(colnames(mean.by.lsn),"mean",sep="_")
  mean.by.lsn=as.matrix(mean.by.lsn)

  median.by.lsn=row.stats.by.group(expr.matrix, lsn.matrix, median)
  colnames(median.by.lsn) = paste(colnames(median.by.lsn),"median",sep="_")
  median.by.lsn=as.matrix(median.by.lsn)

  sd.by.lsn=row.stats.by.group(expr.matrix, lsn.matrix, sd)
  colnames(sd.by.lsn) = paste(colnames(sd.by.lsn),"sd",sep="_")
  sd.by.lsn=as.matrix(sd.by.lsn)

  kw.res.final=cbind(kw.res.annotated, count.by.lesion, mean.by.lsn, median.by.lsn, sd.by.lsn)
  kw.res.final=kw.res.final[!duplicated(kw.res.final[ , "gene.name"]),]
  res=kw.res.final

  return(res)

}
