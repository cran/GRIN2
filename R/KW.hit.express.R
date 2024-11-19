
#' Associate Lesion with Expression Data Using Kruskal-Wallis Test
#'
#' @description
#' Function uses Kruskal-Wallis test to evaluate the association between lesion groups and expression level of the same corresponding gene.
#'
#' @param alex.data output of the alex.prep.lsn.expr function. It's a list of three data tables that include "row.mtch", "alex.expr" with expression data, "alex.lsn" with lesion data. Rows of alex.expr, and "alex.lsn" matrices are ordered by gene ensembl IDs and columns are ordered by patient ID.
#' @param gene.annotation Gene annotation data either provided by the user or retrieved from ensembl BioMart database using get.ensembl.annotation function included in the GRIN2.0 library. Data.frame should has four columns: "gene" which is the ensembl ID of annotated genes, "chrom" which is the chromosome on which the gene is located, "loc.start" which is the gene start position, and "loc.end" the gene end position.
#' @param min.grp.size Minimum number of subjects in a lesion group to be included in the KW test (there should be at least two groups with number of patients > min.grp.size) to run the KW test for a certain gene.
#'
#' @details
#' The function uses the ensembl IDs in each row of the row.mtch file and run the Kruskal-Wallis test for association between lesion groups of the gene in the "hit.row" column with expression level of the gene in the "expr.row" column. IDs in the two columns should be the same if the KW test will be used to evaluate association between lesion groups and expression level of the same corresponding gene. If the same patient is affected with multiple types of lesions in the same gene for example gain AND mutations, the entry will be denoted as "multiple" and patients without any type of lesions will be coded as "none".
#'
#' @return
#' A data table with multiple columns that include:
#' \item{gene}{ensembl ID of the gene of interest.}
#' \item{gene.name}{Gene name of the gene of interest.}
#' \item{p.KW}{Kruskal-Wallis test p-value.}
#' \item{q.KW}{Kruskal-Wallis test FDR adjusted q-value.}
#' \item{_n.subjects}{Multiple columns with number of subjects with each type of lesion affecting the gene, number of subjects without any lesion and number of subjects with multiple types of lesions.}
#' \item{_mean}{Multiple columns with mean expression level of the gene in subjects with each type of lesion, mean expression in subjects without any lesion and mean expression in subjects with multiple types of lesions.}
#' \item{_median}{Multiple columns with median expression of the gene in subjects with each type of lesion, median expression in subjects without any lesion and median expression in subjects with multiple types of lesions.}
#' \item{_sd}{Multiple columns with standard deviation of the expression level of the gene in subjects with each type of lesion, standard deviation in subjects without any lesion and standard deviation in subjects with multiple types of lesions.}
#'
#' @export
#'
#' @importFrom stats kruskal.test p.adjust
#'
#' @references
#' Myles Hollander and Douglas A. Wolfe (1973), Nonparametric Statistical Methods. New York: John Wiley & Sons. Pages 115–120.
#'
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author {Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [alex.prep.lsn.expr()]
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
#' # run Kruskal-Wallis test for association between lesion groups and expression level of the
#' # same corresponding gene:
#' alex.kw.results=KW.hit.express(alex.data, hg19.gene.annotation, min.grp.size=5)

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
