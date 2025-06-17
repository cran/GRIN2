
#' Order and Index Lesion Data
#'
#' @description
#' This function orders and indexes lesion data by lesion type, chromosome, and subject ID. It prepares lesion data for downstream GRIN analysis by structuring it in a way that facilitates efficient access and overlap computations.
#'
#' @param lsn.data A `data.frame` containing lesion data formatted for GRIN. It must include the following five columns:
#' \describe{
#'   \item{"ID"}{Patient identifier.}
#'   \item{"chrom"}{Chromosome on which the lesion is located.}
#'   \item{"loc.start"}{Start position of the lesion.}
#'   \item{"loc.end"}{End position of the lesion.}
#'   \item{"lsn.type"}{Lesion type (e.g., gain, loss, mutation, fusion, etc...).}
#' }
#'
#' @return
#' A list with two elements:
#' \item{lsn.data}{The input lesion data, ordered by lesion type, chromosome, and subject.}
#' \item{lsn.index}{A `data.frame` with two columns, `row.start` and `row.end`, indicating the index range of lesions for each subject-lesion type-chromosome combination. For example, if a patient has a single deletion on chromosome 5, `row.start` will equal `row.end`. If there are four deletions, the range will span four rows.}
#'
#' @export
#'
#' @references
#' Pounds, S., et al. (2013). A genomic random interval model for statistical analysis of genomic lesion data.
#'
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @examples
#' data(lesion_data)
#'
#' ordered.lsn <- order.index.lsn.data(lesion_data)

order.index.lsn.data=function(lsn.data)  # lesion data provided by the user in a GRIN compatible format

{
  l=nrow(lsn.data)
  lsn.ord=order(lsn.data[,"lsn.type"],
                lsn.data[,"chrom"],
                lsn.data[,"ID"])
  lsn.data=lsn.data[lsn.ord,]

  lsn.chng=which((lsn.data[-1,"lsn.type"]!=lsn.data[-l,"lsn.type"])|
                   (lsn.data[-1,"chrom"]!=lsn.data[-l,"chrom"])|
                   (lsn.data[-1,"ID"]!=lsn.data[-l,"ID"]))
  lsn.start=c(1,lsn.chng+1)
  lsn.end=c(lsn.chng,l)
  lsn.index=cbind.data.frame(lsn.type=lsn.data[lsn.start,"lsn.type"],
                             chrom=lsn.data[lsn.start,"chrom"],
                             ID=lsn.data[lsn.start,"ID"],
                             row.start=lsn.start,
                             row.end=lsn.end)
  lsn.data$lsn.row=1:l

  res=list(lsn.data=lsn.data,
           lsn.index=lsn.index)
  return(res)
}
