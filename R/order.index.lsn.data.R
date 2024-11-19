
#' Order Index Lesion Data
#'
#' @description
#' This function order and index lesion data by lesion type, the chromosome on which the lesion is located , and subject.
#'
#' @param lsn.data data.frame with lesion data prepared by the user in a GRIN compatible format. The data.frame should has five columns that include "ID" which is a column with id of the patient affected by the lesion, "chrom" which is the chromosome on which the lesion is located, "loc.start" which is the lesion start position, "loc.end" the lesion end position and "lsn.type" which is the lesion type for example gain, loss, mutation, fusion, etc...
#'
#' @return
#' A list with the following components:
#' \item{lsn.data}{Input lesion data}
#' \item{lsn.index}{data.frame with row start and row end for each type of lesions affecting each subject on a certain chromosome. For example, if a certain patient is affected by 1 deletion on chromosome 5, row start wil be equal to row end for loss on chromosome 5. However, if the patient is affected by 4 deletions, difference between row.start and row.end will be 3.}
#'
#' @export
#'
#' @references
#' Pounds, Stan, et al. (2013) A genomic random interval model for statistical analysis of genomic lesion data.
#'
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author {Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @examples
#' data(lesion.data)
#'
#' ordered.lsn=order.index.lsn.data(lesion.data)

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
