
#' Prepare Binary Lesion Matrix
#'
#' @description
#' Prepares a lesion matrix with each gene affected by a certain type of lesion as a row and each patient as a column.
#'
#' @param ov.data list of six data.frames that represent the output results of the find.gene.lsn.overlaps function.
#' @param min.ngrp if specified, rows with number of patients affected by a specific type of lesion that's less than the specified number will be discarded (default is 0; function will return all genes affected by a lesion in at least one patient), for example if only one patient is affected by gain in MYB gene.
#'
#' @details
#' The function uses the output results of the find.gene.lsn.overlaps function and create a binary lesion matrix with each gene affected by certain lesion type as a row and each patient as a column. Rownames are labelled as gene.ID_lesion.type (for example: ENSG00000118513_gain for gains affecting MYB gene). The entry for each patient in the table will be denoted as 1 if the patient is affected by this specific type of lesion in the gene, for example gain in MYB gene (ENSG00000118513) or 0 otherwise.
#'
#' @return
#' The function returns a binary lesion matrix with each row labelled as gene.ID_lesion.type and each column is a patient. Entry for each patient in the table will be denoted as 1 if the gene is affected by this specific type of lesion or 0 otherwise.
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
#' @seealso [prep.gene.lsn.data()], [find.gene.lsn.overlaps()]
#'
#' @examples
#' data(lesion.data)
#' data(hg19.gene.annotation)
#'
#' # prepare gene and lesion data for later computations:
#' prep.gene.lsn=prep.gene.lsn.data(lesion.data,
#'                                  hg19.gene.annotation)
#'
#' # determine lesions that overlap each gene (locus):
#' gene.lsn.overlap=find.gene.lsn.overlaps(prep.gene.lsn)
#'
#' # prepare the lesion binary matrix with a minimum of 5 patients affected by the lesion to be
#' # included in the final matrix:
#' lsn.binary.mtx=prep.binary.lsn.mtx(gene.lsn.overlap, min.ngrp=5)

prep.binary.lsn.mtx=function(ov.data,     # output of the find.gene.lsn.overlaps function
                             min.ngrp=0)  # if specified, rows with number of patients affected by this type of lesion that's less than the specified number will be discarded (default is 0; will return all genes affected by this specified type of lesion in at least one patient).

{
  gene.lsn=ov.data$gene.lsn.hits

  # Order and index data by gene and lesion type
  gene.lsn.type=paste0(gene.lsn$gene,"_",gene.lsn$lsn.type)
  ord=order(gene.lsn$gene,gene.lsn$lsn.type)
  gene.lsn=gene.lsn[ord,]
  m=nrow(gene.lsn)
  new.gene.lsn=which((gene.lsn$gene[-1]!=gene.lsn$gene[-m])|
                       (gene.lsn$lsn.type[-1]!=gene.lsn$lsn.type[-m]))

  row.start=c(1,new.gene.lsn+1)
  row.end=c(new.gene.lsn,m)
  gene.lsn.indx=gene.lsn$gene.lsn[row.start]
  gene.lsn.index=cbind.data.frame(gene=gene.lsn$gene[row.start],
                                  lsn.type=gene.lsn$lsn.type[row.start],
                                  row.start=row.start,
                                  row.end=row.end)
  k=nrow(gene.lsn.index)
  uniq.ID=unique(gene.lsn$ID)
  n=length(uniq.ID)
  gene.lsn.mtx=matrix(0,k,n)
  colnames(gene.lsn.mtx)=as.character(uniq.ID)
  rownames(gene.lsn.mtx)=paste0(gene.lsn.indx$gene,"_",
                                gene.lsn.index$lsn.type)

  for (i in 1:k)
  {
    rows=(gene.lsn.index$row.start[i]:gene.lsn.index$row.end[i])
    ids=as.character(gene.lsn$ID[rows])
    gene.lsn.mtx[i,ids]=1
  }

  n.hit=rowSums(gene.lsn.mtx)
  min.n=pmin(n.hit,n-n.hit)

  keep.row=which(min.n>=min.ngrp)
  gene.lsn.mtx=matrix(gene.lsn.mtx[keep.row,],
                      length(keep.row),n)
  colnames(gene.lsn.mtx)=as.character(uniq.ID)
  rownames(gene.lsn.mtx)=paste0(gene.lsn.index$gene,"_",
                                gene.lsn.index$lsn.type)[keep.row]


  return(gene.lsn.mtx)
}
