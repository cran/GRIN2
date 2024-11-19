
#' Prepare Lesion Type Matrix
#'
#' @description
#' The function prepare a lesion matrix with all types of lesions affecting certain gene as a row and each patient as a column.
#'
#' @param ov.data list of six data.frames that represent the output results of the find.gene.lsn.overlaps function.
#' @param min.ngrp if specified, rows with number of patients affected by all different types of lesions that's less than the specified number will be discarded (default is 0; will return all genes affected by any type of lesions in at least one patient).
#'
#' @details
#' The function returns a lesion matrix with each row as a gene and each column is a patient. If a gene is affected by one type of lesions in a certain patient, the entry will be labelled by lesion type (for example: gain OR mutation). However, if the same gene is affected by more than one type of lesions in a certain patient (for example: gain AND mutation), the entry will be labelled as "multiple". If the gene is not affected by any lesion, the entry for this patient will be labelled as "none".
#'
#' @return
#' The function returns a lesion matrix with all types of lesions affecting certain gene as a row and each patient as a column.
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
#' # prepare the lesion matrix with a minimum of 5 patients affected by any type of lesion in the
#' # gene to be included in the final matrix
#' lsn.type.mtx=prep.lsn.type.matrix(gene.lsn.overlap, min.ngrp=5)

prep.lsn.type.matrix=function(ov.data,    # output of the find.gene.lsn.overlaps function
                              min.ngrp=0) # if specified, genes with number of patients affected by all different types of lesions that's less than the specified number will be discarded (default is 0; will return all genes affected by any type of lesions in at least one patient).
{
  gene.lsn=ov.data$gene.lsn.hits

  # order and index gene-lesion overlaps by subject ID and gene
  ord=order(gene.lsn$ID,
            gene.lsn$gene)
  gene.lsn=gene.lsn[ord,]
  m=nrow(gene.lsn)
  new.block=which((gene.lsn$ID[-1]!=gene.lsn$ID[-m])|
                    (gene.lsn$gene[-1]!=gene.lsn$gene[-m]))
  row.start=c(1,new.block+1)
  row.end=c(new.block,m)
  ID.gene.index=cbind.data.frame(ID=gene.lsn$ID[row.start],
                                 gene=gene.lsn$gene[row.start],
                                 row.start=row.start,
                                 row.end=row.end)

  # initialize the result matrix
  lsn.types=sort(unique(gene.lsn$lsn.type))
  uniq.genes=unique(ID.gene.index$gene)
  uniq.IDs=unique(ID.gene.index$ID)

  grp.mtx=matrix("none",
                 length(uniq.genes),
                 length(uniq.IDs))
  rownames(grp.mtx)=uniq.genes
  colnames(grp.mtx)=uniq.IDs

  # fill in the matrix
  n.index=nrow(ID.gene.index)
  for (i in 1:n.index)
  {
    gene.lsn.rows=(ID.gene.index$row.start[i]:ID.gene.index$row.end[i])
    block.ID=ID.gene.index$ID[i]
    block.gene=ID.gene.index$gene[i]
    block.lsns=gene.lsn$lsn.type[gene.lsn.rows]
    block.lsns=unique(block.lsns)
    if (length(block.lsns)==1) grp.mtx[block.gene,block.ID]=block.lsns
    else grp.mtx[block.gene,block.ID]="multiple"
  }

  n.none=rowSums(grp.mtx== "none")
  lsn.grp.keep=(ncol(grp.mtx)-n.none)>=min.ngrp
  grp.mtx=grp.mtx[lsn.grp.keep,]

  return(grp.mtx)
}
