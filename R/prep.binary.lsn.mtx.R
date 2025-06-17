
#' Prepare Binary Lesion Matrix
#'
#' @description
#' Constructs a binary matrix representing the presence or absence of specific lesion types affecting individual genes across patients. Each row corresponds to a gene-lesion type combination, and each column corresponds to a patient.
#'
#' @param ov.data A list of six \code{data.frame} objects representing the output from the \code{\link{find.gene.lsn.overlaps}} function.
#' @param min.ngrp Optional integer specifying the minimum number of patients that must be affected by a given gene-lesion combination to be retained in the output matrix. The default is \code{0}, which includes all combinations affecting at least one patient.
#'
#' @details
#' The function processes the overlap results from \code{\link{find.gene.lsn.overlaps}} and constructs a binary matrix with dimensions: (gene and lesion type) by patient.
#'
#' Each row is labeled using the format \code{<gene.ID>_<lesion.type>} (e.g., \code{ENSG00000118513_gain} for a gain affecting the MYB gene). For each gene-lesion combination, a patient receives a value of \code{1} if affected by that specific lesion type in the corresponding gene, and \code{0} otherwise.
#'
#' Rows representing rare lesions (i.e., affecting fewer patients than \code{min.ngrp}) are excluded from the final matrix if \code{min.ngrp > 0}.
#'
#' @return
#' A binary matrix (as a \code{data.frame}) where:
#' \itemize{
#'   \item Rows correspond to gene-lesion combinations (\code{gene.ID_lesion.type}).
#'   \item Columns correspond to patient IDs.
#'   \item Entries are binary: \code{1} if the patient is affected by such a specific type of  lesion in that gene, \code{0} otherwise.
#' }
#'
#' @export
#'
#' @references
#' Pounds, S., et al. (2013). A genomic random interval model for statistical analysis of genomic lesion data.
#'
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link{prep.gene.lsn.data}}, \code{\link{find.gene.lsn.overlaps}}
#'
#' @examples
#' data(lesion_data)
#' data(hg38_gene_annotation)
#'
#' # 1) Prepare gene-lesion input data:
#' prep.gene.lsn <- prep.gene.lsn.data(lesion_data,
#'                                     hg38_gene_annotation)
#'
#' # 2) Identify gene-lesion overlaps:
#' gene.lsn.overlap <- find.gene.lsn.overlaps(prep.gene.lsn)
#'
#' # 3) Create binary lesion matrix including only lesion-gene pairs affecting >= 5 patients:
#' lsn.binary.mtx <- prep.binary.lsn.mtx(gene.lsn.overlap, min.ngrp = 5)

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
