
#' Prepare Lesion Type Matrix
#'
#' @description
#' Constructs a matrix that summarizes the type(s) of lesions affecting each gene across patients. Each row represents a gene, and each column represents a patient.
#'
#' @param ov.data A list of six \code{data.frame} objects returned by the \code{\link{find.gene.lsn.overlaps}} function, containing gene-lesion overlap results.
#' @param min.ngrp Optional integer specifying the minimum number of patients that must be affected by any lesion in a given gene for that gene to be retained in the final matrix. The default is \code{0}, which includes all genes affected by any lesion in at least one patient.
#'
#' @details
#' This function produces a matrix with genes as rows and patients as columns. For each gene-patient pair:
#' \itemize{
#'   \item If the patient has no lesion in the gene, the entry is \code{"none"}.
#'   \item If the gene is affected by exactly one type of lesion in the patient (e.g., gain OR mutation), the entry is labeled with the corresponding lesion type.
#'   \item If the gene is affected by more than one lesion type in the patient (e.g., gain AND mutation), the entry is labeled as \code{"multiple"}.
#' }
#'
#' Genes affected in fewer than \code{min.ngrp} patients across all lesion types will be excluded if \code{min.ngrp > 0}.
#'
#' @return
#' A character matrix where:
#' \itemize{
#'   \item Rows correspond to genes (identified by Ensembl gene IDs).
#'   \item Columns correspond to patient IDs.
#'   \item Entries are \code{"none"}, a specific lesion type (e.g., \code{"gain"}, \code{"mutation"}), or \code{"multiple"}.
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
#' # 1) Prepare gene and lesion data:
#' prep.gene.lsn <- prep.gene.lsn.data(lesion_data, hg38_gene_annotation)
#'
#' # 2) Identify gene-lesion overlaps:
#' gene.lsn.overlap <- find.gene.lsn.overlaps(prep.gene.lsn)
#'
#' # 3) Create lesion type matrix for genes affected in >= 5 patients:
#' lsn.type.mtx <- prep.lsn.type.matrix(gene.lsn.overlap, min.ngrp = 5)

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
