
#' Count Gene Lesion Hits
#'
#' @description
#' The function computes the number of hits affecting each gene by lesion category. It also compute the number of subjects with a hit in each annotated gene by lesion category as well.
#'
#' @param ov.data a list of six data.frames that represent the output results of the find.gene.lsn.overlaps function.
#'
#' @details
#' The function use the output of the find.gene.lsn.overlaps function and return the number of unique subjects affected by each lesion category in the provided list of annotated genes and regulatory features (nsubj stats). It also count the number of hits affecting each loci per lesion category (nhits stats). For example, if NOTCH1 gene was found affected by three different mutations in the same subject, this patient will be considered as one subject in the nsubj stats but in the nhits stats for this event will be counted as 3 mutations that affect NOTCH1 gene.
#'
#' @return
#' A list with the following components:
#' \item{lsn.data}{Input lesion data}
#' \item{lsn.index}{data.frame that shows the overlapped gene-lesion data rows taht belong to each lesion in the gene.lsn.data table.}
#' \item{gene.data}{Input gene annotation data}
#' \item{gene.index}{data.frame with overlapped gene-lesion data rows that belong to each chromosome in the gene.lsn.data table.}
#' \item{nhit.mtx}{A data.frame with number of hits in each gene by lesion type (number of columns will be equal to the number of lesion types in the lsn.type column).}
#' \item{nsubj.mtx}{A data.frame with number of affected subjects by lesion type in each annotated gene.}
#' \item{gene.lsn.data}{Each row represent a gene overlapped by a certain lesion. Column "gene" shows the overlapped gene and "ID" column has the patient ID.}
#' \item{glp.data}{data.frame ordered by gene and lesions start position. Gene start position is coded as 1 in the cty column and gene end position is coded as 4. Lesion start position is coded as 2 in the cty column and lesion end position is coded as 3.}
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
#' # count number of subjects affected by different types of lesions and number of hits that affect
#' # each locus:
#' count.nsubj.nhits=count.hits(gene.lsn.overlap)

count.hits=function(ov.data) # output results of find.gene.lsn.overlaps function
{
  lsn.data=ov.data$lsn.data
  lsn.index=ov.data$lsn.index
  gene.lsn.hits=ov.data$gene.lsn.hits
  gene.lsn.data=ov.data$gene.lsn.data
  gene.data=ov.data$gene.data
  gene.index=ov.data$gene.index

  g=nrow(gene.data)

  # Compute the number of hits matrix
  lsn.types=sort(unique(lsn.index[,"lsn.type"]))
  k=length(lsn.types)

  nhit.mtx=matrix(0,g,k)
  colnames(nhit.mtx)=lsn.types

  nhit.tbl=table(gene.lsn.hits$gene.row,
                 gene.lsn.hits$lsn.type)
  nhit.rows=as.numeric(rownames(nhit.tbl))

  for (i in 1:ncol(nhit.tbl))
    nhit.mtx[nhit.rows,colnames(nhit.tbl)[i]]=nhit.tbl[,i]

  # Compute the matrix of the number of subjects with a hit
  gene.subj.type=paste0(gene.lsn.hits$gene.row,"_",
                        gene.lsn.hits$ID,"_",
                        gene.lsn.hits$lsn.type)
  dup.gene.subj.type=duplicated(gene.subj.type)

  subj.gene.hits=gene.lsn.hits[!dup.gene.subj.type,]

  nsubj.mtx=matrix(0,g,k)
  colnames(nsubj.mtx)=lsn.types

  nsubj.tbl=table(subj.gene.hits$gene.row,
                  subj.gene.hits$lsn.type)

  nsubj.rows=as.numeric(rownames(nsubj.tbl))

  for (i in 1:ncol(nsubj.tbl))
    nsubj.mtx[nsubj.rows,colnames(nsubj.tbl)[i]]=nsubj.tbl[,i]

  res=list(lsn.data=lsn.data,  # Input lesion data
           lsn.index=lsn.index, # data.frame that shows row start and row end for each lesion in the gene.lsn.data table
           gene.data=gene.data, # Input gene annotation data
           gene.index=gene.index, # data.frame that shows ordered row start and row end for each chromosome in the gene.lsn.data table
           nhit.mtx=nhit.mtx, # A data matrix with number of hits in each gene by lesion type
           nsubj.mtx=nsubj.mtx, # A data matrix with number of affected subjects by lesion type
           gene.lsn.data=gene.lsn.hits, # Each row represent a gene overlapped by a certain lesion. Column "gene" shows the overlapped gene and ID column has the patient ID
           glp.data=gene.lsn.data) # data.frame ordered by gene and lesions start position. Gene start position is coded as 1 in the cty column and gene end position is coded as 4. Lesion start position is coded as 2 in the cty column and lesion end position is coded as 3

  return(res)

}
