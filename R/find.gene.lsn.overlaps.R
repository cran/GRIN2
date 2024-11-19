
#' Find Gene Lesion Overlaps
#'
#' @description
#' The function use the output of the prep.gene.lsn.data function to find lesion-gene overlaps.
#'
#' @param gl.data a list of five data.frames that represent the output results of the prep.gene.lsn.data function.
#'
#' @return
#' A list with the following components:
#' \item{lsn.data}{Input lesion data}
#' \item{gene.data}{Input gene annotation data}
#' \item{gene.lsn.data}{data.frame ordered by gene and lesions start position. Gene start position is coded as 1 in the cty column and gene end position is coded as 4. Lesion start position is coded as 2 in the cty column and lesion end position is coded as 3.}
#' \item{gene.lsn.hits}{data.frame on which each row represent a gene overlapped by a certain lesion. The data.frame has 11 columns that include "gene" with ensembl ID of the overlapped gene, "gene.chrom", "gene.loc.start" and "gene.loc.end" with data for the chromosome on which the gene is located, start and end positions of the gene. In addition, column "ID" has the ID of the patient with a lesion that overlapped this gene, "lsn.chrom", "lsn.loc.start", "lsn.loc.end" and "lsn.type" have data for the chromosome, lesion start, lesion end positions and the lesion type respectively.}
#' \item{gene.index}{data.frame that shows row start and row end for each chromosome in the gene.lsn.data table}
#' \item{lsn.index}{data.frame that shows row start and row end for each lesion in the gene.lsn.data table}
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
#' @seealso [prep.gene.lsn.data()]
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

find.gene.lsn.overlaps=function(gl.data) # A list of five data.frames that represent the output results of the prep.gene.lsn.data function

{
  gene.data=gl.data$gene.data
  lsn.data=gl.data$lsn.data
  lsn.index=gl.data$lsn.index
  gene.index=gl.data$gene.index
  gene.lsn.data=gl.data$gene.lsn.data
  m=nrow(gene.lsn.data)

  message(paste0("Scanning through combined lesion and gene data to find gene-lesion overlaps: ",date()))


  gene.row.mtch=NULL  # initialize vector for rows of gene data matched to rows of lesion data
  lsn.row.mtch=NULL   # initialize vector for rows of lesion data matched to rows of gene data
  current.genes=NULL  # initialize vector of genes overlapping this point of the scan
  current.lsns=NULL   # initialize vector of lesions overlapping this point of the scan
  for (i in 1:m)      # loop over rows of gene.lsn.data
  {
    # enter a gene
    if (gene.lsn.data$cty[i]==1)
    {
      # add this gene to the set of current genes
      current.genes=c(current.genes,
                      gene.lsn.data$gene.row[i])

      # match this gene to set of current lesions
      lsn.row.mtch=c(lsn.row.mtch,current.lsns)          # add current lesions to lsn.row.mtch
      gene.row.mtch=c(gene.row.mtch,                     # add this gene for each current lesion
                      rep(gene.lsn.data$gene.row[i],
                          length(current.lsns)))

    }

    # exit a gene
    if (gene.lsn.data$cty[i]==4)
    {
      # drop this gene from the set of current genes
      current.genes=setdiff(current.genes,
                            gene.lsn.data$gene.row[i])
    }


    if (gene.lsn.data$cty[i]==2)       # enter a lesion
    {
      lsn.row.mtch=c(lsn.row.mtch,
                     rep(gene.lsn.data$lsn.row[i],length(current.genes)))
      gene.row.mtch=c(gene.row.mtch,current.genes)
      current.lsns=c(current.lsns,
                     gene.lsn.data$lsn.row[i])
    }

    if (gene.lsn.data$cty[i]==3)
    {
      current.lsns=setdiff(current.lsns,
                           gene.lsn.data$lsn.row[i])
    }
  }

  message(paste0("Completed scan of combined gene-lesion data: ",date()))

  # Generate the gene-lesion hit data
  gene.lsn.hits=cbind.data.frame(gene.data[gene.row.mtch,
                                           c("gene.row","gene","chrom","loc.start","loc.end")],
                                 lsn.data[lsn.row.mtch,
                                          c("lsn.row","ID","chrom","loc.start","loc.end","lsn.type")])

  colnames(gene.lsn.hits)=c("gene.row","gene","gene.chrom","gene.loc.start","gene.loc.end",
                            "lsn.row","ID","lsn.chrom","lsn.loc.start","lsn.loc.end","lsn.type")

  res=list(lsn.data=lsn.data, # Input lesion data
           gene.data=gene.data, # Input gene annotation data
           gene.lsn.data=gene.lsn.data, # Data.frame ordered by gene and lesions start position. Gene start position is coded as 1 in the cty column and gene end position is coded as 4. Lesion start position is coded as 2 in the cty column and lesion end position is coded as 3
           gene.lsn.hits=gene.lsn.hits, # Each row represent a gene overlapped by a certain lesion. Gene column shows the overlapped gene and ID column has the patient ID
           gene.index=gene.index, # Data.frame that shows row start and row end for each chromosome in the gene.lsn.data table
           lsn.index=lsn.index) # Data.frame that shows row start and row end for each lesion in the gene.lsn.data table

  return(res)

}
