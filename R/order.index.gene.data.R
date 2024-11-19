
#' Order Index Gene Data
#'
#' @description
#' This function order and index gene annotation data by chromosome on which the gene is located, gene start, and end positions.
#'
#' @param gene.data data.frame with gene annotation data either provided by the user or retrieved from ensembl BioMart database using get.ensembl.annotation function included in the GRIN2.0 library. data.frame should has four columns that include "gene" which is the ensembl ID of the annotated genes to which the lesion data will be overlapped, "chrom" which is the chromosome on which the gene is located, "loc.start" which is the gene start position, and "loc.end" the gene end position.
#'
#' @return
#' A list with the following components:
#' \item{gene.data}{Input gene annotation data}
#' \item{gene.index}{data.frame with two columns of ordered row start and row end based on the number of genes annotated to each chromosome.}
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
#' data(hg19.gene.annotation)
#'
#' ordered.genes=order.index.gene.data(hg19.gene.annotation)

order.index.gene.data=function(gene.data)  # gene annotation data
{
  g=nrow(gene.data)
  gene.ord=order(gene.data[,"chrom"],
                 gene.data[,"loc.start"],
                 gene.data[,"loc.end"])

  gene.data=gene.data[gene.ord,]
  new.chrom=which(gene.data[-1,"chrom"]!=gene.data[-g,"chrom"])
  chr.start=c(1,new.chrom+1)
  chr.end=c(new.chrom,g)
  gene.index=cbind.data.frame(chrom=gene.data[chr.start,"chrom"],
                              row.start=chr.start,
                              row.end=chr.end)
  gene.data$gene.row=1:g

  res=list(gene.data=gene.data,
           gene.index=gene.index)
  return(res)
}
