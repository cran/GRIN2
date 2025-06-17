
#' Order and Index Gene Annotation Data
#'
#' @description
#' This function orders and indexes gene annotation data by chromosome, gene start, and gene end positions. It is typically used to prepare gene data for overlap analyses with lesion data.
#'
#' @param gene.data A `data.frame` containing gene annotation information, either provided by the user or retrieved using the `get.ensembl.annotation` function from the GRIN2.0 package. The `data.frame` must contain four columns:
#' \describe{
#'   \item{"gene"}{Ensembl gene ID.}
#'   \item{"chrom"}{Chromosome on which the gene is located.}
#'   \item{"loc.start"}{Start position of the gene.}
#'   \item{"loc.end"}{End position of the gene.}
#' }
#'
#' @return
#' A list with two components:
#' \item{gene.data}{The input gene annotation data, ordered by chromosome and genomic coordinates.}
#' \item{gene.index}{A `data.frame` with two columns (`row.start` and `row.end`) indicating the row indices for genes on each chromosome.}
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
#' data(hg38_gene_annotation)
#'
#' ordered.genes <- order.index.gene.data(hg38_gene_annotation)

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
