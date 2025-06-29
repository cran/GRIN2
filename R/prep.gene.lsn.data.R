
#' Prepare Gene and Lesion Data for GRIN Analysis
#'
#' @description
#' Prepares and indexes gene and lesion data for downstream GRIN (Genomic Random Interval) analysis. This function merges and orders gene and lesion coordinates to support efficient computation of overlaps between genes and all different types of genomic lesions (structural or sequence lesions).
#'
#' @param lsn.data A `data.frame` containing lesion data in GRIN-compatible format. Must include the following five columns:
#' \describe{
#'   \item{ID}{Unique patient identifier.}
#'   \item{chrom}{Chromosome on which the lesion is located.}
#'   \item{loc.start}{Start position of the lesion in base pairs.}
#'   \item{loc.end}{End position of the lesion in base pairs.}
#'   \item{lsn.type}{Type of lesion (e.g., gain, loss, mutation, fusion, etc...).}
#' }
#' @param gene.data A `data.frame` containing gene annotation data with the following four required columns:
#' \describe{
#'   \item{gene}{Ensembl gene ID.}
#'   \item{chrom}{Chromosome on which the gene is located.}
#'   \item{loc.start}{Start position of the gene in base pairs.}
#'   \item{loc.end}{End position of the gene in base pairs.}
#' }
#' @param mess.freq Integer specifying the frequency at which progress messages are displayed. Messages are printed every \code{mess.freq}-th lesion block processed (default is 10).
#'
#' @details
#' This function performs pre-processing by ordering and indexing both gene and lesion data. It combines gene and lesion coordinates into a unified structure, marking each with a specific code (`cty`) that identifies whether the row represents a gene or lesion. This merged data is then used in the `find.gene.lsn.overlaps()` function to detect gene-lesion overlaps.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{lsn.data}{Original lesion data.}
#'   \item{gene.data}{Original gene annotation data.}
#'   \item{gene.lsn.data}{Combined and ordered data.frame of gene and lesion intervals. The `cty` column encodes position type: 1 = gene start, 2 = lesion start, 3 = lesion end, 4 = gene end.}
#'   \item{gene.index}{Index data.frame indicating the start and end rows for each chromosome within `gene.lsn.data` for genes.}
#'   \item{lsn.index}{Index data.frame indicating the start and end rows for each lesion (grouped by type, chromosome, and subject) within `gene.lsn.data`.}
#' }
#'
#' @export
#'
#' @references
#' Pounds, S., et al. (2013). A genomic random interval model for statistical analysis of genomic lesion data.
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}

#'
#' @seealso \code{\link{order.index.gene.data}}, \code{\link{order.index.lsn.data}}, \code{\link{find.gene.lsn.overlaps}}
#' @examples
#' data(lesion_data)
#' data(hg38_gene_annotation)
#'
#' # Prepare gene and lesion data for GRIN analysis:
#' prep.gene.lsn <- prep.gene.lsn.data(lesion_data, hg38_gene_annotation)

prep.gene.lsn.data=function(lsn.data,       # lesion data file with five columns: "ID" which is patient ID, "chrom" chromosome on which the lesion is located, "loc.start" lesion start position, "loc.end" lesion end position, "lsn.type" lesion category such as gain, loss, mutation, etc...
                            gene.data,      # gene annotation data with four required columns: "gene" has ensembl ID, "chrom" chromosome on which the gene is located, "loc.start" gene start position, "loc.end" which is the gene end position
                            mess.freq=10)   # message frequency: display message every mess.freq^{th} lesion block
{
  old_opt <- options(stringsAsFactors = FALSE)
  on.exit(options(old_opt), add = TRUE)

  # order lesion data by type, chromosome, and subject
  lsn.dset=order.index.lsn.data(lsn.data)
  lsn.data=lsn.dset$lsn.data
  lsn.index=lsn.dset$lsn.index

  # order and index gene locus data by chromosome and position
  gene.dset=order.index.gene.data(gene.data)
  gene.data=gene.dset$gene.data
  gene.index=gene.dset$gene.index

  # Extract some basic information
  g=nrow(gene.data) # number of genes
  l=nrow(lsn.data)  # number of lesions

  # Create gene position data
  message(paste0("Formatting gene position data for counting: ",date()))
  gene.pos.data=rbind.data.frame(cbind.data.frame(ID="",  # gene start data
                                                  lsn.type="",
                                                  lsn.row=NA,
                                                  gene=gene.data[,"gene"],
                                                  gene.row=gene.data[,"gene.row"],
                                                  chrom=gene.data[,"chrom"],
                                                  pos=gene.data[,"loc.start"],
                                                  cty=1),
                                 cbind.data.frame(ID="", # gene end data
                                                  lsn.type="",
                                                  lsn.row=NA,
                                                  gene=gene.data[,"gene"],
                                                  gene.row=gene.data[,"gene.row"],
                                                  chrom=gene.data[,"chrom"],
                                                  pos=gene.data[,"loc.end"],
                                                  cty=4))
  # order gene position data
  ord=order(gene.pos.data[,"chrom"],
            gene.pos.data[,"pos"],
            gene.pos.data[,"cty"])
  gene.pos.data=gene.pos.data[ord,]

  # Create lesion position data with one row for each edge of each lesion
  message(paste0("Formatting lesion position data for counting:  ",date()))
  lsn.pos.data=rbind.data.frame(cbind.data.frame(ID=lsn.data[,"ID"],
                                                 lsn.type=lsn.data[,"lsn.type"],
                                                 lsn.row=lsn.data[,"lsn.row"],
                                                 gene="",
                                                 gene.row=NA,
                                                 chrom=lsn.data[,"chrom"],
                                                 pos=lsn.data[,"loc.start"],
                                                 cty=2),
                                cbind.data.frame(ID=lsn.data[,"ID"],
                                                 lsn.type=lsn.data[,"lsn.type"],
                                                 lsn.row=lsn.data[,"lsn.row"],
                                                 gene="",
                                                 gene.row=NA,
                                                 chrom=lsn.data[,"chrom"],
                                                 pos=lsn.data[,"loc.end"],
                                                 cty=3))
  # order lesion position data
  ord=order(lsn.pos.data[,"chrom"],
            lsn.pos.data[,"pos"],
            lsn.pos.data[,"cty"])
  lsn.pos.data=lsn.pos.data[ord,]

  # Combine gene & lesion data
  message(paste0("Combining formatted gene and lesion postion data: ",date()))
  gene.lsn.data=rbind.data.frame(gene.pos.data,
                                 lsn.pos.data)

  # Order and index gene & lesion data
  ord=order(gene.lsn.data[,"chrom"],
            gene.lsn.data[,"pos"],
            gene.lsn.data[,"cty"])
  gene.lsn.data=gene.lsn.data[ord,]
  m=nrow(gene.lsn.data)
  gene.lsn.data[,"glp.row"]=1:m

  # compute vector to order gene.lsn.data by lsn.row and gene.row
  ord=order(gene.lsn.data[,"lsn.row"],
            gene.lsn.data[,"gene.row"],
            gene.lsn.data[,"cty"])

  # use that vector to add gene.lsn.data row.start and row.end indices to lsn.data
  lsn.pos=gene.lsn.data[ord[1:(2*l)],]
  lsn.data[,"glp.row.start"]=lsn.pos[2*(1:l)-1,"glp.row"]
  lsn.data[,"glp.row.end"]=lsn.pos[2*(1:l),"glp.row"]


  # use that vector to add gene.lsn.data row.start and row.end indices to gene.data
  gene.pos=gene.lsn.data[ord[-(1:(2*l))],]
  gene.data[,"glp.row.start"]=gene.pos[2*(1:g)-1,"glp.row"]
  gene.data[,"glp.row.end"]=gene.pos[2*(1:g),"glp.row"]

  # Double-check table pointers from lsn.data and gene.data to gene.lsn.data

  message(paste0("Verifying structure of combined gene and lesion data: ",date()))
  glp.gene.start=gene.lsn.data[gene.data$glp.row.start,c("gene","chrom","pos")]
  colnames(glp.gene.start)=c("gene","chrom","loc.start")
  ok.glp.gene.start=all(glp.gene.start==gene.data[,c("gene","chrom","loc.start")])

  glp.gene.end=gene.lsn.data[gene.data$glp.row.end,c("gene","chrom","pos")]
  colnames(glp.gene.end)=c("gene","chrom","loc.end")
  ok.glp.gene.end=all(glp.gene.end==gene.data[,c("gene","chrom","loc.end")])

  glp.lsn.start=gene.lsn.data[lsn.data$glp.row.start,c("ID","chrom","pos","lsn.type")]
  colnames(glp.lsn.start)=c("ID","chrom","loc.start","lsn.type")
  ok.glp.lsn.start=all(glp.lsn.start==lsn.data[,c("ID","chrom","loc.start","lsn.type")])

  glp.lsn.end=gene.lsn.data[lsn.data$glp.row.end,c("ID","chrom","pos","lsn.type")]
  colnames(glp.lsn.end)=c("ID","chrom","loc.end","lsn.type")
  ok.glp.lsn.end=all(glp.lsn.end==lsn.data[,c("ID","chrom","loc.end","lsn.type")])

  # Double-check table pointers from gene.lsn.data to gene.data and lsn.data
  glp.gene.start=gene.lsn.data[gene.lsn.data$cty==1,c("gene.row","gene","chrom","pos")]
  ok.gene.start=all(glp.gene.start[,c("gene","chrom","pos")]==gene.data[glp.gene.start$gene.row,c("gene","chrom","loc.start")])

  glp.gene.end=gene.lsn.data[gene.lsn.data$cty==4,c("gene.row","gene","chrom","pos")]
  ok.gene.end=all(glp.gene.end[,c("gene","chrom","pos")]==gene.data[glp.gene.end$gene.row,c("gene","chrom","loc.end")])

  glp.lsn.start=gene.lsn.data[gene.lsn.data$cty==2,c("lsn.row","ID","chrom","pos","lsn.type")]
  ok.lsn.start=all(glp.lsn.start[,c("ID","chrom","pos","lsn.type")]==lsn.data[glp.lsn.start$lsn.row,c("ID","chrom","loc.start","lsn.type")])

  glp.lsn.end=gene.lsn.data[gene.lsn.data$cty==3,c("lsn.row","ID","chrom","pos","lsn.type")]
  ok.lsn.end=all(glp.lsn.end[,c("ID","chrom","pos","lsn.type")]==lsn.data[glp.lsn.end$lsn.row,c("ID","chrom","loc.end","lsn.type")])

  all.ok=all(c(ok.glp.gene.start,ok.glp.gene.end,
               ok.glp.lsn.start,ok.glp.lsn.end,
               ok.gene.start,ok.gene.end,
               ok.lsn.start,ok.lsn.end))

  if (!all.ok)
    stop("Error in constructing and indexing combined lesion and gene data.")

  message(paste0("Verified correct construction and indexing of combined lesion and gene data: ",date()))

  return(list(lsn.data=lsn.data, # Input lesion data
              gene.data=gene.data, # Input gene annotation data
              gene.lsn.data=gene.lsn.data, # data.frame ordered by gene and lesions start position. Gene start position is coded as 1 in the cty column and gene end position is coded as 4. Lesion start position is coded as 2 in the cty column and lesion end position is coded as 3.
              gene.index=gene.index, # data.frame that shows ordered row start and row end for each chromosome in the gene.lsn.data table
              lsn.index=lsn.index)) # data.frame that shows row start and row end for each lesion in the gene.lsn.data table

}
