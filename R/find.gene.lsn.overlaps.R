
#' Find Gene Lesion Overlaps
#'
#' @description
#' Identifies overlaps between genes and genomic lesions using the output from the `prep.gene.lsn.data()` function.
#' This function detects all instances where a lesion spans or intersects the genomic coordinates of a gene.
#'
#' @param gl.data A list of five `data.frame` objects returned by the `prep.gene.lsn.data()` function.
#' These include processed and indexed gene and lesion data ready for overlap analysis.
#'
#' @return
#' A list containing the following components:
#' \item{lsn.data}{Original input lesion data.}
#' \item{gene.data}{Original input gene annotation data.}
#' \item{gene.lsn.data}{A `data.frame` ordered by chromosome and start position, containing both gene and lesion entries.
#' Each row includes a `cty` code indicating the type and boundary of the interval:
#' 1 = gene start, 2 = lesion start, 3 = lesion end, 4 = gene end.}
#' \item{gene.lsn.hits}{A `data.frame` where each row corresponds to a gene overlapped by a lesion.
#' Includes 11 columns:
#' \code{"gene"} (Ensembl gene ID),
#' \code{"gene.chrom"}, \code{"gene.loc.start"}, \code{"gene.loc.end"} (chromosome and coordinates of the gene),
#' \code{"ID"} (patient/sample ID),
#' \code{"lsn.chrom"}, \code{"lsn.loc.start"}, \code{"lsn.loc.end"} (chromosome and coordinates of the lesion),
#' and \code{"lsn.type"} (type of lesion).}
#' \item{gene.index}{A `data.frame` indexing the rows in `gene.lsn.data` that correspond to each chromosome's genes (row start and row end per chromosome).}
#' \item{lsn.index}{A `data.frame` indexing the rows in `gene.lsn.data` that correspond to each lesion (row start and row end per lesion).}
#'
#' @export
#'
#' @references
#' Pounds, S. et al. (2013). A genomic random interval model for statistical analysis of genomic lesion data.
#'
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link{prep.gene.lsn.data}}
#'
#' @examples
#' data(lesion_data)
#' data(hg38_gene_annotation)
#'
#' # Prepare gene and lesion data for GRIN-based computations:
#' prep.gene.lsn <- prep.gene.lsn.data(lesion_data, hg38_gene_annotation)
#'
#' # Identify genes that are overlapped by lesions:
#' gene.lsn.overlap <- find.gene.lsn.overlaps(prep.gene.lsn)

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
