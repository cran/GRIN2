
#' Execute GRIN Statistical Framework
#'
#' @description
#' Executes the Genomic Random Interval (GRIN) statistical framework to determine whether a specific genomic locus (gene or regulatory region) is significantly affected by either individual or a constellation of multiple lesion types.
#'
#' @param lsn.data A `data.frame` containing lesion data formatted for GRIN analysis. It must include the following five columns:
#'   \itemize{
#'     \item{\code{ID}}: Sample or patient identifier.
#'     \item{\code{chrom}}: Chromosome on which the lesion is located.
#'     \item{\code{loc.start}}: Genomic start coordinate of the lesion.
#'     \item{\code{loc.end}}: Genomic end coordinate of the lesion.
#'     \item{\code{lsn.type}}: Lesion type (e.g., gain, loss, mutation, fusion, etc...).
#'   }
#'   For Single Nucleotide Variants (SNVs), loc.start and loc.end should be the same. For Copy Number Alterations (CNAs) such as gains and deletions, these fields represent the lesion start and end positions (lesion boundary). Structural rearrangements (e.g., translocations, inversions) should be represented by two entries (two separate rows), one for each breakpoint. An example dataset is available in the GRIN2 package (`lesion_data.rda`).
#'
#' @param gene.data A `data.frame` containing gene annotation data. Must include the following columns:
#'   \itemize{
#'     \item{\code{gene}}: Ensembl gene ID.
#'     \item{\code{chrom}}: Chromosome where the gene is located.
#'     \item{\code{loc.start}}: Gene start position.
#'     \item{\code{loc.end}}: Gene end position.
#'   }
#'   This data can be user-provided or retrieved automatically via `get.ensembl.annotation()` if \code{genome.version} is specified.
#'
#' @param chr.size A `data.frame` specifying chromosome sizes. Must contain:
#'   \itemize{
#'     \item{\code{chrom}}: Chromosome number.
#'     \item{\code{size}}: Chromosome length in base pairs.
#'   }
#'   The data can be user-provided or directly retrieved using `get.chrom.length()` if \code{genome.version} is specified.
#'
#' @param genome.version Optional. If gene annotation and chromosome size files are not provided, users can specify a supported genome assembly to retrieve these files automatically. Currently, the package only support "Human_GRCh38" genome assembly.
#'
#' @details
#' The GRIN algorithm evaluates each locus to determine whether the observed frequency and distribution of lesions is greater than expected by chance. This is modeled using a convolution of independent, non-identical Bernoulli distributions, accounting for lesion type, locus size, and chromosome context.
#'
#' For each gene, the function calculates:
#' \itemize{
#'   \item{A p-value for the enrichment of lesion events}
#'   \item{An FDR-adjusted q-value using the Pounds & Cheng (2006) method}
#'   \item{Significance of multi-lesion constellation patterns (e.g., p-value for a locus being affected by 1, 2, etc., lesion types)}
#' }
#'
#' @return
#' A list containing:
#' \item{gene.hits}{A `data.frame` of GRIN results for each gene, including annotation, subject/hit counts by lesion type, and p/q-values for individual and multi-lesion constellation significance.}
#' \item{lsn.data}{The original lesion input data.}
#' \item{gene.data}{The original gene annotation input data.}
#' \item{gene.lsn.data}{A `data.frame` where each row represents a gene-lesion overlap. Includes columns \code{"gene"} (Ensembl ID) and \code{"ID"} (sample ID).}
#' \item{chr.size}{The chromosome size reference table used in computations.}
#' \item{gene.index}{Indexes linking genes to rows in \code{gene.lsn.data} by chromosome.}
#' \item{lsn.index}{Indexes linking lesions to rows in \code{gene.lsn.data}.}
#'
#' @export
#'
#' @importFrom utils select.list
#'
#' @references
#' Pounds, S. et al. (2013). A genomic random interval model for statistical analysis of genomic lesion data.
#'
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link{prep.gene.lsn.data}}, \code{\link{find.gene.lsn.overlaps}}, \code{\link{count.hits}}, \code{\link{prob.hits}}
#'
#' @examples
#' data(lesion_data)
#' data(hg38_gene_annotation)
#' data(hg38_chrom_size)
#'
#' # Example1: Run GRIN with user-supplied annotation and chromosome size:
#' grin.results <- grin.stats(lesion_data,
#'                            hg38_gene_annotation,
#'                            hg38_chrom_size)
#'
#' # Example 2: User can specify genome version to automatically retrieve annotation
#' # and chromosome size data:
#' # grin.results <- grin.stats(lesion_data,
#' #                            genome.version = "Human_GRCh38")

grin.stats=function(lsn.data,            # data.frame with five columns "ID" which is the subject identifier), "chrom" which is the chromosome on which the lesion is located), "loc.start" the lesion start position, "loc.end" the lesion end position, and "lsn.type" which is the lesion category for example gain, mutation, etc..)
                    gene.data=NULL,      # data.frame with four required columns "gene" which is the gene ensembl ID, "chrom" which is the chromosome on which the gene is located, "loc.start" the gene start position, and "loc.end" which is the gene end position
                    chr.size=NULL,       # data.frame with two columns "chrom" which is the chromosome number and "size" which is the size of the chromosome in base pairs)
                    genome.version=NULL) # character string with genome version. Currently, four genome assemblies are supported including "Human_GRCh38", "Human_GRCh37", "Mouse_HGCm39", and "Mouse_HGCm38"

{
  if (is.null(genome.version)&&(is.null(gene.data)||is.null(chr.size)))
  {
    genome.version=utils::select.list(c("Human_GRCh38"))
  }

  if (is.character(genome.version))
  {
    ensembl.data=get.ensembl.annotation(genome.version)
    if (is.null(gene.data))
    {
      gene.data=ensembl.data$gene.annotation
    }
    chrom.size=get.chrom.length(genome.version)
    if (is.null(chr.size))
    {
      chr.size=chrom.size
    }
  }

  prep.data=prep.gene.lsn.data(lsn.data,
                               gene.data)
  find.overlap=find.gene.lsn.overlaps(prep.data)
  hit.cnt=count.hits(find.overlap)
  hit.pvals=prob.hits(hit.cnt,
                      chr.size)
  return(hit.pvals)
}
