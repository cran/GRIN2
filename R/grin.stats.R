
#' GRIN Statistics Output
#'
#' @description
#' The function run the Genomic Random Interval (GRIN) analysis to determine whether a certain locus has an abundance of lesions or a constellation of multiple types of lesions that is statistically significant.
#'
#' @param lsn.data data.frame with lesion data prepared by the user in a GRIN compatible format. Object should has five columns that include "ID" with patient ID, "chrom" which is the chromosome on which the lesion is located, "loc.start" which is the lesion start position, "loc.end" the lesion end position and "lsn.type" which is the lesion category for example gain, loss, mutation, fusion, etc... For Single Nucleotide Variants (SNVs), loc.start will be the same as loc.end. For Copy Number Alterations (CNAs) such as gain and deletions, loc.start and loc.end should be the gain or deletion start and end positions respectively. For structural rearrangements such as inversions and translocations, each rearrangement should be coded in two different lines, one line for chromosome A involved in the translocation break-point and the second line for chromosome B break-point. For inversions on the same chromosome, the two lines will include the two breakpoints of the inversion. An example lesion data in a GRIN compatible format can be found at the GRIN2.0 package data folder (lesion.data.rda).
#' @param gene.data data.frame with the gene annotation data either provided by the user or directly retreived from ensembl BioMart database using get.ensembl.annotation function included in the GRIN2.0 library if the genome.version is specified. Object should has four columns "gene" which is the ensembl ID of annotated genes to which the lesion data will be overlapped, "chrom" which is the chromosome on which the gene is located, "loc.start" which is the gene start position, and "loc.end" the gene end position.
#' @param chr.size data.frame with the size of the 22 autosomes, in addition to X and Y chromosomes in base pairs. It should has two columns that include "chrom" with the chromosome number and "size" for the size of the chromosome in base pairs. Chromosome size data can be either provided by the user or directly retreived from UCSC genome browser using get.chrom.length function included in the GRIN2.0 library if genome.version is specified.
#' @param genome.version Genome assembly should be only specified if the user selected not to provide gene annotation, chromosome size files, and directly retrieve those files from ensembl BioMart database, and UCSC genome browsers using get.ensembl.annotation and get.chrom.length functions respectively. Currently, the function support four genome assemblies that include "Human_GRCh38", "Human_GRCh37", "Mouse_HGCm39", and "Mouse_HGCm38".
#'
#' @details
#' The function run the Genomic Random Interval (GRIN) analysis and evaluates the probability of each gene locus to be affected by different types of lesions based on a convolution of independent but non-identical Bernoulli distributions to determine whether this locus has an abundance of lesions that is statistically significant.In addition, FDR-adjusted q value is computed for each locus based on Pounds & Cheng (2006) estimator of the proportion of tests with a true null (pi.hat). The function also evaluates if a certain locus is affected by a constellation of multiple types of lesions and return the GRIN results table.
#'
#' @return
#' A list with the following components:
#' \item{gene.hits}{data table of GRIN results that include gene annotation, number of subjects affected by each lesion type for example gain, loss, mutation, etc.., and number of hits affecting each locus. The GRIN results table will also include P and FDR adjusted q-values showing the probability of each locus of being affected by one or a constellation of multiple types of lesions.}
#' \item{lsn.data}{input lesion data}
#' \item{gene.data}{input gene annotation data}
#' \item{gene.lsn.data}{each row represent a gene overlapped by a certain lesion. Column "gene" shows the overlapped gene ensembl ID and "ID"" column has the patient ID.}
#' \item{chr.size}{data table showing the size of the 22 autosomes, in addition to X and Y chromosomes in base pairs.}
#' \item{gene.index}{data.frame with overlapped gene-lesion data rows that belong to each chromosome in the gene.lsn.data table.}
#' \item{lsn.index}{data.frame that shows the overlapped gene-lesion data rows taht belong to each lesion in the gene.lsn.data table.}
#'
#' @export
#'
#' @importFrom utils select.list
#'
#' @references
#' Pounds, Stan, et al. (2013) A genomic random interval model for statistical analysis of genomic lesion data.
#'
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author {Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [prep.gene.lsn.data()], [find.gene.lsn.overlaps()], [count.hits()], [prob.hits()]
#'
#' @examples
#' data(lesion.data)
#' data(hg19.gene.annotation)
#' data(hg19.chrom.size)
#'
#' # if gene annotation and chromosome size files will be provided by the user:
#' grin.results=grin.stats(lesion.data,
#'                         hg19.gene.annotation,
#'                         hg19.chrom.size)
#'
#' # to directly retrieve gene annotation and chromosome size files from Ensembl BioMart database,
#' # and UCSC genome browsers using get.ensembl.annotation and get.chrom.length functions respectively,
#' # users can select to specify certain genome assembly using the 'genome.version' argument:
#' # "Human_GRCh37" can be used for the GRCH37 (hg19) genome assembly, and "Human_GRCh38" can be used
#' # for the GRCH38 (hg38) genome assembly

grin.stats=function(lsn.data,            # data.frame with five columns "ID" which is the subject identifier), "chrom" which is the chromosome on which the lesion is located), "loc.start" the lesion start position, "loc.end" the lesion end position, and "lsn.type" which is the lesion category for example gain, mutation, etc..)
                    gene.data=NULL,      # data.frame with four required columns "gene" which is the gene ensembl ID, "chrom" which is the chromosome on which the gene is located, "loc.start" the gene start position, and "loc.end" which is the gene end position
                    chr.size=NULL,       # data.frame with two columns "chrom" which is the chromosome number and "size" which is the size of the chromosome in base pairs)
                    genome.version=NULL) # character string with genome version. Currently, four genome assemblies are supported including "Human_GRCh38", "Human_GRCh37", "Mouse_HGCm39", and "Mouse_HGCm38"

{
  if (is.null(genome.version)&&(is.null(gene.data)||is.null(chr.size)))
  {
    genome.version=utils::select.list(c("Human_GRCh38",
                                        "Human_GRCh37",
                                        "Mouse_HGCm39",
                                        "Mouse_HGCm38"))
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
