
#' Get Ensembl Gene and Regulatory Features Annotation Data
#'
#' @description
#' Retrieves gene and regulatory feature annotation data from the Ensembl BioMart database
#' for the Human GRCh38 (hg38) genome assembly using the biomaRt package.
#'
#' @param genome.assembly Character string. Currently, only `"Human_GRCh38"` is supported.
#'
#' @details
#' This function retrieves:
#' \itemize{
#'   \item Gene annotation: Ensembl ID, chromosome, start/end positions, gene name, description, biotype, strand, and cytogenetic band.
#'   \item Predicted regulatory features: Promoters, enhancers, CTCF binding sites, etc., from the Ensembl regulatory build.
#'   \item Validated regulatory features: Experimentally confirmed enhancers and TSSs from the FANTOM5 project.
#' }
#'
#' @return A list with three components:
#' \describe{
#'   \item{gene.annotation}{Data frame of gene-level annotation.}
#'   \item{reg.annotation.predicted}{Data frame of predicted regulatory features.}
#'   \item{reg.annotation.validated}{Data frame of validated regulatory features (FANTOM5).}
#' }
#'
#' @export
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' Zerbino, Daniel R., et al. (2015). The ensembl regulatory build.
#'
#' Kinsella, Rhoda J., et al. (2011). Ensembl BioMarts: a hub for data retrieval across taxonomic space.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link[biomaRt]{useEnsembl}}, \code{\link[biomaRt]{getBM}}
#'
#' @examples
#' \donttest{
#' hg38.ann <- get.ensembl.annotation("Human_GRCh38")
#' # gene level annotations:
#' hg38.genes <- hg38.ann$gene.annotation
#' # regulatory sequences from the ensembl genome build:
#' hg38.reg.pred <- hg38.ann$reg.annotation.predicted
#' # regulatory sequences from the FANTOM5 project:
#' hg38.reg.val <- hg38.ann$reg.annotation.validated
#' }

get.ensembl.annotation <- function(genome.assembly) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    message("Package 'biomaRt' is not installed. Please install it using BiocManager::install('biomaRt')")
    return(NULL)
  }

  if (genome.assembly != "Human_GRCh38") {
    stop("Unsupported genome assembly. Only 'Human_GRCh38' is currently supported.")
  }

  res <- tryCatch({
    message(paste0("Retrieving gene annotation data from Ensembl BioMart: ", date()))

    ensembl_GRCh38 <- biomaRt::useEnsembl(
      biomart = "genes",
      dataset = "hsapiens_gene_ensembl",
      version = "110"
    )

    chromosomes <- c(1:22, "X", "Y")
    gene.attrs <- c(
      "ensembl_gene_id", "chromosome_name", "start_position", "end_position",
      "description", "external_gene_name", "gene_biotype", "strand", "band"
    )

    gene.data <- biomaRt::getBM(
      attributes = gene.attrs,
      filters = "chromosome_name",
      values = chromosomes,
      mart = ensembl_GRCh38
    )

    colnames(gene.data) <- c(
      "gene", "chrom", "loc.start", "loc.end", "description",
      "gene.name", "biotype", "chrom.strand", "chrom.band"
    )

    ## Predicted regulatory features
    message(paste0("Retrieving computationally predicted regulatory features: ", date()))
    reg.mart <- biomaRt::useEnsembl(
      biomart = "regulation",
      dataset = "hsapiens_regulatory_feature",
      version = "110"
    )

    reg.features <- c("Promoter", "Promoter Flanking Region", "Enhancer",
                      "CTCF Binding Site", "TF binding site", "Open chromatin")

    reg.data <- biomaRt::getBM(
      attributes = c("regulatory_stable_id", "chromosome_name", "chromosome_start",
                     "chromosome_end", "feature_type_description"),
      filters = c("regulatory_feature_type_name", "chromosome_name"),
      values = list(reg.features, chromosomes),
      mart = reg.mart
    )

    colnames(reg.data) <- c("gene", "chrom", "loc.start", "loc.end", "description")

    ## Validated regulatory features (FANTOM)
    message(paste0("Retrieving experimentally validated regulatory features (FANTOM): ", date()))
    val.mart <- biomaRt::useEnsembl(
      biomart = "regulation",
      dataset = "hsapiens_external_feature",
      version = "110"
    )

    val.data <- biomaRt::getBM(
      attributes = c("display_label", "chromosome_name", "chromosome_start",
                     "chromosome_end", "feature_type_description"),
      filters = c("external_feature_set_name", "chromosome_name"),
      values = list("FANTOM predictions", chromosomes),
      mart = val.mart
    )

    colnames(val.data) <- c("gene", "chrom", "loc.start", "loc.end", "description")

    ## Final result
    list(
      gene.annotation = gene.data,
      reg.annotation.predicted = reg.data,
      reg.annotation.validated = val.data
    )

  }, error = function(e) {
    message("Failed to connect to Ensembl: ", conditionMessage(e), "\n",
            "Please check http://status.ensembl.org or try setting a mirror (see ?useEnsembl).")
    return(NULL)
  })

  return(res)
}
