#' Generate Waterfall Plots for Top Significant Genes
#'
#' @description
#' Generates waterfall plots for genes with significant associations between lesion status and expression level, based on the Kruskal Wallis (KW) test results. Only genes with q-values below a user-specified threshold will be plotted.
#'
#' @param out.dir A character string specifying the output directory where the waterfall plots for selected genes will be saved. Directory must exist or be created by the user prior to running the function.
#' @param alex.data A list of three data tables returned by \code{\link{alex.prep.lsn.expr}}:
#' \itemize{
#'   \item \code{"row.mtch"} matching expression and lesion rows by Ensembl gene ID.
#'   \item \code{"alex.expr"} matrix of gene expression.
#'   \item \code{"alex.lsn"} matrix of lesion status.
#'   }
#' All matrices must have rows ordered by Ensembl gene ID and columns ordered by patient ID.
#'
#' @param alex.kw.results A data table of Kruskal Wallis test results, returned by the \code{\link{KW.hit.express}} function.
#' @param q A numeric threshold indicating the maximum allowed KW q-value for a gene to be included in the waterfall plots.
#' @param lsn.data Lesion data provided in GRIN-compatible format (as used in \code{\link{alex.prep.lsn.expr}}).
#'
#' @details
#' For each gene in the \code{alex.kw.results} table with a q-value less than or equal to the user-specified \code{q} threshold, the function generates a waterfall plot displaying the relationship between lesion status and gene expression level. Each plot is saved as a separate PDF file in the \code{out.dir} folder.
#'
#' Internally, this function relies on helper functions such as \code{\link{alex.waterfall.prep}} and \code{\link{alex.waterfall.plot}} to prepare and render the plots.
#'
#' @return
#' No object is returned to the R environment. The function creates a set of PDF files (one per gene) in the specified \code{out.dir} directory. Each file contains a labeled waterfall plot illustrating gene expression across different lesion groups.
#'
#' @export
#'
#' @importFrom grDevices pdf dev.off
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link{alex.prep.lsn.expr}}, \code{\link{KW.hit.express}}, \code{\link{alex.waterfall.prep}}, \code{\link{alex.waterfall.plot}}
#'
#' @examples
#' data(expr_data)
#' data(lesion_data)
#' data(hg38_gene_annotation)
#'
#' # 1) Prepare matched expression and lesion matrices:
#' alex.data <- alex.prep.lsn.expr(expr_data, lesion_data,
#'                                 hg38_gene_annotation,
#'                                 min.expr = 5, min.pts.lsn = 5)
#' # 2) Run Kruskal Wallis test:
#' alex.kw.results <- KW.hit.express(alex.data,
#'                                   hg38_gene_annotation,
#'                                   min.grp.size = 5)
#'
#' # 3) Create temporary output folder and generate waterfall plots:
#' dir.create(resultsFolder <- file.path(tempdir(), "temp.out"))
#' waterfall.plts <- top.alex.waterfall.plots(out.dir = resultsFolder,
#'                                            alex.data = alex.data,
#'                                            alex.kw.results = alex.kw.results,
#'                                            q = 1e-15,
#'                                            lsn.data = lesion_data)
#'
#' # Clean up:
#' unlink(resultsFolder, recursive = TRUE)

top.alex.waterfall.plots=function(out.dir,            # path to the folder to add waterfall plots
                                  alex.data,          # output of the alex.prep.lsn.expr function (list of three data table "alex.expr" with expression data, "alex.lsn" with lesion data and row.mtch)
                                  alex.kw.results,    # ALEX Kruskal-Wallis test results (output of the KW.hit.express function)
                                  q,                  # minimum q value for a gene to be plotted
                                  lsn.data)           # Lesion data in a GRIN compatible format


{
  # To extract genes below specified q.value for waterfall plots:
  KW.results=alex.kw.results[!is.na(alex.kw.results$q.KW),]
  selected.genes=KW.results[KW.results$q.KW<q,]
  selected.genes=selected.genes[!is.na(selected.genes$gene.name),]
  selected.genes=selected.genes[!duplicated(selected.genes[ , "gene.name"]),]
  top.genes=selected.genes$gene.name

  for (i in 1:length(top.genes))
  {
    grDevices::pdf(paste0(out.dir,top.genes[i], "_waterfall_plot", ".pdf"),
        height=8,width=10)
    temp=alex.waterfall.prep(alex.data, alex.kw.results,top.genes[i],lsn.data)
    plot.gene=alex.waterfall.plot(temp, lsn.data)
    grDevices::dev.off()
  }
}
