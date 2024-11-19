
#' Waterfall Plots for Lesion and Expression Data of Top Significant Genes
#'
#' @description
#' Function return waterfall plots for top significant genes in the KW results table based on the specified q value.
#'
#' @param out.dir Path to the folder where the waterfall plots of selected genes based on the specified q value of the KW results table will be added.
#' @param alex.data output of the alex.prep.lsn.expr function. It's a list of three data tables that include "row.mtch", "alex.expr" with expression data, "alex.lsn" with lesion data. Rows of alex.expr, and "alex.lsn" matrices are ordered by gene ensembl IDs and columns are ordered by patient ID.
#' @param alex.kw.results ALEX Kruskal-Wallis test results (output of the KW.hit.express function).
#' @param q Maximum allowed KW q-value threshold for a gene to be plotted based on the output of the KW.hit.express function.
#' @param lsn.data Lesion data in a GRIN compatible format.
#'
#' @details
#' Function will return waterfall plots for top significant genes in the KW results table based on the user specified q-value threshold of the KW test.The plots will be added to the user specified outdir folder.
#'
#' @return
#' Function will return waterfall plots for top significant genes in the KW results table.
#'
#' @export
#'
#' @importFrom grDevices dev.off pdf
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author {Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [alex.prep.lsn.expr()], [KW.hit.express()], [alex.waterfall.prep()], [alex.waterfall.plot()]
#'
#' @examples
#' data(expr.data)
#' data(lesion.data)
#' data(hg19.gene.annotation)
#'
#' # prepare expression, lesion data and return the set of genes with both types of data available
#' # ordered by gene IDs in rows and patient IDs in columns:
#' alex.data=alex.prep.lsn.expr(expr.data, lesion.data,
#'                              hg19.gene.annotation, min.expr=5,
#'                               min.pts.lsn=5)
#'
#' # run KW test for association between lesion groups and expression level of the same gene:
#' alex.kw.results=KW.hit.express(alex.data, hg19.gene.annotation, min.grp.size=5)
#'
#' # return waterfall plots for a list of top significant genes to a pre-specified folder:
#' dir.create(resultsFolder <- file.path(tempdir(), "temp.out"))
#'
#' waterfall.plts=top.alex.waterfall.plots(out.dir=resultsFolder,
#'                                         alex.data, alex.kw.results,
#'                                         1e-15, lesion.data)
#'
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
