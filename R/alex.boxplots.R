
#' Prepare Box Plots of Expression Data by Lesion Groups
#'
#' @description
#' Function return box plots for expression data by lesion groups for selected number of genes based on a specified q-value of the kruskal-wallis test results.
#'
#' @param out.dir Path to the folder where the boxplots of selected genes based on the specified q value of the KW results table will be added.
#' @param alex.data Output of the alex.prep.lsn.expr function. It's a list of three data tables that include "row.mtch", "alex.expr" with expression data, "alex.lsn" with lesion data. Rows of alex.expr, and "alex.lsn" matrices are ordered by gene ensembl IDs and the columns are ordered by patient ID.
#' @param alex.kw.results ALEX Kruskal-Wallis test results (output of the KW.hit.express function).
#' @param q minimum q value for a gene to be included in output PDF file of box plots.
#' @param gene.annotation Gene annotation data either provided by the user or retrieved from ensembl BioMart database using get.ensembl.annotation function included in the GRIN2.0 library. Data.frame should has four columns: "gene" which is the ensembl ID of annotated genes, "chrom" which is the chromosome on which the gene is located, "loc.start" which is the gene start position, and "loc.end" the gene end position.
#'
#' @return
#' Function return a PDF file with box plots for expression data by lesion groups for selected number of genes based on a specified q-value of the kruskal-wallis test results (one gene per page).
#'
#' @export
#'
#' @importFrom ggplot2 aes geom_boxplot geom_jitter position_jitter xlab ylab theme_bw scale_fill_discrete guide_legend theme element_text ggplot
#' @importFrom forcats fct_reorder
#' @importFrom grDevices dev.off pdf
#'
#' @author {Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023) Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @seealso [alex.prep.lsn.expr()], [KW.hit.express()]
#'
#' @examples
#' data(expr.data)
#' data(lesion.data)
#' data(hg19.gene.annotation)
#'
#' # prepare expression, lesion data and return the set of genes with both types of data available
#' # ordered by gene IDs in rows and patient IDs in columns:
#' alex.data=alex.prep.lsn.expr(expr.data, lesion.data,
#'                              hg19.gene.annotation, min.expr=5, min.pts.lsn=5)
#'
#' # run KW test for association between lesion groups and expression level of the same gene:
#' alex.kw.results=KW.hit.express(alex.data, hg19.gene.annotation, min.grp.size=5)
#'
#' # return boxplots for a list of top significant genes to a pre-specified folder using 'out.dir':
#' dir.create(resultsFolder <- file.path(tempdir(), "temp.out"))
#'
#' boxplots=alex.boxplots(out.dir=resultsFolder,
#'                        alex.data, alex.kw.results,
#'                        1e-15, hg19.gene.annotation)
#'
#' unlink(resultsFolder, recursive = TRUE)

alex.boxplots=function(out.dir,            # Path to the folder where the boxplots of selected genes based on the specified q value of the KW results table will be added
                       alex.data,          # output of the alex.prep.lsn.expr function (list of three data tables "alex.expr" with expression data, "alex.lsn" with overlapped gene lesion data and row.mtch)
                       alex.kw.results,    # ALEX Kruskal-Wallis test results (output of the KW.hit.express function)
                       q,                  # minimum q value for a gene to be included in output PDF file of box plots
                       gene.annotation)    # gene annotation data

{
  # To extract genes below specified q.value for boxplots:
  KW.results=alex.kw.results[!is.na(alex.kw.results$q.KW),]
  selected.genes=KW.results[KW.results$q.KW<q,]

  selected.IDS=selected.genes$gene
  selected.lsns=alex.data$alex.lsn[rownames(alex.data$alex.lsn) %in% selected.IDS,]
  selected.lsns=t(selected.lsns)
  selected.expr=alex.data$alex.expr[rownames(alex.data$alex.expr) %in% selected.IDS,]
  selected.expr=t(selected.expr)

  # To make sure that both lesion and expression matrices have same patients and genes order
  check.rownames=all(rownames(selected.lsns)==rownames(selected.expr))
  if (check.rownames==FALSE)
    stop("Gene-lesion matrix rownames must match expression matrix rownames.")

  check.colnames=all(colnames(selected.lsns)==colnames(selected.expr))
  if (check.colnames==FALSE)
    stop("Gene-lesion matrix patient IDs must match expression matrix patient IDs.")

  # add ensembl.ID to gene.name column if empty
  #gene.annotation=gene.annotation %>% mutate_all(na_if,"")
  gene.annotation$gene.name <- ifelse(is.na(gene.annotation$gene.name),gene.annotation$gene, gene.annotation$gene.name)

  n=ncol(selected.lsns)

  for (i in 1:n)
  {
    ensembl.id=colnames(selected.lsns)[i]
    ensembl.df=as.data.frame(ensembl.id)
    colnames(ensembl.df) <- ("gene")
    genes.df.merged=merge(gene.annotation,ensembl.df,by="gene", all.y=TRUE)
    gene.names=as.character(genes.df.merged$gene.name)
    lsn.box=unlist(selected.lsns[,i])
    expr.box=unlist(selected.expr[,i])
    df=as.data.frame(cbind(lsn.box, expr.box))
    df$expr.box=as.numeric(df$expr.box)

    {
      grDevices::pdf(paste0(out.dir,gene.names,"_boxplot", ".pdf"),
                     height=8,width=10)
      p=ggplot2::ggplot(df, ggplot2::aes(x = forcats::fct_reorder(lsn.box, expr.box, .desc =TRUE), y = expr.box)) +
        ggplot2::geom_boxplot(ggplot2::aes(fill = forcats::fct_reorder(lsn.box, expr.box,  .desc =TRUE))) +
        ggplot2::geom_jitter(position=ggplot2::position_jitter(0.02)) + # Increasing the number from 0.02 increase the space between the dots but add some extra dots if number of patients is small
        ggplot2::xlab(paste0(gene.names, '_lsn'))+
        ggplot2::ylab(paste0(gene.names, '_Expression'))+
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::scale_fill_discrete(guide = ggplot2::guide_legend(title = "Lesion")) +
        ggplot2::theme(axis.text=ggplot2::element_text(size=10), axis.title=ggplot2::element_text(size=12,face="bold"))+
        ggplot2::theme(legend.text=ggplot2::element_text(size=8))
      print(p)
      grDevices::dev.off()
    }
  }
}
