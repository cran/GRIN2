
#' Generate Box Plots of Gene Expression by Lesion Groups
#'
#' @description
#' Generates box plots of gene expression levels stratified by lesion groups for a subset of genes selected based on a user-specified q-value threshold from Kruskal Wallis test results.
#'
#' @param out.dir Path to the directory where the resulting PDF file containing the box plots will be saved.
#' @param alex.data A list returned by the \code{alex.prep.lsn.expr} function, containing three components: \code{"row.mtch"}, \code{"alex.expr"} (gene expression matrix), and \code{"alex.lsn"} (lesion matrix). Rows of \code{alex.expr} and \code{alex.lsn} is ordered by Ensembl gene IDs; columns represent patient IDs.
#' @param alex.kw.results Output of the \code{KW.hit.express} function containing Kruskal Wallis test results.
#' @param q Q-value threshold. Only genes with q-values below this threshold will be included in the output plots.
#' @param gene.annotation A \code{data.frame} containing gene annotation data, either provided by the user or retrieved using \code{get.ensembl.annotation}. Must contain four columns: \code{"gene"} (Ensembl gene ID), \code{"chrom"} (chromosome), \code{"loc.start"} (start position), and \code{"loc.end"} (end position).
#'
#' @return
#' A PDF file saved in \code{out.dir}, containing one box plot per page for each selected gene, showing gene expression distribution across lesion groups.
#'
#' @export
#'
#' @importFrom ggplot2 aes geom_boxplot geom_jitter position_jitter xlab ylab theme_bw scale_fill_discrete guide_legend theme element_text ggplot
#' @importFrom forcats fct_reorder
#' @importFrom grDevices dev.off pdf
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @seealso \code{\link{alex.prep.lsn.expr}}, \code{\link{KW.hit.express}}
#'
#' @examples
#' data(expr_data)
#' data(lesion_data)
#' data(hg38_gene_annotation)
#'
#' # Prepare expression and lesion data
#' alex.data <- alex.prep.lsn.expr(expr_data, lesion_data,
#'                                 hg38_gene_annotation, min.expr = 5, min.pts.lsn = 5)
#'
#' # Run Kruskal Wallis test
#' alex.kw.results <- KW.hit.express(alex.data, hg38_gene_annotation, min.grp.size = 5)
#'
#' # Generate box plots for significant genes
#' dir.create(resultsFolder <- file.path(tempdir(), "temp.out"))
#' alex.boxplots(out.dir = resultsFolder,
#'               alex.data = alex.data,
#'               alex.kw.results = alex.kw.results,
#'               q = 1e-15,
#'               gene.annotation = hg38_gene_annotation)
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
