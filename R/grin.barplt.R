
#' GRIN Lesion Stacked Bar Plot
#'
#' @description
#' Generates a stacked bar plot showing the number of patients affected by different types of genomic lesions in a user-specified list of genes of interest, based on GRIN analysis results.
#'
#' @param grin.res A data frame of GRIN results, typically the output from the \code{\link{grin.stats}} function.
#' @param count.genes A character vector of gene names to include in the bar plot. Only genes present in the GRIN results table will be used.
#' @param lsn.colors (Optional) A named vector specifying colors for each lesion type. If not provided, default lesion colors will be automatically assigned using the \code{\link{default.grin.colors}} function.
#'
#' @details
#' The function subsets the GRIN results to the genes specified in \code{count.genes}, extracts the number of patients affected by each lesion type for each gene, and visualizes the data as a horizontal stacked bar plot. Each bar represents a gene and is segmented by lesion type (e.g., gain, loss, mutation), with segment size proportional to the number of affected patients.
#'
#' This visualization is useful for highlighting the burden and distribution of different lesion types across key driver genes or other genes of interest.
#'
#' @return
#' A ggplot2-based stacked bar plot showing lesion type distribution across the selected genes.
#'
#' @export
#'
#' @importFrom dplyr select group_by summarize arrange mutate
#' @importFrom ggplot2 ggplot geom_bar geom_text scale_fill_manual position_stack labs coord_flip theme_classic theme element_text
#' @importFrom utils stack
#' @importFrom stats na.omit reorder
#' @importFrom tidyselect starts_with
#' @importFrom magrittr %>%
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link{grin.stats}}, \code{\link{default.grin.colors}}
#'
#' @examples
#' data(lesion_data)
#' data(hg38_gene_annotation)
#' data(hg38_chrom_size)
#'
#' # Run GRIN analysis
#' grin.results <- grin.stats(lesion_data,
#'                            hg38_gene_annotation,
#'                            hg38_chrom_size)
#'
#' # Define a list of genes to be included in the bar plot (e.g., candidate driver genes)
#' count.genes <- c("TAL1", "FBXW7", "PTEN", "IRF8", "NRAS",
#'                  "BCL11B", "MYB", "LEF1", "RB1", "MLLT3",
#'                  "EZH2", "ETV6", "CTCF", "JAK1", "KRAS",
#'                  "RUNX1", "IKZF1", "KMT2A", "RPL11", "TCF7",
#'                  "WT1", "JAK2", "JAK3", "FLT3")
#'
#' # Generate the bar plot
#' grin.barplt(grin.results, count.genes)

grin.barplt=function(grin.res,        # GRIN results (output of the grin.stats function)
                     count.genes,     # vector with gene names of a list of genes to be added to the bar plot
                     lsn.colors=NULL) # Lesion colors (If not provided by the user, colors will be automatically assigned using default.grin.colors function).


{
  hits=grin.res$gene.hits

  # extract count data for all lesion types from GRIN results table
  count.clms=hits %>% dplyr::select(tidyselect::starts_with("nsubj"))
  gene.name=hits$gene.name
  count.data=as.data.frame(cbind(gene.name, count.clms))
  count.genes=count.genes

  # limit the bar plot to the list of genes of interest
  selected.count=count.data[count.data$gene.name%in%count.genes,]
  names(selected.count) <- sub('^nsubj.', '', names(selected.count))
  count.clms=selected.count[,-1]
  count.clms=utils::stack(count.clms)

  nsubj=count.clms$values
  lsn.type=count.clms$ind

  nsubj.data=cbind.data.frame(gene.name=selected.count$gene.name,
                              nsubj=count.clms$values,
                              lsn.type=count.clms$ind)

  nsubj.data[nsubj.data == 0] <- NA
  nsubj.data=stats::na.omit(nsubj.data)

  # assign colors for lesion groups automatically if not provided by the user
  lsn.types=sort(unique(grin.res$lsn.index$lsn.type))

  if (is.null(lsn.colors))
  {
    lsn.colors=default.grin.colors(lsn.types)
  }

  label_pos = cumsum(nsubj) - (nsubj / 2)
  # summarize the df
  nsubj.data <- nsubj.data %>%
    dplyr::group_by(gene.name, lsn.type) %>%
    dplyr::summarize(nsubj = sum(nsubj), .groups = "drop") %>%
    dplyr::group_by(gene.name) %>%
    dplyr::arrange(gene.name, lsn.type) %>%
    dplyr::mutate(label_pos = cumsum(nsubj) - (nsubj / 2))

  # generate stacked bar plot for number of patients affected by each type of lesions in the list of genes of interest
  plot=ggplot2::ggplot(data = nsubj.data,
              aes(x = stats::reorder(gene.name, nsubj, sum)), cex.axis=10) +
    ggplot2::geom_bar(aes(y = nsubj, fill = lsn.type),
                      position = ggplot2::position_stack(reverse = T),
                      stat="identity", width = .5)+
    ggplot2::scale_fill_manual(values=c(lsn.colors))+
    ggplot2::geom_text(aes(y = label_pos, label = nsubj),
                       color = "white", size = 3.5) +
    ggplot2::labs(x="Gene Name",y="Count", fill = "Lesion Type")+
    ggplot2::coord_flip() + ggplot2::theme_classic()

  final.plt=plot+ggplot2::theme(text = ggplot2::element_text(size = 15))

  return(final.plt)

}
