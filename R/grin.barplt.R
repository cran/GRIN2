
#' GRIN Bar Plot
#'
#' @description
#' Function return a stacked bar plot with number of patients affected by all different types of lesions in a pre-specified list of genes of interest.
#'
#' @param grin.res GRIN results (output of the grin.stats function).
#' @param count.genes vector with gene names of a list of genes to be added to the bar plot.
#' @param lsn.colors Lesion colors (If not provided by the user, colors will be automatically assigned using default.grin.colors function).
#'
#' @details
#' Function will use the input list of gene names and extract the number of patients affected by all different types of lesions in those genes from the GRIN results table (output of the grin.stats function).
#'
#' @return
#' Function return a stacked bar plot with number of patients affected by all different types of lesions in the pre-specified list of genes of interest.
#'
#' @export
#'
#' @importFrom dplyr select group_by summarize arrange mutate
#' @importFrom ggplot2 geom_bar position_stack scale_fill_manual geom_text labs coord_flip theme_classic element_text theme ggplot
#' @importFrom utils stack
#' @importFrom stats na.omit reorder
#' @importFrom tidyselect starts_with
#' @importFrom magrittr %>%
#'
#' @author {Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [grin.stats()]
#'
#' @examples
#' data(lesion.data)
#' data(hg19.gene.annotation)
#' data(hg19.chrom.size)
#'
#'# run GRIN analysis using grin.stats function
#' grin.results=grin.stats(lesion.data,
#'                         hg19.gene.annotation,
#'                         hg19.chrom.size)
#'
#' # specify a list of genes to be included in the bar plot (driver genes)
#' count.genes=as.vector(c("TAL1", "FBXW7", "PTEN", "IRF8","NRAS",
#'                         "BCL11B", "MYB", "LEF1","RB1", "MLLT3", "EZH2", "ETV6", "CTCF",
#'                         "JAK1", "KRAS", "RUNX1", "IKZF1", "KMT2A", "RPL11", "TCF7",
#'                         "WT1", "JAK2", "JAK3", "FLT3"))
#' # return the stacked barplot
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
