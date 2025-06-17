
#' Oncoprint Proportions by Lesion Type
#'
#' @description
#' Calculates and assigns the proportion of each oncoprint rectangle to be color-filled based on the average size of lesion types. Lesion types are ordered by their average genomic size, and proportions are either computed automatically or manually specified by the user.
#' @param lsn.data A data frame with five columns:
#' \itemize{
#'   \item \code{ID}: Subject or patient identifier
#'   \item \code{chrom}: Chromosome on which the lesion is located
#'   \item \code{loc.start}: Start genomic position of the lesion
#'   \item \code{loc.end}: End genomic position of the lesion
#'   \item \code{lsn.type}: Lesion category (e.g., gain, mutation, fusion)
#' }
#' @param clr Optional. A named vector of colors for each lesion type. If not provided, default colors will be assigned using \code{\link{default.grin.colors}}.
#' @param hgt Optional. A named numeric vector specifying the proportion (height) of the oncoprint rectangle to be filled for each lesion type. If not provided, proportions will be determined automatically based on average lesion sizes.
#'
#' @details
#' In cases where a patient has multiple types of lesions (e.g., gain and mutation) in the same gene, this function ensures that all lesion types are visually represented within a single oncoprint rectangle.
#'
#' If \code{hgt} is not specified, lesion types are ranked by their average genomic size (calculated as \code{loc.end - loc.start + 1}), and the oncoprint proportions are derived accordingly. Smaller lesions (like point mutations) will occupy a smaller portion of the rectangle, while larger lesions (like CNVs) will occupy a larger portion.
#'
#' Alternatively, the user can manually define the fill proportions via the \code{hgt} parameter.
#'
#' @return
#' A list with the following components:
#' \itemize{
#'   \item \code{col}: Named vector of colors assigned to each lesion type
#'   \item \code{hgt}: Named vector of normalized proportions to fill for each lesion type
#'   \item \code{legend}: Legend parameters for use in oncoprint plots
#' }
#'
#' @export
#'
#' @importFrom tibble rownames_to_column
#' @importFrom data.table setorder
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Lakshmi Patibandla \email{LakshmiAnuhya.Patibandla@stjude.org},
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org},
#' Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @examples
#' data(lesion_data)
#'
#' # Automatically assign oncoprint proportions based on average lesion size:
#' onco.props <- onco.print.props(lesion_data)
#'
#' # Alternatively, manually specify the oncoprint fill proportions for each lesion type:
#' onco.props <- onco.print.props(lesion_data,
#'                         hgt = c("gain" = 4, "loss" = 3, "mutation" = 2, "fusion" = 1))

onco.print.props=function (lsn.data, clr = NULL, hgt = NULL)
{
  results <- list()
  lsn.types = unique(lsn.data$lsn.type)
  if (is.null(clr)) {
    clr = default.grin.colors(sort(lsn.types))
  }
  size = lsn.data$loc.end - lsn.data$loc.start + 1
  lsn.data$size = lsn.data$loc.end - lsn.data$loc.start + 1
  ave.lsn.size = data.table::setorder(stats::aggregate(size ~
                                                         lsn.type, data = lsn.data, mean), -size)
  if (is.null(hgt)) {
    ave.lsn.size$hgt = length(ave.lsn.size$lsn.type):1
  }
  else {
    hgt = as.data.frame(hgt)
    hgt = tibble::rownames_to_column(hgt, "lsn.type")
    ave.lsn.size = merge(hgt, ave.lsn.size, by = "lsn.type",
                         all.y = TRUE)
  }
  ave.lsn.size$wdth = rep(1, length(ave.lsn.size$lsn.type))
  ave.lsn.size = ave.lsn.size[order(ave.lsn.size$lsn.type),
  ]
  twhc = cbind.data.frame(type = ave.lsn.size$lsn.type, clr = clr,
                          hgt = ave.lsn.size$hgt, wdth = ave.lsn.size$wdth)
  res = onco.print.alter.func(twhc)
  res = gsub("\"", "'", res, fixed = T)
  res2 = lapply(FUN = eval, X = parse(text = res))
  names(res2) = names(res)
  results$alter_func <- res2
  names(clr) <- ave.lsn.size$lsn.type
  results$col <- sort(clr)
  results$heatmap_legend_param <- oncoprint.legend(ave.lsn.size$lsn.type)
  return(results)
}
