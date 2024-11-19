
#' Default GRIN Colors
#'
#' @description
#' Function assigns default colors for each lesion group in the whole set of GRIN plots.
#'
#' @param lsn.types Unique lesion types as specified in the lesion data file.
#'
#' @details
#' The function specifies 10 colors for different lesion types. If the number of lesion types is more than 10, the user will be asked to specify the colors manually.
#'
#' @return
#' Function return a vector of colors assigned to each unique lesion type.
#'
#' @export
#'
#' @references
#' Pounds, Stan, et al. (2013) A genomic random interval model for statistical analysis of genomic lesion data.
#'
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author {Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @examples
#' data(lesion.data)
#'
#' lsn.types=unique(lesion.data$lsn.type)
#' # assign colors for different lesion categories using default.grin.colors function:
#' default.grin.colors(lsn.types)

default.grin.colors=function(lsn.types) # Unique lesion types as specified in the lesion data file
{
  message("Computationally assigning lesion type colors for GRIN plots.")
  uniq.types=sort(unique(lsn.types))
  n.types=length(uniq.types)
  default.colors=c("black", "red","blue",
                   "olivedrab", "purple",
                   "cyan", "brown", "gold",
                   "orange","steelblue")
  if (length(n.types)>length(default.colors))
    stop(paste0("Too many lesion types for default grin colors; please assign colors manually."))

  res=default.colors[1:n.types]
  names(res)=uniq.types
  return(res)

}
