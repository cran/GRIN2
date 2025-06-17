
#' Assign Default GRIN Colors
#'
#' @description
#' Assigns a default set of colors to lesion types for use in GRIN plots.
#'
#' @param lsn.types A character vector of unique lesion types, typically derived from the lesion data.
#'
#' @details
#' This function provides a predefined palette of up to 10 distinct colors for lesion types used in GRIN visualizations. If more than 10 lesion types are provided, the function will prompt the user to manually define custom colors to ensure visual distinction.
#'
#' @return
#' A named character vector of colors corresponding to each lesion type.
#'
#' @export
#'
#' @references
#' Pounds, S. B., et al. (2013). A genomic random interval model for statistical analysis of genomic lesion data.
#'
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @examples
#' data(lesion_data)
#'
#' # Extract unique lesion types
#' lsn.types <- unique(lesion_data$lsn.type)
#'
#' # Assign default colors to lesion types
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
