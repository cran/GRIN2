
#' Oncoprint proportions
#'
#' The function order lesion types based on their average size and assign the proportion of the oncoprint rectangle that should be color filled based on the average size of each lesion type.
#'
#' @param lsn.data data.frame with 5 columns including "ID" which is the subject identifier, "chrom" which is the chromosome on which the lesion is located, "loc.start" with lesion start position, "loc.end" which is the lesion end position), and "lsn.type" which is the lesion category for example gain, mutation, etc..).
#' @param clr Lesion colors (If not provided by the user, colors will be automatically assigned using default.grin.colors function).
#' @param hgt Manually assign the proportion of the oncoprint rectangle that should be color filled for each lesion group.
#'
#' @details
#' Some patients might be affected by two or more lesion types in the same gene for example gain AND mutations. To get all lesion types represented in the same rectangle in the oncoprint, this function order lesion types based on the average size of each type and assign the proportion of the oncoprint rectangle that should be color filled based on the average size of each lesion type. Color filled proportion of the oncoprint rectangles can be also specified by the user for  each lesion type based on the hgt parameter.
#'
#' @return
#' Function return a list of three lists specifying the color assigned to each lesion type, the proportion of the rectangle that should be color filled in the oncoprint based on the average size of each lesion type, and the legend parameters.
#'
#' @export
#'
#' @importFrom tibble rownames_to_column
#' @importFrom data.table setorder
#' @import ComplexHeatmap
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author {Lakshmi Patibandla \email{LakshmiAnuhya.Patibandla@stjude.org}, Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @examples
#' data(lesion.data)
#' onco.props=onco.print.props(lesion.data, hgt = c("gain"=4, "loss"=3, "mutation"=2, "fusion"=1))
#' # if hgt argument is not specified, the lesion category "mutation" for single point mutations will
#' # be assigned size=1 because it has the smallest average lesion size and will have the smallest
#' # proportion of the filled oncoprint rectangles 1/4=0.25

onco.print.props=function(lsn.data,  # data.frame with columns ID (subject identifier), chrom (chromosome on which the lesion is located), loc.start (lesion start position), loc.end (lesion end position), lsn.type (lesion category for example gain, mutation, etc..)
                          clr=NULL,  # Lesion colors (If not provided by the user, colors will be automatically assigned using default.grin.colors function).
                          hgt=NULL)  # manually assign the proportion of the oncoprint rectangle that should be color filled for each lesion group
{
  # get unique lesion types
  results <- list()
  lsn.types=unique(lsn.data$lsn.type)

  if (is.null(clr))
  {
    clr=default.grin.colors(sort(lsn.types))
  }
  size=lsn.data$loc.end-lsn.data$loc.start+1
  lsn.data$size=lsn.data$loc.end-lsn.data$loc.start+1
  ave.lsn.size=data.table::setorder(stats::aggregate(size~lsn.type,data=lsn.data,mean),-size)

  # specify proportion of the oncoprint recatangle that should be color filled based on the lesion size
  if (is.null(hgt))
  {
    ave.lsn.size$hgt = length(ave.lsn.size$lsn.type):1
  }
  else{
    hgt=as.data.frame(hgt)
    hgt=tibble::rownames_to_column(hgt, "lsn.type")
    ave.lsn.size=merge(hgt,ave.lsn.size,by="lsn.type", all.y=TRUE)
  }

  ave.lsn.size$wdth = rep(1,length(ave.lsn.size$lsn.type))
  ave.lsn.size=ave.lsn.size[order(ave.lsn.size$lsn.type),]

  twhc=cbind.data.frame(type=ave.lsn.size$lsn.type,
                        clr=default.grin.colors(sort(lsn.types)),
                        hgt=ave.lsn.size$hgt,
                        wdth=ave.lsn.size$wdth)

  res=onco.print.alter.func(twhc)
  res=gsub('\"',"'",res,fixed=T)
  res2=lapply(FUN=eval,X=parse(text=res))
  names(res2) = names(res)
  results$alter_func <- res2
  names(clr) <- ave.lsn.size$lsn.type
  results$col <- sort(clr)
  results$heatmap_legend_param <- oncoprint.legend(ave.lsn.size$lsn.type)
  return(results)

}

