
#' Example T-ALL Dataset Lesion Data
#'
#' Lesion data file showing copy number variations, single nucleotide variations and structrual rearrangments affecting 265 newly diagnosed T-cell Acute Lymphoblastic Leukemia (T-ALL) patients that was reported by Liu, Yu, et al. (2017).
#'
#' @format ## `lesion.data`
#' A data frame with 6,887 rows and 5 columns:
#' \describe{
#'   \item{ID}{patient identifier for the patient affected by the genomic lesion}
#'   \item{chrom}{the chromosome on which the lesion is located}
#'   \item{loc.start}{the lesion start position in base pairs}
#'   \item{loc.end}{the lesion end position in base pairs}
#'   \item{lsn.type}{the lesion type for example gain, loss, mutation, fusion, etc...}
#' }
#' @source extracted from the supplementary material tables of the published Liu, Yu, et al. (2017) manuscript <https://www.nature.com/articles/ng.3909#Sec27>
"lesion.data"
