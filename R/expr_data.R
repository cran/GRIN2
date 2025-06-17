
#' Example T-ALL Gene Expression Dataset
#'
#' Log2-normalized gene expression data for 417 genes across 265 newly diagnosed T-cell Acute Lymphoblastic Leukemia (T-ALL) patients, as reported by Liu, Yu, et al. (2017).
#'
#' @format ## `expr_data`
#' A data frame with 417 rows and 265 columns:
#' \describe{
#'   \item{gene}{Ensembl gene IDs of the 417 selected genes included in the dataset.}
#'   \item{...}{Each remaining column represents a T-ALL patient with log2-normalized expression values.}
#' }
#'
#' @source Data extracted from the supplementary materials of Liu, Yu, et al. (2017), *Nature Genetics*. <https://www.nature.com/articles/ng.3909#Sec27>
"expr_data"
