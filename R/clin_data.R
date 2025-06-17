
#' Example Clinical Dataset for T-cell Acute Lymphoblastic Leukemia (T-ALL)
#'
#' This dataset contains clinical and demographic information for 265 newly diagnosed T-ALL patients.
#' The data originates from Liu, Yu, et al. (2017) and includes variables relevant to patient characteristics
#' and clinical outcomes.
#'
#' @format A data frame with 265 rows and 11 columns:
#' \describe{
#'   \item{ID}{Unique patient identifier}
#'   \item{Sex}{Patient gender}
#'   \item{Race}{Patient race}
#'   \item{Age_Days}{Age at diagnosis, in days}
#'   \item{WBC}{White Blood Cell count at diagnosis}
#'   \item{MRD29}{Minimal Residual Disease percentage at day 29 post-treatment}
#'   \item{MRD.binary}{Categorical MRD status (0 = MRD <= 0.1, 1 = MRD > 0.1)}
#'   \item{os.time}{Overall survival time in years (from diagnosis to last follow-up or death)}
#'   \item{os.censor}{Overall survival status (0 = alive at last follow-up, 1 = deceased)}
#'   \item{efs.time}{Event-free survival time in years}
#'   \item{efs.censor}{Event status for event-free survival (0 = censored, 1 = event occurred)}
#' }
#'
#' @source Liu, Yu, et al. (2017), *Nature Genetics*. Supplementary tables:
#' <https://www.nature.com/articles/ng.3909#Sec27>. Additional clinical variables were integrated from the
#' TARGET database (https://ocg.cancer.gov/programs/target).
"clin_data"
