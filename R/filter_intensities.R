#' Filter Intensities
#'
#' @param data Dataframe containing results from an LC-MS experiment.
#' @param prevalence Percent of samples required to keep a metabolite.
#'
#' @return A logical vector indicating which rows of the original dataset
#' should be retained
#' @export
filter_intensities <- function(data, prevalence = .5) {
  row_nas <- data |>
    is.na() |>
    rowSums()
  keep <- (row_nas / ncol(data)) > prevalence
  return(keep)
}
