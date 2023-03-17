#' Class to represent mass spectroscopy data.
#' @slot name A character indicating the name of the experiment.
#' @slot raw A data frame containing the raw data.
#' @slot scaled A data frame containing the scaled data.
#' @slot consolidated A logical indicating whether or not the data has
#' been consolidated.
#' @slot metadata A data frame containing the metadata.
#' @rdname MSObject-class
#' @name MSObject-class
#' @concept objects
#' @exportClass MSObject
setClass(
  Class = "MSObject",
  slots = list(
    name = "character",
    raw = "tbl_df",
    scaled = "tbl_df",
    consolidated = "logical",
    metadata = "tbl_df"
  )
)

#' Class to represent merged mass spectroscopy data.
#' @slot ms1 A character indicating the name of the experiment.
#' @slot ms2 A data frame containing the raw data.
#' @slot merged A data frame containing the scaled data.
#' @slot metadata A logical indicating whether or not the data has
#' been consolidated.
#' @rdname MergedMSObject-class
#' @name MergedMSObject-class
#' @concept objects
#' @exportClass MergedMSObject
setClass(
  Class = "MergedMSObject",
  slots = list(
    ms1 = "MSObject",
    ms2 = "MSObject",
    merged = "tbl_df",
    metadata = "tbl_df"
  )
)

create_ms_obj <- function(df,
                          name,
                          id_name,
                          rt_name = NULL,
                          mz_name = NULL,
                          int_name = NULL,
                          has_metadata = FALSE) {
  ms <- new("MSObject",
    name = name,
    consolidated = FALSE, raw = NULL,
    scaled = NULL, metadata = NULL
  )

  raw_data <- df |>
    dplyr::select(id_name, rt_name, mz_name, int_name)
  colnames(raw_data) <- c("id", "rt", "mz", "int")
  ms@raw <- raw_data

  if (has_metadata) {
    meta_data <- df |>
      dplyr::select(id_name, -rt_name, -mz_name, -int_name)
    ms@metadata <- meta_data
  }
  return(ms)
}
