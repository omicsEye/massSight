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
    raw_df = "tbl_df",
    scaled_df = "tbl_df",
    consolidated = "logical",
    metadata = "tbl_df"
  )
)

setMethod("name", signature = "MSObject", definition = function(x) x@name)
setMethod("name<-", signature = "MSObject", definition = function(x, value) {
  x@name <- value
  return(x)
})
setMethod("raw_df", signature = "MSObject", definition = function(x) x@raw)
setMethod("raw_df<-", signature = "MSObject", definition = function(x, value) {
  x@raw_df <- value
  return(x)
})
setMethod("scaled_df", signature = "MSObject", definition = function(x) x@scaled)
setMethod("scaled_df<-", signature = "MSObject", definition = function(x, value) {
  x@scaled_df <- value
  return(x)
})
setMethod("consolidated",
  signature = "MSObject",
  definition = function(x) x@consolidated
)
setMethod("consolidated<-",
  signature = "MSObject",
  definition = function(x, value) {
    x@consolidated <- value
    return(x)
  }
)
setMethod("metadata",
  signature = "MSObject",
  definition = function(x) x@metadata
)
setMethod("metadata<-",
  signature = "MSObject",
  definition = function(x, value) {
    x@metadata <- value
    return(x)
  }
)

#' Class to represent merged mass spectroscopy data.
#' @slot ms1 A character indicating the name of the experiment.
#' @slot ms2 A data frame containing the raw data.
#' @slot merged A data frame containing the scaled data.
#' @slot metadata A logical indicating whether or not the data has
#' been consolidated.
#' @slot smooth A string indicating the smoothing method used.
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
    metadata = "tbl_df",
    smooth = "character"
  )
)

create_ms_obj <- function(df,
                          path,
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

  name(ms) <- name
  consolidated(ms) <- FALSE

  raw_data <- df |>
    dplyr::select(id_name, rt_name, mz_name, int_name)
  colnames(raw_data) <- c("id", "rt", "mz", "int")
  raw_df(ms) <- raw_data

  if (has_metadata) {
    meta_data <- df |>
      dplyr::select(id_name, -rt_name, -mz_name, -int_name)
    metadata(ms) <- meta_data
  }
  return(ms)
}
