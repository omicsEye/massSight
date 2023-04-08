#' Class to represent mass spectroscopy data.
#' @slot name A character indicating the name of the experiment.
#' @slot raw A data frame containing the raw data.
#' @slot isolated A data frame containing the isolated data.
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
    isolated = "tbl_df",
    scaled_df = "tbl_df",
    consolidated = "logical",
    metadata = "tbl_df"
  )
)

#' @export
setGeneric("name", function(x) standardGeneric("name"))
setMethod("name", signature = "MSObject", definition = function(x) x@name)

#' @export
setGeneric("raw_df", function(x) standardGeneric("raw_df"))
setMethod("raw_df", signature = "MSObject", definition = function(x) x@raw_df)

#' @export
setGeneric("isolated", function(x) standardGeneric("isolated"))
setMethod("isolated", signature = "MSObject", definition = function(x) x@isolated)

#' @export
setGeneric("scaled_df", function(x) standardGeneric("scaled_df"))
setMethod("scaled_df", signature = "MSObject", definition = function(x) x@scaled_df)

#' @export
setGeneric("consolidated", function(x) standardGeneric("consolidated"))
setMethod("consolidated", signature = "MSObject", definition = function(x) x@consolidated)

#' @export
setGeneric("metadata", function(x) standardGeneric("metadata"))
setMethod("metadata", signature = "MSObject", definition = function(x) x@metadata)

#' @export
setGeneric("name<-", function(x, value) standardGeneric("name<-"))
setMethod("name<-", signature = "MSObject", definition = function(x, value) {
  x@name <- value
  return(x)
})

#' @export
setGeneric("raw_df<-", function(x, value) standardGeneric("raw_df<-"))
setMethod("raw_df<-", signature = "MSObject", definition = function(x, value) {
  x@raw_df <- value
  return(x)
})

#' @export
setGeneric("isolated<-", function(x, value) standardGeneric("isolated<-"))
setMethod("isolated<-", signature = "MSObject", definition = function(x, value) {
  x@isolated <- value
  return(x)
})

#' @export
setGeneric("scaled_df<-", function(x, value) standardGeneric("scaled_df<-"))
setMethod("scaled_df<-", signature = "MSObject", definition = function(x, value) {
  x@scaled_df <- value
  return(x)
})

#' @export
setGeneric("consolidated<-", function(x, value) standardGeneric("consolidated<-"))
setMethod("consolidated<-", signature = "MSObject", definition = function(x, value) {
  x@consolidated <- value
  return(x)
})

#' @export
setGeneric("metadata<-", function(x, value) standardGeneric("metadata<-"))
setMethod("metadata<-", signature = "MSObject", definition = function(x, value) {
  x@metadata <- value
  return(x)
})

#' Class to represent merged mass spectroscopy data.
#' @slot ms1 A character indicating the name of the experiment.
#' @slot ms2 A data frame containing the raw data.
#' @slot all_matched A data frame containing the scaled data.
#' @slot iso_matched A data frame containing the scaled data.
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
    all_matched = "tbl_df",
    iso_matched = "tbl_df",
    scaled_values = "tbl_df",
    cutoffs = "numeric",
    aligned = "tbl_df",
    metadata = "tbl_df",
    smooth = "list"
  )
)

#' @export
setGeneric("ms1", function(x) standardGeneric("ms1"))
setMethod("ms1", signature = "MergedMSObject", definition = function(x) x@ms1)

#' @export
setGeneric("ms2", function(x) standardGeneric("ms2"))
setMethod("ms2", signature = "MergedMSObject", definition = function(x) x@ms2)

#' @export
setGeneric("all_matched", function(x) standardGeneric("all_matched"))
setMethod("all_matched", signature = "MergedMSObject", definition = function(x) x@all_matched)

#' @export
setGeneric("iso_matched", function(x) standardGeneric("iso_matched"))
setMethod("iso_matched", signature = "MergedMSObject", definition = function(x) x@iso_matched)

#' @export
setGeneric("scaled_values", function(x) standardGeneric("scaled_values"))
setMethod("scaled_values", signature = "MergedMSObject", definition = function(x) x@scaled_values)

#' @export
setGeneric("cutoffs", function(x) standardGeneric("cutoffs"))
setMethod("cutoffs", signature = "MergedMSObject", definition = function(x) x@cutoffs)

#' @export
setGeneric("aligned", function(x) standardGeneric("aligned"))
setMethod("aligned", signature = "MergedMSObject", definition = function(x) x@aligned)

#' @export
setGeneric("metadata", function(x) standardGeneric("metadata"))
setMethod("metadata", signature = "MergedMSObject", definition = function(x) x@metadata)

#' @export
setGeneric("smooth", function(x) standardGeneric("smooth"))
setMethod("smooth", signature = "MergedMSObject", definition = function(x) x@smooth)

#' @export
setGeneric("ms1<-", function(x, value) standardGeneric("ms1<-"))
setMethod("ms1<-", signature = "MergedMSObject", definition = function(x, value) {
  x@ms1 <- value
  return(x)
})

#' @export
setGeneric("ms2<-", function(x, value) standardGeneric("ms2<-"))
setMethod("ms2<-", signature = "MergedMSObject", definition = function(x, value) {
  x@ms2 <- value
  return(x)
})

#' @export
setGeneric("all_matched<-", function(x, value) standardGeneric("all_matched<-"))
setMethod("all_matched<-", signature = "MergedMSObject", definition = function(x, value) {
  x@all_matched <- value
  return(x)
})

#' @export
setGeneric("iso_matched<-", function(x, value) standardGeneric("iso_matched<-"))
setMethod("iso_matched<-", signature = "MergedMSObject", definition = function(x, value) {
  x@iso_matched <- value
  return(x)
})

#' @export
setGeneric("scaled_values<-", function(x, value) standardGeneric("scaled_values<-"))
setMethod("scaled_values<-", signature = "MergedMSObject", definition = function(x, value) {
  x@scaled_values <- value
  return(x)
})

#' @export
setGeneric("cutoffs<-", function(x, value) standardGeneric("cutoffs<-"))
setMethod("cutoffs<-", signature = "MergedMSObject", definition = function(x, value) {
  x@cutoffs <- value
  return(x)
})

#' @export
setGeneric("aligned<-", function(x, value) standardGeneric("aligned<-"))
setMethod("aligned<-", signature = "MergedMSObject", definition = function(x, value) {
  x@aligned <- value
  return(x)
})

#' @export
setGeneric("metadata<-", function(x, value) standardGeneric("metadata<-"))
setMethod("metadata<-",
  signature = "MergedMSObject",
  definition = function(x, value) {
    x@metadata <- value
    return(x)
  }
)

#' @export
setGeneric("smooth<-", function(x, value) standardGeneric("smooth<-"))
setMethod("smooth<-",
  signature = "MergedMSObject",
  definition = function(x, value) {
    x@smooth <- value
    return(x)
  }
)

#' @export
#' @title Create MS Object
#' @description Create an MSObject from a data frame.
#' @param df A data frame containing the raw data.
#' @param name A character indicating the name of the experiment.
#' @param id_name A character indicating the name of the column containing the
#' compound IDs.
#' @param rt_name A character indicating the name of the column containing the
#' retention times.
#' @param mz_name A character indicating the name of the column containing the
#' m/z values.
#' @param int_name A character indicating the name of the column containing the
#' intensities.
#' @param has_metadata A logical indicating whether or not the data has
#' metadata.
#' @return An MSObject.
create_ms_obj <- function(df,
                          name,
                          id_name,
                          rt_name = NULL,
                          mz_name = NULL,
                          int_name = NULL,
                          has_metadata = FALSE) {
  ms <- new("MSObject")
  name(ms) <- name
  consolidated(ms) <- FALSE

  raw_data <- df |>
    dplyr::as_tibble() |>
    dplyr::select(id_name, rt_name, mz_name, int_name)
  colnames(raw_data) <- c("Compound_ID", "RT", "MZ", "Intensity")
  raw_df(ms) <- raw_data

  if (has_metadata) {
    meta_data <- df |>
      dplyr::select(id_name, -rt_name, -mz_name, -int_name)
    metadata(ms) <- meta_data
  }
  return(ms)
}

create_aligned_ms_obj <- function(ms1, ms2) {
  ms <- new("MergedMSObject",
    ms1 = ms1,
    ms2 = ms2,
    all_matched = NULL,
    iso_matched = NULL,
    scaled_values = NULL,
    cutoffs = NULL,
    aligned = NULL,
    metadata = NULL,
    smooth = NULL
  )
  return(ms)
}
