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
    pre_iso_matched = "data.frame",
    iso_matched = "data.frame",
    all_matched = "data.frame",
    scaled_values = "data.frame",
    adjusted_df = "data.frame",
    cutoffs = "numeric",
    aligned = "data.frame",
    metadata = "data.frame",
    smooth = "character"
  )
)

#' @export
setGeneric("ms1", function(x) {
  standardGeneric("ms1")
})
setMethod(
  "ms1",
  signature = "MergedMSObject",
  definition = function(x) {
    x@ms1
  }
)

#' @export
setGeneric("ms2", function(x) {
  standardGeneric("ms2")
})
setMethod(
  "ms2",
  signature = "MergedMSObject",
  definition = function(x) {
    x@ms2
  }
)


#' @export
setGeneric("pre_iso_matched", function(x) {
  standardGeneric("pre_iso_matched")
})
setMethod(
  "pre_iso_matched",
  signature = "MergedMSObject",
  definition = function(x) {
    x@pre_iso_matched
  }
)

#' @export
setGeneric("all_matched", function(x) {
  standardGeneric("all_matched")
})
setMethod(
  "all_matched",
  signature = "MergedMSObject",
  definition = function(x) {
    x@all_matched
  }
)

#' @export
setGeneric("iso_matched", function(x) {
  standardGeneric("iso_matched")
})
setMethod(
  "iso_matched",
  signature = "MergedMSObject",
  definition = function(x) {
    x@iso_matched
  }
)

#' @export
setGeneric("scaled_values", function(x) {
  standardGeneric("scaled_values")
})
setMethod(
  "scaled_values",
  signature = "MergedMSObject",
  definition = function(x) {
    x@scaled_values
  }
)

#' @export
setGeneric("adjusted_df", function(x) {
  standardGeneric("adjusted_df")
})
setMethod(
  "adjusted_df",
  signature = "MergedMSObject",
  definition = function(x) {
    x@adjusted_df
  }
)

#' @export
setGeneric("cutoffs", function(x) {
  standardGeneric("cutoffs")
})
setMethod(
  "cutoffs",
  signature = "MergedMSObject",
  definition = function(x) {
    x@cutoffs
  }
)

#' @export
setGeneric("aligned", function(x) {
  standardGeneric("aligned")
})
setMethod(
  "aligned",
  signature = "MergedMSObject",
  definition = function(x) {
    x@aligned
  }
)

#' @export
setGeneric("metadata", function(x) {
  standardGeneric("metadata")
})
setMethod(
  "metadata",
  signature = "MergedMSObject",
  definition = function(x) {
    x@metadata
  }
)

#' @export
setGeneric("smooth", function(x) {
  standardGeneric("smooth")
})
setMethod(
  "smooth",
  signature = "MergedMSObject",
  definition = function(x) {
    x@smooth
  }
)

#' @export
setGeneric("ms1<-", function(x, value) {
  standardGeneric("ms1<-")
})
setMethod(
  "ms1<-",
  signature = "MergedMSObject",
  definition = function(x, value) {
    x@ms1 <- value
    return(x)
  }
)

#' @export
setGeneric("ms2<-", function(x, value) {
  standardGeneric("ms2<-")
})
setMethod(
  "ms2<-",
  signature = "MergedMSObject",
  definition = function(x, value) {
    x@ms2 <- value
    return(x)
  }
)

#' @export
setGeneric("pre_iso_matched<-", function(x, value) {
  standardGeneric("pre_iso_matched<-")
})
setMethod(
  "pre_iso_matched<-",
  signature = "MergedMSObject",
  definition = function(x, value) {
    x@pre_iso_matched <- value
    return(x)
  }
)

#' @export
setGeneric("all_matched<-", function(x, value) {
  standardGeneric("all_matched<-")
})
setMethod(
  "all_matched<-",
  signature = "MergedMSObject",
  definition = function(x, value) {
    x@all_matched <- value
    return(x)
  }
)

#' @export
setGeneric("iso_matched<-", function(x, value) {
  standardGeneric("iso_matched<-")
})
setMethod(
  "iso_matched<-",
  signature = "MergedMSObject",
  definition = function(x, value) {
    x@iso_matched <- value
    return(x)
  }
)

#' @export
setGeneric("scaled_values<-", function(x, value) {
  standardGeneric("scaled_values<-")
})
setMethod(
  "scaled_values<-",
  signature = "MergedMSObject",
  definition = function(x, value) {
    x@scaled_values <- value
    return(x)
  }
)

#' @export
setGeneric("adjusted_df<-", function(x, value) {
  standardGeneric("adjusted_df<-")
})
setMethod(
  "adjusted_df<-",
  signature = "MergedMSObject",
  definition = function(x, value) {
    x@adjusted_df <- value
    return(x)
  }
)

#' @export
setGeneric("cutoffs<-", function(x, value) {
  standardGeneric("cutoffs<-")
})
setMethod(
  "cutoffs<-",
  signature = "MergedMSObject",
  definition = function(x, value) {
    x@cutoffs <- value
    return(x)
  }
)

#' @export
setGeneric("aligned<-", function(x, value) {
  standardGeneric("aligned<-")
})
setMethod(
  "aligned<-",
  signature = "MergedMSObject",
  definition = function(x, value) {
    x@aligned <- value
    return(x)
  }
)

#' @export
setGeneric("metadata<-", function(x, value) {
  standardGeneric("metadata<-")
})
setMethod(
  "metadata<-",
  signature = "MergedMSObject",
  definition = function(x, value) {
    x@metadata <- value
    return(x)
  }
)

#' @export
setGeneric("smooth<-", function(x, value) {
  standardGeneric("smooth<-")
})
setMethod(
  "smooth<-",
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
                          annot = FALSE,
                          id_name = "Compound_ID",
                          rt_name = "RT",
                          mz_name = "MZ",
                          int_name = "Intensity",
                          metab_name = NULL,
                          has_metadata = FALSE) {
  if (annot == TRUE && is.null(metab_name)) {
    stop("If annot is TRUE, then metabolite name must be supplied")
  }
  ms <- new("MSObject")
  name(ms) <- name
  consolidated(ms) <- FALSE

  raw_data <- df |>
    dplyr::select(id_name, metab_name, rt_name, mz_name, int_name) |>
    dplyr::rename(
      Compound_ID = id_name,
      Metabolite = metab_name,
      RT = rt_name,
      MZ = mz_name,
      Intensity = int_name
    )
  raw_df(ms) <- raw_data

  if (has_metadata) {
    meta_data <- df |>
      dplyr::select(id_name, -rt_name, -mz_name, -int_name, -metab_name)
    metadata(ms) <- meta_data
  }

  return(ms)
}

create_aligned_ms_obj <- function(ms1, ms2) {
  ms <- new(
    "MergedMSObject",
    ms1 = ms1,
    ms2 = ms2,
    all_matched = NULL,
    iso_matched = NULL,
    scaled_values = NULL,
    adjusted_df = NULL,
    cutoffs = NULL,
    aligned = NULL,
    metadata = NULL,
    smooth = NULL
  )
  return(ms)
}
