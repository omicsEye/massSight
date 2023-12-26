#' Class to represent mass spectroscopy data.
#' @slot name A character indicating the name of the experiment.
#' @slot raw_df A data frame containing the raw data.
#' @slot isolated A data frame containing the isolated data.
#' @slot scaled_df A data frame containing the scaled data.
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
    raw_df = "data.frame",
    isolated = "data.frame",
    scaled_df = "data.frame",
    consolidated = "logical",
    metadata = "data.frame"
  )
)

#' @export
setGeneric("name", function(x) {
  standardGeneric("name")
})
setMethod(
  "name",
  signature = "MSObject",
  definition = function(x) {
    x@name
  }
)

#' @export
setGeneric("raw_df", function(x) {
  standardGeneric("raw_df")
})
setMethod(
  "raw_df",
  signature = "MSObject",
  definition = function(x) {
    x@raw_df
  }
)

#' @export
setGeneric("isolated", function(x) {
  standardGeneric("isolated")
})
setMethod(
  "isolated",
  signature = "MSObject",
  definition = function(x) {
    x@isolated
  }
)

#' @export
setGeneric("scaled_df", function(x) {
  standardGeneric("scaled_df")
})
setMethod(
  "scaled_df",
  signature = "MSObject",
  definition = function(x) {
    x@scaled_df
  }
)

#' @export
setGeneric("consolidated", function(x) {
  standardGeneric("consolidated")
})
setMethod(
  "consolidated",
  signature = "MSObject",
  definition = function(x) {
    x@consolidated
  }
)

#' @export
setGeneric("metadata", function(x) {
  standardGeneric("metadata")
})
setMethod(
  "metadata",
  signature = "MSObject",
  definition = function(x) {
    x@metadata
  }
)

#' @export
setGeneric("name<-", function(x, value) {
  standardGeneric("name<-")
})
setMethod(
  "name<-",
  signature = "MSObject",
  definition = function(x, value) {
    x@name <- value
    return(x)
  }
)

#' @export
setGeneric("raw_df<-", function(x, value) {
  standardGeneric("raw_df<-")
})
setMethod(
  "raw_df<-",
  signature = "MSObject",
  definition = function(x, value) {
    x@raw_df <- value
    return(x)
  }
)

#' @export
setGeneric("isolated<-", function(x, value) {
  standardGeneric("isolated<-")
})
setMethod(
  "isolated<-",
  signature = "MSObject",
  definition = function(x, value) {
    x@isolated <- value
    return(x)
  }
)

#' @export
setGeneric("scaled_df<-", function(x, value) {
  standardGeneric("scaled_df<-")
})
setMethod(
  "scaled_df<-",
  signature = "MSObject",
  definition = function(x, value) {
    x@scaled_df <- value
    return(x)
  }
)

#' @export
setGeneric("consolidated<-", function(x, value) {
  standardGeneric("consolidated<-")
})
setMethod(
  "consolidated<-",
  signature = "MSObject",
  definition = function(x, value) {
    x@consolidated <- value
    return(x)
  }
)

#' @export
setGeneric("metadata<-", function(x, value) {
  standardGeneric("metadata<-")
})
setMethod(
  "metadata<-",
  signature = "MSObject",
  definition = function(x, value) {
    x@metadata <- value
    return(x)
  }
)

# MergedMSObject ----------------------------------------------------------


#' Class to represent merged mass spectroscopy data.
#'
#' @slot ms1 A character indicating the name of the experiment.
#' @slot ms2 A data frame containing the raw data.
#' @slot all_matched A data frame containing the scaled data.
#' @slot iso_matched A data frame containing the scaled data.
#' @slot metadata A logical indicating whether or not the data has
#' been consolidated.
#' @slot pre_iso_matched A data frame containing the matched pairs before
#' initial isolation
#' @slot scaled_values A vector of scaled values
#' @slot adjusted_df A data frame containing the drift corrected metabolites
#' @slot cutoffs A vector of cutoffs
#' @slot aligned A dataframe containing the final metabolite mathces
#' @slot smooth_method A string indicating the smooth_methoding method used.
#'
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
    smooth_method = "character"
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
setMethod(
  "metadata",
  signature = "MergedMSObject",
  definition = function(x) {
    x@metadata
  }
)

#' @export
setGeneric("smooth_method", function(x) {
  standardGeneric("smooth_method")
})
setMethod(
  "smooth_method",
  signature = "MergedMSObject",
  definition = function(x) {
    x@smooth_method
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
setMethod(
  "metadata<-",
  signature = "MergedMSObject",
  definition = function(x, value) {
    x@metadata <- value
    return(x)
  }
)

#' @export
setGeneric("smooth_method<-", function(x, value) {
  standardGeneric("smooth_method<-")
})
setMethod(
  "smooth_method<-",
  signature = "MergedMSObject",
  definition = function(x, value) {
    x@smooth_method <- value
    return(x)
  }
)
