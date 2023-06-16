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
