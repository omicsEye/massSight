#' @title Consolidate
#' @description Consolidate adducts/isomers with similar features into one
#' observation
#' @param massSight_obj a `massSight` object
#' @param rt_threshold user defined retention time threshold for defining
#' adducts/isomers
#' @param mz_threshold user defined mass to charge ratio threshold for defining
#' adducts/isomers
#' @param int_threshold user defined intensity for defining adducts/isomers
#'
#' @return a `massSight` object with consolidated metabolites
#' @export
#'
#' @examples
consolidate <-
  function(massSight_obj,
           rt_threshold,
           mz_threshold,
           int_threshold) {

  }
