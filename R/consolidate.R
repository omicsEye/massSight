#' @title Consolidate
#' @description Consolidate adducts/isomers with similar features into one
#' observation
#'
#' @param massSight_obj a `massSight` object
#' @param rt_threshold user defined retention time threshold for defining
#' adducts/isomers
#' @param mz_threshold user defined mass to charge ratio threshold for defining
#' adducts/isomers
#' @param use_rt boolean if retention time should be used in identifying
#' metabolites that should be consolidated
#' @param use_mz boolean if mass to charge ration should be used in identifying
#' metabolites that should be consolidated
#'
#' @return a `massSight` object with consolidated metabolites
#' @export
#'
#' @examples
consolidate <-
  function(massSight_obj,
           use_rt = T,
           use_mz = T,
           rt_threshold = NULL,
           mz_threshold = NULL) {
    df <- massSight_obj@raw
    if (use_rt) {
      stopifnot(!is.null(rt_threshold),
                "If `use_rt` is true, `rt_threshold` must be defined")
      df <- df |> arrange(RT)
      i <- 1
      rt_adducts <- c()
      while (i <= nrow(df)) {
        j <- i + 1
        while (df[j, "RT"] - df[i, "RT"] <= rt_threshold) {
          rt_adducts <- rt_adducts |> c(j)
          j <- j + 1
          if (j > nrow(df)) {
            break
          }
        }
        i <- j
      }
    }
    if (use_mz) {
      stopifnot(!is.null(mz_threshold),
                "If `use_mz` is true, `mz_threshold` must be defined")
      df <- df |> arrange(MZ)
      i <- 1
      mz_adducts <- c()
      while (i < nrow(df)) {
        j <- i + 1
        while (df[j, "MZ"] - df[i, "MZ"] <= mz_threshold) {
          mz_adducts <- mz_adducts |> c(j)
          j <- j + 1
          if (j > nrow(df)) {
            break
          }
        }
        i <- j
      }
    }
    if (use_mz & use_rt) {
      df <- df[-(intersect(rt_adducts, mz_adducts))]
    } else if (use_rt) {
      df <- df[-rt_adducts]
    } else if (use_mz) {
      df <- df[-mz_adducts]
    }
    massSight_obj@consolidated <- df
    return(massSight_obj)
  }
