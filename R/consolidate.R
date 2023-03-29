#' @title Consolidate
#' @description Consolidate adducts/isomers with similar features into one
#' observation
#'
#' @param ms_obj a `massSight` object
#' @param rt_threshold user defined retention time threshold for defining
#' adducts/isomers
#' @param mz_threshold user defined mass to charge ratio threshold for defining
#' adducts/isomers
#' @param use_rt boolean if retention time should be used in identifying
#' metabolites that should be consolidated
#' @param use_mz boolean if mass to charge ratio should be used in identifying
#' metabolites that should be consolidated
#'
#' @return a `massSight` object with consolidated metabolites
#' @export
consolidate <-
  function(ms_obj,
           use_rt = T,
           use_mz = T,
           rt_threshold = NULL,
           mz_threshold = NULL) {
    stopifnot(!consolidated(ms), "Data is already consolidated")
    df <- ms_obj@raw
    if (use_rt) {
      stopifnot(
        !is.null(rt_threshold),
        "If `use_rt` is true, `rt_threshold` must be defined"
      )
      df <- df |> arrange(rt)
      i <- 1
      rt_adducts <- c()
      while (i <= nrow(df)) {
        j <- i + 1
        while (df[j, "rt"] - df[i, "rt"] <= rt_threshold) {
          rt_adducts <- rt_adducts |> c(df[j, "id"])
          j <- j + 1
          if (j > nrow(df)) {
            break
          }
        }
        i <- j
      }
    }
    if (use_mz) {
      stopifnot(
        !is.null(mz_threshold),
        "If `use_mz` is true, `mz_threshold` must be defined"
      )
      df <- df |> arrange(mz)
      i <- 1
      mz_adducts <- c()
      while (i < nrow(df)) {
        j <- i + 1
        while (df[j, "mz"] - df[i, "mz"] <= mz_threshold) {
          mz_adducts <- mz_adducts |> c(df[j, "id"])
          j <- j + 1
          if (j > nrow(df)) {
            break
          }
        }
        i <- j
      }
    }
    if (use_mz & use_rt) {
      df <- df |>
        dplyr::filter(!(id %in% intersect(rt_adducts, mz_adducts)))
    } else if (use_rt) {
      df <- df |>
        dplyr::filter(!(id %in% rt_adducts))
    } else if (use_mz) {
      df <- df |>
        dplyr::filter(!(id %in% mz_adducts))
    }
    ms_obj@raw <- df
    ms_obj@consolidated <- T

    if (nrow(ms_obj@metadata) != 0) {
      metadata <- ms_obj@metadata
      metadata <- metadata |>
        dplyr::filter(id %in% df$id)
      ms_obj@metadata <- metadata
    }

    return(ms_obj)
  }
