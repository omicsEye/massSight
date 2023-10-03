#' @export
#' @title Auto Scale
#' @description This function will automatically scale your data based on the
#'  normalization method you choose. It will also calculate the CVs for each
#' sample and each metabolite.
#' @param is_to_use A vector of the internal standards to use for normalization.
#' If you do not want to use internal standards, leave this blank.
#' @param pref_to_use A string of the preferred reference to use for
#' normalization. If you do not want to use a preferred reference, leave this
#' blank.
#' @param prefs_to_remove A vector of the preferred references to remove from
#' the data. If you do not want to remove any preferred references, leave this
#' blank.
#' @param normalization A string of the normalization method to use. This can be
#' "IS", "NN", or "SMOOTH".
#' @param data A data frame of the data to be normalized. This should be the
#' output of the \code{read_data} function.
#' @param sample_information A data frame of the sample information. This should
#' be the output of the \code{read_sample_information} function.
#' @param pool_missing_p A boolean indicating whether or not to pool missing
#' values. If TRUE, missing values will be pooled. If FALSE, missing values will
#' be filled with the half-minimum.
#' @param fill_method A string indicating how to fill missing values. This can
#' be "half-min" or "mean". If "half-min", missing values will be filled with
#' the half-minimum. If "mean", missing values will be filled with the mean.
#' @param smooth_method A string indicating the smoothing method to use. This
#' can be "lowess" or "loess". If "lowess", a local polynomial regression will
#' be used. If "loess", a local polynomial regression will be used.
#' @return A data frame of the normalized data.
auto_scale <-
  function(is_to_use,
           pref_to_use,
           prefs_to_remove,
           normalization,
           data,
           sample_information,
           pool_missing_p,
           fill_method,
           smooth_method = "lowess") {
    normalization <- toupper(normalization)
    error <-
      verify_df(is_to_use,
                normalization,
                data,
                sample_information,
                pool_missing_p)
    if (length(error) > 0) {
      return(error)
    }

    warnings <- c()
    metadata <- data |>
      dplyr::select(.data$Compound_ID, .data$MZ, .data$RT, .data$Metabolite)

    data <- data |>
      dplyr::select(-.data$MZ,-.data$RT,-.data$Metabolite)
    data <- data[data %in% c(0, 1, 2, 3, 4)] <- NA
    ndata <- data
    if ("IS" %in% normalization) {
      message("Starting IS...")
      if (length(is_to_use) > 0) {
        scalars <- data |>
          dplyr::filter(.data$Compound_ID %in% is_to_use) |>
          dplyr::select(-.data$Compound_ID) |>
          t()
        scalars <- scalars / apply(scalars, 2, stats::median)
        scalars <- colMeans(scalars)
        if (any(is.na(scalars))) {
          warnings <-
            c(
              warnings,
              paste(
                "There are missing or zero values in your internal standards.",
                "The following samples were not IS normalized,",
                "and their IS values were filled with the half-minimumm:"
              )
            )
          # TODO line 516
        }
        ndata <- ndata / scalars
        message("IS complete")
      } else {
        return("Please specify internal standards")
      }
    }
    if ("NN" %in% normalization) {
      message("Starting NN...")
      nn_normalize_res <- nn_normalize(ndata,
                                       sample_information,
                                       pref_to_use,
                                       prefs_to_remove,
                                       pool_missing_p)
      ndata <- nn_normalize_res[1]
      normalized <- nn_normalize_res[2]
      rm(nn_normalize_res)
      message("NN Complete")
    } else if ("SMOOTH" %in% normalization) {
      message("Starting smooth normalization...")
      # TODO
    } else {
      if (fill_method == "half-min") {
        row_mins <- ndata |>
          apply(1, na.rm = T, FUN = min)
        row_mins <- row_mins / 2
        na_inds <- which(is.na(ndata), arr.ind = T)
        # TODO
      }
    }
    # calculate our CVs
    raw_cvs <- calc_cvs(data, "all")
    raw_prefa <- calc_cvs(data, "PREFA")
    raw_prefb <- calc_cvs(data, "PREFB")
    n_cvs <- calc_cvs(ndata, "all")
    n_prefa <- calc_cvs(ndata, "PREFA")
    n_prefb <- calc_cvs(ndata, "PREFB")
    if ("NN" %in% normalization | "SMOOTH" %in% normalization) {
      ndata$normalization_status <- normalized
    }
    ndata <- ndata |>
      dplyr::mutate(
        "raw_cvs" = raw_cvs,
        "raw_prefa" = raw_prefa,
        "raw_prefb" <- raw_prefb,
        "normalized_cvs" <- n_cvs,
        "normalized_prefa" <- n_prefa,
        "normalized_prefb" <- n_prefb
      )
    ndata <- rbind(metadata, ndata)
    ndata <- ndata |>
      dplyr::mutate(RT = round(.data$RT, 2),
                    MZ = round(.data$MZ, 4))
    return(ndata)
  }

get_normalization_ind <- function(sample_injections,
                                  pool_injections,
                                  pool_names,
                                  ref_to_use) {
  pool_names_set <- unique(pool_names)

  pool_bool <-
    (!is.null(ref_to_use)) & (ref_to_use %in% pool_names_set)
  pool_to_use <- pool_names[pool_bool]
  # TODO
  return(pool_to_use)
}
