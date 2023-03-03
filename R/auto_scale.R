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
      dplyr::select("Compound_ID", "MZ", "RT", "Metabolite")

    data <- data |>
      dplyr::select(-"MZ", -"RT", -"Metabolite")
    data <- data[data %in% c(0, 1, 2, 3, 4)] <- NA
    ndata <- data
    if ("IS" %in% normalization) {
      message("Starting IS...")
      if (length(is_to_use) > 0) {
        scalars <- data |>
          dplyr::filter("Compound_ID" %in% is_to_use) |>
          dplyr::select(-"Compound_ID") |>
          t()
        scalars <- scalars / apply(scalars, 2, median)
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
      dplyr::mutate(RT = round("RT", 2),
                    MZ = round("MZ", 4))
    return(ndata)
  }

smooth_normalize <-
  function(ndata,
           sample_information,
           pref_to_use,
           prefs_to_remove,
           pool_missing_p,
           smooth_method = "lowess") {
    skipped <- 0

    # rename prefs we want to ignore
    prefs_to_remove %>%
      purrr::walk(function(pool) {
        sample_information <- sample_information %>%
          dplyr::mutate(Broad_name = dplyr::case_when(
            stringr::str_detect(Broad_name, pool) ~ "do not use",
            TRUE ~ Broad_name
          ))
        message(paste("Removing pool reference:", pool))
      })

    # look in the metadata, and shorten to wherever the short name
    # has PREFA or PREFB in it
    pref_info <- sample_information %>%
      filter(stringr::str_detect(Collaborator_ID, pref_to_use))
    if (!(stringr::str_detect(sample_information$Collaborator_ID,
                              pref_to_use))) {
      # TODO
    }
    pref_info <- pref_info %>%
      dplyr::arrange(Injection_order)
    all_p_names <- pref_info$Collaborator_ID
    all_p_injections <- pref_info$Injection_order
    sample_injection_order <- sample_information$Injection_order
    ref_to_use <- sample_information$Ref_to_use
  }
