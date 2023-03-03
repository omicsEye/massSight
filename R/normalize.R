nn_normalize <- function(ndata,
                         sample_information,
                         pref_to_use,
                         prefs_to_remove,
                         pool_missing_p) {
  skipped <- 0

  # TODO add missing value case for sample_info

  sample_information |>
    dplyr::mutate(Broad_name = case_when(Broad_name %in% prefs_to_remove ~ "do not use",
                                         T ~ Broad_name))
  prefs_information <- sample_information |>
    filter(stringr::str_detect(Collaborator_ID, pref_to_use))

  prefs_raw_data <- ndata[, colnames(prefs_information)]
  prefs_present <- !is.na(prefs_raw_data)
  for (i in 1:nrow(ndata)) {
    key <- prefs_present[i,]
    # TODO pool_names <-
  }
}

smooth_normalize <- function(ndata,
                             sample_information,
                             pref_to_use,
                             prefs_to_remove,
                             pool_missing_p) {
  prefs_information <-
    check_prefs(sample_information, pref_to_use, prefs_to_remove)
  prefs_raw_data <- ndata |>
    dplyr::select(prefs_information[, 1])
  prefs_present <- !is.na(prefs_raw_data)

  prefs_information <- prefs_information |>
    dplyr::arrange(Injection_order)
  sample_injection_order <- sample_information$Injection_order
  ref_to_use <- sample_information$Ref_to_use

  normalized_data <- data.frame()
  nn_normalized <- c()

  for (i in 1:nrow(ndata)) {
    key <- prefs_present[i,]
    pool_names <- prefs_information$Collaborator_ID[key]
    pool_injections <- prefs_information$Injection_order[key]
    pool_values <- prefs_raw_data[i, key]

    if ((1 - length(pool_names) / length(all_p_names)) * 100 >= pool_missing_p) {
      normalization_scalars <- ""
      skipped <- skipped + 1
    } else {
      pools_to_use_ind <-
        get_norm_indices(sample_injection_order,
                         pool_injections,
                         pool_names,
                         ref_to_use)
    }
    if (length(pool_values) < 2 | length(pool_injections) < 2) {
      if (length(pool_values) == 0 | length(pool_injections) == 0) {
        normalized_data <- normalized_data |>
          rbind(ndata[i,])
        nn_normalized <- nn_normalized |>
          c("Not smooth normalized")
        message(
          paste(
            "Not smooth normalized, pool values:",
            pool_values,
            "pool injections:",
            pool_injections
          )
        )
        next
      } else {
        normalization_scalars <-
          pool_values[pool_to_use_ind] / median(pool_values[unique(pool_to_use_ind)])
        normalized_data <- normalized_data |>
          rbind(ndata[i,] / normalization_scalars)
        nn_normalized <- nn_normalized |>
          c("NN Normalized")
        next
      }
    }
    # if (smooth)
  }
}
