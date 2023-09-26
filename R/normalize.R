nn_normalize <- function(ndata,
                         sample_information,
                         pref_to_use,
                         prefs_to_remove,
                         pool_missing_p) {
  skipped <- 0

  # TODO add missing value case for sample_info
  for (pool in prefs_to_remove) {
    message(paste("Removing pool reference:", pool))
    sample_information <- sample_information |>
      dplyr::mutate(name = dplyr::case_when(
        grepl(pool, name) ~ "do not use",
        TRUE ~ name
      ))
  }

  prefs_information <- sample_information |>
    filter(grepl(.data$Collaborator_ID, pref_to_use))

  prefs_raw_data <- ndata[, colnames(prefs_information)]
  prefs_present <- !is.na(prefs_raw_data)

  for (i in seq_len(nrow(ndata))) {
    key <- prefs_present[i, ]
    pool_names <- all_p_names[key]
    pool_injections <- prefs_information$Injection_order[key]
    pool_values <- prefs_raw_data[i, key]
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
    dplyr::arrange(.data$Injection_order)
  sample_injection_order <- sample_information$Injection_order
  ref_to_use <- sample_information$Ref_to_use

  normalized_data <- data.frame()
  nn_normalized <- c()

  for (i in seq_len(nrow(ndata))) {
    key <- prefs_present[i, ]
    pool_names <- prefs_information$Collaborator_ID[key]
    pool_injections <- prefs_information$Injection_order[key]
    pool_values <- prefs_raw_data[i, key]

    if ((1 - length(pool_names) / length(prefs_information$Collaborator_ID)) *
      100 >= pool_missing_p) {
      normalization_scalars <- ""
      skipped <- skipped + 1
    } else {
      pools_to_use_ind <-
        get_norm_indices(
          sample_injection_order,
          pool_injections,
          pool_names,
          ref_to_use
        )
    }
    if (length(pool_values) < 2 | length(pool_injections) < 2) {
      if (length(pool_values) == 0 | length(pool_injections) == 0) {
        normalized_data <- normalized_data |>
          rbind(ndata[i, ])
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
          pool_values[pools_to_use_ind] /
            median(pool_values[unique(pools_to_use_ind)])
        normalized_data <- normalized_data |>
          rbind(ndata[i, ] / normalization_scalars)
        nn_normalized <- nn_normalized |>
          c("NN Normalized")
        next
      }
    }
    if (smooth_method == "line") {
      smooth_x <- pool_injections[unique(pools_to_use_ind)]
      smooth_y <- pool_values[unique(pools_to_use_ind)]
    } else if (smooth_method == "lowess") {
      lowess_fun <-
        stats::lowess(
          pool_injections[unique(pools_to_use_ind)],
          pool_values[unique(pools_to_use_ind)],
          .2
        )
      smooth_x <- lowess_fun$x
      smooth_y <- lowess_fun$y
    }
  }
}
