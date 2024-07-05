nn_normalize <- function(ndata,
                         sample_information,
                         pref_to_use,
                         prefs_to_remove,
                         pool_missing_p) {
  skipped <- 0

  # TODO add missing value case for sample_info
  # for (pool in prefs_to_remove) {
  #   message(paste("Removing pool reference:", pool))
  #   browser()
  #   sample_information <- sample_information |>
  #     dplyr::mutate(name = dplyr::case_when(
  #       grepl(pool, name) ~ "do not use",
  #       TRUE ~ name
  #     ))
  # }

  prefs_information <- sample_information |>
    dplyr::filter(grepl(pref_to_use, .data$Collaborator_ID, fixed = TRUE))

  prefs_raw_data <- ndata[, prefs_information$Collaborator_ID]
  prefs_present <- !is.na(prefs_raw_data)

  for (i in seq_len(nrow(ndata))) {
    key <- prefs_present[i, ]
    # pool_names <- all_p_names[key]
    pool_injections <- prefs_information$Injection_order[key]
    pool_values <- prefs_raw_data[i, key]
  }
}

smooth_normalize <- function(ndata,
                             sample_information,
                             pref_to_use,
                             prefs_to_remove,
                             pool_missing_p,
                             smooth_method = "lowess") {
  skipped <- 0
  prefs_information <-
    check_prefs(sample_information, pref_to_use, prefs_to_remove)
  prefs_raw_data <- ndata |>
    dplyr::select(prefs_information[, 1])
  prefs_present <- !is.na(prefs_raw_data)

  prefs_information <- prefs_information |>
    dplyr::arrange(.data$Injection_order)
  sample_injection_order <- sample_information$Injection_order
  ref_to_use <- sample_information$Ref_to_use

  normalized_data <- list()
  nn_normalized <- c()

  for (i in seq_len(nrow(ndata))) {
    key <- prefs_present[i, ]
    pool_names <- prefs_information$Collaborator_ID[key]
    pool_injections <- prefs_information$Injection_order[key]
    pool_values <- as.numeric(prefs_raw_data[i, key])

    if ((1 - length(pool_names) / length(prefs_information$Collaborator_ID)) *
      100 >= pool_missing_p) {
      normalization_scalars <- NULL
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
        normalized_data[[i]] <- as.numeric(ndata[i, ])
        nn_normalized <- c(nn_normalized, "Not smooth normalized")
        message(
          paste(
            "Not smooth normalized, pool values:",
            paste(pool_values, collapse = ", "),
            "pool injections:",
            paste(pool_injections, collapse = ", ")
          )
        )
        next
      } else {
        normalization_scalars <-
          pool_values[pools_to_use_ind] /
            stats::median(pool_values[unique(pools_to_use_ind)])
        normalized_data[[i]] <- as.numeric(ndata[i, ]) / normalization_scalars
        nn_normalized <- c(nn_normalized, "NN normalized")
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
          f = 0.2
        )
      smooth_x <- lowess_fun$x
      smooth_y <- lowess_fun$y
    } else if (smooth_method == "spline") {
      spline_fun <- stats::smooth.spline(
        pool_injections[unique(pools_to_use_ind)],
        pool_values[unique(pools_to_use_ind)]
      )
      smooth_x <- pool_injections
      smooth_y <- stats::predict(spline_fun, smooth_x)$y
    } else if (smooth_method == "gaussian") {
      # Note: Gaussian Process Regression is not a standard function in R
      # You might need to use an external package like 'kernlab' for this
      message("Gaussian smoothing method not implemented in this R version")
      next
    }

    smooth_y_dropna <- smooth_y[!is.na(smooth_y)]
    smooth_x_dropna <- smooth_x[!is.na(smooth_y)]

    smooth_min <- smooth_y_dropna[which.min(smooth_x_dropna)]
    smooth_max <- smooth_y_dropna[which.max(smooth_x_dropna)]

    f <- stats::approxfun(smooth_x_dropna, smooth_y_dropna,
      method = "linear",
      yleft = smooth_min,
      yright = smooth_max
    )

    smooth_x <- sample_injection_order
    smooth_y <- f(smooth_x)
    normalization_scalars <- smooth_y / stats::median(smooth_y)

    normalized_data[[i]] <- as.numeric(ndata[i, ]) / normalization_scalars
    nn_normalized <- c(nn_normalized, "Smooth normalized")
  }

  return(list(
    normalized_data = do.call(rbind, normalized_data),
    nn_normalized = nn_normalized
  ))
}
