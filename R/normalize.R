nn_normalize <- function(ndata,
                         sample_information,
                         pref_to_use,
                         prefs_to_remove,
                         pool_missing_p) {
  skipped <- 0

  prefs_information <- check_prefs(sample_information, pref_to_use, prefs_to_remove)

  prefs_raw_data <- ndata[, prefs_information$Collaborator_ID]
  prefs_present <- !is.na(prefs_raw_data)

  all_p_names <- prefs_information$Collaborator_ID
  all_p_injections <- prefs_information$Injection_order
  sample_injection_order <- sample_information$Injection_order
  ref_to_use <- sample_information$Ref_to_use

  normalized_data <- matrix(nrow = nrow(ndata), ncol = ncol(ndata))
  nn_normalized <- character(nrow(ndata))

  for (i in seq_len(nrow(ndata))) {
    key <- prefs_present[i, ]
    pool_names <- all_p_names[key]
    pool_injections <- all_p_injections[key]
    pool_values <- as.numeric(prefs_raw_data[i, key])

    if ((1 - length(pool_names) / length(all_p_names)) * 100 >= pool_missing_p) {
      normalized_data[i, ] <- as.numeric(ndata[i, ])
      nn_normalized[i] <- "Not NN normalized"
      skipped <- skipped + 1
      next
    }

    pools_to_use_indices <- get_normalization_indices(sample_injection_order, pool_injections, pool_names, ref_to_use)
    normalization_scalars <- pool_values[pools_to_use_indices] / stats::median(pool_values[unique(pools_to_use_indices)])

    normalized_data[i, ] <- as.numeric(ndata[i, ]) / normalization_scalars
    nn_normalized[i] <- "NN normalized"
  }

  return(list(normalized_data = normalized_data, nn_normalized = nn_normalized))
}

smooth_normalize <- function(ndata,
                             sample_information,
                             pref_to_use,
                             prefs_to_remove,
                             pool_missing_p,
                             smooth_method = "lowess") {
  skipped <- 0
  prefs_information <- check_prefs(sample_information, pref_to_use, prefs_to_remove)
  prefs_raw_data <- ndata[, prefs_information$Collaborator_ID]
  prefs_present <- !is.na(prefs_raw_data)

  prefs_information <- prefs_information[order(prefs_information$Injection_order), ]
  sample_injection_order <- sample_information$Injection_order
  ref_to_use <- sample_information$Ref_to_use

  normalized_data <- matrix(nrow = nrow(ndata), ncol = ncol(ndata))
  nn_normalized <- character(nrow(ndata))

  for (i in seq_len(nrow(ndata))) {
    key <- prefs_present[i, ]
    pool_names <- prefs_information$Collaborator_ID[key]
    pool_injections <- prefs_information$Injection_order[key]
    pool_values <- as.numeric(prefs_raw_data[i, key])

    if ((1 - length(pool_names) / length(prefs_information$Collaborator_ID)) * 100 >= pool_missing_p) {
      normalized_data[i, ] <- as.numeric(ndata[i, ])
      nn_normalized[i] <- "Not smooth normalized"
      skipped <- skipped + 1
      next
    }

    pools_to_use_ind <- get_normalization_indices(sample_injection_order, pool_injections, pool_names, ref_to_use)

    if (length(pool_values) < 2 || length(pool_injections) < 2) {
      if (length(pool_values) == 0 || length(pool_injections) == 0) {
        normalized_data[i, ] <- as.numeric(ndata[i, ])
        nn_normalized[i] <- "Not smooth normalized"
        message(paste("Not smooth normalized, pool values:", paste(pool_values, collapse = ", "),
                      "pool injections:", paste(pool_injections, collapse = ", ")))
        next
      } else {
        normalization_scalars <- pool_values[pools_to_use_ind] / stats::median(pool_values[unique(pools_to_use_ind)])
        normalized_data[i, ] <- as.numeric(ndata[i, ]) / normalization_scalars
        nn_normalized[i] <- "NN normalized"
        next
      }
    }

    smooth_x <- pool_injections[unique(pools_to_use_ind)]
    smooth_y <- pool_values[unique(pools_to_use_ind)]

    if (smooth_method == "line") {
      f <- stats::lm(smooth_y ~ smooth_x)
      smooth_y <- predict(f, data.frame(smooth_x = sample_injection_order))
    } else if (smooth_method == "lowess") {
      lowess_res <- stats::lowess(smooth_x, smooth_y, f = 0.2)
      f <- stats::approxfun(lowess_res$x, lowess_res$y, rule = 2)
      smooth_y <- f(sample_injection_order)
    } else if (smooth_method == "spline") {
      f <- stats::smooth.spline(smooth_x, smooth_y)
      smooth_y <- predict(f, sample_injection_order)$y
    } else if (smooth_method == "gaussian") {
      if (!requireNamespace("kernlab", quietly = TRUE)) {
        stop("Package \"kernlab\" needed for Gaussian smoothing. Please install it.", call. = FALSE)
      }
      gp <- kernlab::gausspr(smooth_x, smooth_y)
      smooth_y <- kernlab::predict(gp, matrix(sample_injection_order, ncol = 1))
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

    normalized_data[i, ] <- as.numeric(ndata[i, ]) / normalization_scalars
    nn_normalized[i] <- "Smooth normalized"
  }

  return(list(normalized_data = normalized_data, nn_normalized = nn_normalized))
}
