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
#' @param pool_missing_p A numeric value indicating the percentage of missing
#' pools allowed before skipping normalization (0-100).
#' @param fill_method A string indicating how to fill missing values. This can
#' be "half-min" or "none".
#' @param smooth_method A string indicating the smoothing method to use. This
#' can be "lowess", "line", "spline", or "gaussian".
#' @return A data frame of the normalized data.
auto_scale <- function(data,
                       sample_information,
                       is_to_use = "",
                       pref_to_use = "PREFA",
                       prefs_to_remove = "",
                       normalization = c("IS", "NN"),
                       pool_missing_p = 100.0,
                       fill_method = "none",
                       smooth_method = "lowess") {
  normalization <- toupper(normalization)
  error <- verify_df(is_to_use, normalization, data, sample_information, pool_missing_p)

  if (nchar(error) > 0) {
    return(error)
  }

  warnings <- ""

  metadata <- data[, c("Compound_ID", "MZ", "RT", "Metabolite")]
  metabolite_col <- which(colnames(data) == "Metabolite")
  data <- data[, (metabolite_col + 1):ncol(data)]
  data[data %in% c(0, 1, 2, 3, 4)] <- NA
  ndata <- data

  if ("IS" %in% normalization) {
    message("Starting IS...")

    compound_ids <- metadata$Metabolite
    is_index <- which(compound_ids %in% is_to_use)

    if (length(is_index) > 0) {
      scalars <- t(data[is_index, ])
      scalars <- scalars / apply(scalars, 2, median, na.rm = TRUE)
      scalars <- rowMeans(scalars, na.rm = TRUE)

      if (any(is.na(scalars))) {
        warnings <- paste0(warnings, "There are missing or zero values in your internal standards. ",
                           "The following samples were not IS normalized, and their IS values were filled with 1.0: ",
                           paste(colnames(data)[is.na(scalars)], collapse = ", "), "\n")
        scalars[is.na(scalars)] <- 1.0
      }

      ndata <- sweep(ndata, 2, scalars, "/")
      message("IS complete")
    } else {
      return("Please specify internal standards if you are using IS normalization.")
    }
  }

  if ("NN" %in% normalization) {
    message("Starting NN...")
    nn_result <- nn_normalize(ndata, sample_information, pref_to_use, prefs_to_remove, pool_missing_p)
    if (is.character(nn_result)) {
      return(nn_result)
    }
    ndata <- nn_result$normalized_data
    normalized <- nn_result$nn_normalized
    message("NN is complete.")
  } else if ("SMOOTH" %in% normalization) {
    message("Starting smooth normalization...")
    smooth_result <- smooth_normalize(ndata, sample_information, pref_to_use, prefs_to_remove, pool_missing_p, smooth_method)
    if (is.character(smooth_result)) {
      return(smooth_result)
    }
    ndata <- smooth_result$normalized_data
    normalized <- smooth_result$nn_normalized
    message("Smooth curve normalization is complete.")
  }

  if (fill_method == "half-min") {
    row_min <- apply(ndata, 1, min, na.rm = TRUE) / 2
    ndata[is.na(ndata)] <- rep(row_min, ncol(ndata))[is.na(ndata)]
  }

  ndata <- round(ndata, 0)

  raw_cvs <- calc_cvs(data, "all")
  raw_prefa <- calc_cvs(data, "PREFA")
  raw_prefb <- calc_cvs(data, "PREFB")

  colnames(ndata) <- colnames(data)
  n_cvs <- calc_cvs(ndata, "all")
  n_prefa <- calc_cvs(ndata, "PREFA")
  n_prefb <- calc_cvs(ndata, "PREFB")

  if ("NN" %in% normalization || "SMOOTH" %in% normalization) {
    ndata <- as.data.frame(ndata)
    ndata$Normalization_status <- normalized
  }

  ndata$Raw_CVs <- raw_cvs
  ndata$Raw_PREFA_CVs <- raw_prefa
  ndata$Raw_PREFB_CVs <- raw_prefb
  ndata$Normalized_CVs <- n_cvs
  ndata$Normalized_PREFA_CVs <- n_prefa
  ndata$Normalized_PREFB_CVs <- n_prefb

  ndata <- cbind(metadata, ndata)

  ndata$RT <- round(as.numeric(ndata$RT), 2)
  ndata$MZ <- round(as.numeric(ndata$MZ), 4)

  return(ndata)
}

smooth_normalize <- function(ndata, sample_information, pref_to_use, prefs_to_remove, pool_missing_p,
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
        message("Not smooth normalized, pool values:", paste(pool_values, collapse = ", "),
                " pool injections:", paste(pool_injections, collapse = ", "))
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

    normalization_scalars <- smooth_y / stats::median(smooth_y)
    normalized_data[i, ] <- as.numeric(ndata[i, ]) / normalization_scalars
    nn_normalized[i] <- "Smooth normalized"
  }

  return(list(normalized_data = normalized_data, nn_normalized = nn_normalized))
}

nn_normalize <- function(ndata, sample_information, pref_to_use, prefs_to_remove, pool_missing_p) {
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

z_score <- function(x) {
  return((x - mean(x)) / stats::sd(x))
}

cvs <- function(x) {
  temp_mean <- mean(x, na.rm = TRUE)
  temp_sd <- stats::sd(x, na.rm = TRUE)
  if (temp_mean == 0 && temp_sd == 0) {
    return(0.0)
  }
  return(temp_sd / temp_mean)
}

calc_cvs <- function(data, pref_to_use) {
  if (pref_to_use == "all") {
    pref_columns <- rep(TRUE, ncol(data))
  } else {
    pref_columns <- grepl(pref_to_use, colnames(data))
  }

  pref_data <- data[, pref_columns, drop = FALSE]
  cvs_result <- apply(pref_data, 2, cvs)
  cvs_result[is.na(cvs_result)] <- 0
  return(cvs_result)
}

find_indx <- function(df, word = "Metabolite") {
  if (word %in% colnames(df)) {
    return(c(0, which(colnames(df) == word)))
  } else {
    return(which(df == word, arr.ind = TRUE)[1, ])
  }
}

verify_df <- function(is_to_use, normalization, data, sample_information, pool_missing_p) {
  error <- character()

  if (!(0.0 <= pool_missing_p && pool_missing_p <= 100.0)) {
    error <- c(error, sprintf("Your pool missingness value must range from 0 to 100. You selected %s\n", pool_missing_p))
  }

  if ("IS" %in% normalization && is_to_use == "") {
    error <- c(error, "You have selected IS normalization but have not specified any internal standards. Either specify a standard or uncheck the IS box in the previous window.\n")
  }

  required_columns <- c("Compound_ID", "MZ", "RT", "Metabolite")
  for (col in required_columns) {
    if (!(col %in% colnames(data))) {
      error <- c(error, sprintf("There must be '%s' in your columns.\n", col))
    }
  }

  required_sample_info_columns <- c("Broad_name", "Collaborator_ID", "Injection_order", "Ref_to_use")
  for (col in required_sample_info_columns) {
    if (!(col %in% colnames(sample_information))) {
      error <- c(error, sprintf("There must be '%s' in your sample info columns.\n", col))
    }
  }

  metabolite_col <- which(colnames(data) == "Metabolite")
  for (i in seq_along(colnames(data)[(metabolite_col + 1):ncol(data)])) {
    column_name <- colnames(data)[metabolite_col + i]
    sample_info <- sample_information$Broad_name[i]
    if (column_name != sample_info) {
      error <- c(error, sprintf("There is at least one mismatch in your sample setup or sample name, starting with '%s' and '%s'.\n", column_name, sample_info))
      break
    }
  }

  for (IS in is_to_use) {
    if (!(IS %in% data$Metabolite)) {
      error <- c(error, sprintf("Your internal standard '%s' was not found in the 'Metabolite' column.\n", IS))
    }
  }

  prefa_count <- sum(grepl("PREFA", sample_information$Collaborator_ID))
  prefb_count <- sum(grepl("PREFB", sample_information$Collaborator_ID))
  if (prefa_count == 0 && prefb_count == 0) {
    error <- c(error, sprintf("PREFs are missing. I counted %s PREFAs and %s PREFBs. There must be at least one pool reference.\n", prefa_count, prefb_count))
  }

  for (pref in unique(sample_information$Ref_to_use)) {
    if (!is.na(pref) && !(pref %in% sample_information$Collaborator_ID)) {
      error <- c(error, sprintf("Your 'Ref_to_use' '%s' was not found in the 'Collaborator_ID' column.\n", pref))
    }
  }

  # Check for duplicate Injection_order
  if (any(duplicated(sample_information$Injection_order))) {
    error <- c(error, "There are duplicate Injection_order values in your sample information.\n")
  }

  # Check for duplicate PREF names
  pref_names <- sample_information$Collaborator_ID[grep("PREF", sample_information$Collaborator_ID)]
  if (any(duplicated(pref_names))) {
    error <- c(error, "There are duplicate PREF names in your Collaborator_ID column.\n")
  }

  return(paste(error, collapse = ""))
}

get_normalization_indices <- function(sample_injection, pool_injection, pool_names, ref_to_use) {
  pool_to_use_indices <- integer(length(sample_injection))

  pool_names_set <- unique(pool_names)

  for (i in seq_along(sample_injection)) {
    injection <- sample_injection[i]
    ref_to_use_single <- ref_to_use[i]

    if (!is.na(ref_to_use_single) && ref_to_use_single %in% pool_names_set) {
      pool_to_use_indices[i] <- which(pool_names == ref_to_use_single)[1]
    } else {
      closest_pool_index <- which.min(abs(pool_injection - injection))
      pool_to_use_indices[i] <- closest_pool_index
    }
  }

  return(pool_to_use_indices)
}

check_prefs <- function(sample_info, pref_to_use, prefs_to_remove) {
  for (pool in prefs_to_remove) {
    sample_info$Broad_name[grepl(pool, sample_info$Broad_name)] <- "do not use"
  }
  prefs_information <- sample_info[grepl(pref_to_use, sample_info$Collaborator_ID), ]
  if (nrow(prefs_information) == 0) {
    stop("No matching PREFs found for the given pref_to_use.")
  }
  return(prefs_information)
}

numeric_dataframe <- function(input) {
  input[, c(1:ncol(input))] <- sapply(input[, c(1:ncol(input))], as.numeric)
  return(input)
}

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

