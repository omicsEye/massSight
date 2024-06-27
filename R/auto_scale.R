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
  data <- data[, !(names(data) %in% c("MZ", "RT", "Metabolite"))]
  data[data %in% c(0, 1, 2, 3, 4)] <- NA
  ndata <- data

  if ("IS" %in% normalization) {
    message("Starting IS...")

    compound_ids <- metadata$Metabolite
    is_index <- which(compound_ids %in% is_to_use)

    if (length(is_index) > 0) {
      scalars <- data[is_index, ]
      scalars <- t(scalars / apply(scalars, 2, stats::median))
      scalars <- rowMeans(scalars, na.rm = TRUE)

      if (any(is.na(scalars))) {
        warnings <- paste0(
          warnings, "There are missing or zero values in your internal standards. ",
          "The following samples were not IS normalized, and their IS values were filled with the half-minimum:"
        )
      }

      ndata <- ndata / scalars
      message("IS complete")
    } else {
      return("Please specify an internal standards if you are using IS normalization.")
    }
  }

  if ("NN" %in% normalization) {
    message("Starting NN...")

    nn_normalize_res <- nn_normalize(ndata, sample_information, pref_to_use, prefs_to_remove, pool_missing_p)
    ndata <- nn_normalize_res[[1]]
    normalized <- nn_normalize_res[[2]]

    message("NN is complete.")

    if (is.character(ndata)) {
      return(ndata)
    }

    ndata <- matrix(unlist(ndata), ncol = ncol(data))

    if (fill_method == "half-min") {
      row_min <- apply(ndata, 1, min, na.rm = TRUE) / 2
      inds <- which(is.na(ndata), arr.ind = TRUE)
      ndata[inds] <- row_min[inds[, 1]]
    }

    ndata <- as.data.frame(ndata)
    colnames(ndata) <- colnames(data)
    rownames(ndata) <- rownames(data)
    ndata <- round(ndata, 0)
  } else if ("SMOOTH" %in% normalization) {
    message("Starting smooth normalization...")

    smooth_normalize_res <- smooth_normalize(ndata, sample_information, pref_to_use, prefs_to_remove, pool_missing_p, smooth_method)
    ndata <- smooth_normalize_res[[1]]
    normalized <- smooth_normalize_res[[2]]

    message("Smooth curve normalization is complete.")

    if (is.character(ndata)) {
      return(ndata)
    }

    ndata <- matrix(unlist(ndata), ncol = ncol(data))

    if (fill_method == "half-min") {
      row_min <- apply(ndata, 1, min, na.rm = TRUE) / 2
      inds <- which(is.na(ndata), arr.ind = TRUE)
      ndata[inds] <- row_min[inds[, 1]]
    }

    ndata <- as.data.frame(ndata)
    colnames(ndata) <- colnames(data)
    rownames(ndata) <- rownames(data)
    ndata <- round(ndata, 0)
  } else {
    ndata <- as.matrix(ndata)

    if (fill_method == "half-min") {
      row_min <- apply(ndata, 1, min, na.rm = TRUE) / 2
      inds <- which(is.na(ndata), arr.ind = TRUE)
      ndata[inds] <- row_min[inds[, 1]]
    }

    ndata <- round(ndata, 0)
    ndata <- as.data.frame(ndata)
  }

  raw_cvs <- calc_cvs(data, "all")
  raw_prefa <- calc_cvs(data, "PREFA")
  raw_prefb <- calc_cvs(data, "PREFB")
  n_cvs <- calc_cvs(ndata, "all")
  n_prefa <- calc_cvs(ndata, "PREFA")
  n_prefb <- calc_cvs(ndata, "PREFB")

  if ("NN" %in% normalization) {
    ndata$Normalization_status <- normalized
  } else if ("SMOOTH" %in% normalization) {
    ndata$Normalization_status <- normalized
  }

  ndata$Raw_CVs <- raw_cvs
  ndata$Raw_PREFA_CVs <- raw_prefa
  ndata$Raw_PREFB_CVs <- raw_prefb
  ndata$Normalized_CVs <- n_cvs
  ndata$Normalized_PREFA_CVs <- n_prefa
  ndata$Normalized_PREFB_CVs <- n_prefb

  ndata <- cbind(metadata, ndata)

  ndata$RT <- round(ndata$RT, 2)
  ndata$MZ <- round(ndata$MZ, 4)

  return(ndata)
}

z_score <- function(x) {
  return((x - mean(x)) / stats::sd(x))
}

cvs <- function(x) {
  temp_mean <- mean(x)
  temp_sd <- stats::sd(x)
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
  return(apply(pref_data, 1, cvs))
}

find_indx <- function(df, word = "Metabolite") {
  if (word %in% colnames(df)) {
    return(c(0, which(colnames(df) == word)))
  }
  a <- which(df == word, arr.ind = TRUE)
  b <- a[, 1]
  row <- b[1]
  column <- a[1, 2]
  return(c(row, column))
}

verify_df <- function(is_to_use, normalization, data, sample_information, pool_missing_p) {
  error <- ""
  error_flag <- FALSE

  if (!(pool_missing_p >= 0.0 && pool_missing_p <= 100.0)) {
    error <- paste0(error, sprintf("Your pool missingness value must range from 0 to 100. You selected %s<BR><BR>\n", pool_missing_p))
  }

  if ("IS" %in% normalization && is_to_use == "") {
    error <- paste0(
      error, "You have selected IS normalization but have not specified any internal standards. Either specify ",
      "a standard or uncheck the IS box in the previous window.<BR><BR>\n"
    )
  }

  if (!("Compound_ID" %in% colnames(data))) {
    error <- paste0(error, "There must be 'Compound_ID' in your columns.<br><BR>\n")
    error_flag <- TRUE
  }
  if (!("MZ" %in% colnames(data))) {
    error <- paste0(error, "There must be 'MZ' in your columns.<br><BR>\n")
    error_flag <- TRUE
  }
  if (!("RT" %in% colnames(data))) {
    error <- paste0(error, "There must be 'RT' in your columns.<br><BR>\n")
    error_flag <- TRUE
  }
  if (!("Metabolite" %in% colnames(data))) {
    error <- paste0(error, "There must be 'Metabolite' in your columns.<br><BR>\n")
    error_flag <- TRUE
  }
  if (!("Broad_name" %in% colnames(sample_information))) {
    error <- paste0(error, "There must be 'Broad_name' in your sample info columns.<br><BR>\n")
    error_flag <- TRUE
  }
  if (!("Collaborator_ID" %in% colnames(sample_information))) {
    error <- paste0(error, "There must be 'Collaborator_ID' in your sample info columns.<br><BR>\n")
    error_flag <- TRUE
  }
  if (!("Injection_order" %in% colnames(sample_information))) {
    error <- paste0(error, "There must be 'Injection_order' in your sample info columns.<br><BR>\n")
    error_flag <- TRUE
  }
  if (!("Ref_to_use" %in% colnames(sample_information))) {
    error <- paste0(error, "There must be 'Ref_to_use' in your sample info columns.<br><BR>\n")
    error_flag <- TRUE
  }

  if (error_flag) {
    return(error)
  }

  metabolite_row <- find_indx(data, word = "Metabolite")[1]
  metabolite_col <- find_indx(data, word = "Metabolite")[2]

  for (i in seq_along(colnames(data)[(metabolite_col + 1):ncol(data)])) {
    column_name <- colnames(data)[metabolite_col + i]
    sample_info <- sample_information$Broad_name[i]

    if (column_name != sample_info) {
      error <- paste0(error, sprintf(
        "There is at least one mismatch in your sample setup or sample name, starting the with '%s' and '%s'.<BR><BR>\n",
        column_name, sample_info
      ))
      break
    }
  }

  for (IS in is_to_use) {
    if (!(IS %in% data$Metabolite)) {
      error <- paste0(error, sprintf("Your internal standard '%s' was not found in the 'Metabolite' column.<BR>\n", IS))
    }
  }

  if (sum(grepl("PREFA", sample_information$Collaborator_ID)) == 0 &&
    sum(grepl("PREFB", sample_information$Collaborator_ID)) == 0) {
    error <- paste0(error, sprintf(
      "PREFs are missing. I counted %s PREFAs and %s PREFBs. There must be at least one pool reference.<BR><BR>\n",
      sum(grepl("PREFA", sample_information$Collaborator_ID)),
      sum(grepl("PREFB", sample_information$Collaborator_ID))
    ))
  }

  for (pref in unique(sample_information$Ref_to_use)) {
    if (!is.na(pref) && !(pref %in% sample_information$Collaborator_ID)) {
      error <- paste0(error, sprintf("Your 'Ref_to_use' '%s' was not found in the 'Collaborator_ID' column.<BR><BR>", pref))
    }
  }

  return(error)
}

get_normalization_indices <- function(sample_injection, pool_injection, pool_names, ref_to_use) {
  pool_to_use_indices <- c()

  pool_names_set <- unique(pool_names)

  for (i in seq_along(sample_injection)) {
    injection <- sample_injection[i]
    ref_to_use_single <- ref_to_use[i]

    if (!is.na(ref_to_use_single) && ref_to_use_single %in% pool_names_set) {
      pool_to_use_indices <- c(pool_to_use_indices, which(pool_names == ref_to_use_single)[1])
    } else {
      tryCatch(
        {
          pool_to_use_indices <- c(pool_to_use_indices, which.min(abs(pool_injection - injection)))
        },
        error = function(e) {
          return(NULL)
        }
      )
    }
  }

  return(pool_to_use_indices)
}

smooth_normalize <- function(ndata, sample_information, pref_to_use, prefs_to_remove, pool_missing_p,
                             smooth_method = "lowess") {
  skipped <- 0

  for (pool in prefs_to_remove) {
    message("Removing pool reference: ", pool)
    sample_information$Broad_name[grepl(pool, sample_information$Broad_name)] <- "do not use"
  }

  prefs_information <- sample_information[grepl(pref_to_use, sample_information$Collaborator_ID), ]

  prefs_raw_data <- ndata[, prefs_information[, 1]]
  prefs_present <- !is.na(prefs_raw_data)
  prefs_information <- prefs_information[order(prefs_information$Injection_order), ]
  all_p_names <- prefs_information$Collaborator_ID
  all_p_injections <- prefs_information$Injection_order
  sample_injection_order <- sample_information$Injection_order
  ref_to_use <- sample_information$Ref_to_use
  normalized_data <- list()
  nn_normalized <- c()
  for (i in seq_len(nrow(ndata))) {
    key <- prefs_present[i, ]
    pool_names <- all_p_names[key]
    pool_injections <- all_p_injections[key]
    pool_values <- as.numeric(prefs_raw_data[i, key])
    if ((1 - (length(pool_names) / length(all_p_names))) * 100 >= pool_missing_p) {
      normalization_scalars <- NULL
      skipped <- skipped + 1
    } else {
      pools_to_use_indices <- get_normalization_indices(
        sample_injection_order,
        pool_injections, pool_names, ref_to_use
      )
    }

    if (length(pool_values) < 2 || length(pool_injections) < 2) {
      if (length(pool_values) == 0 || length(pool_injections) == 0) {
        normalized_data[[i]] <- as.numeric(ndata[i, ])
        nn_normalized <- c(nn_normalized, "Not smooth normalized")
        message("Not smooth normalized", " pool values:", pool_values, " pool injections:", pool_injections)
        next
      } else {
        normalization_scalars <- pool_values[pools_to_use_indices] / stats::median(pool_values[unique(pools_to_use_indices)])
        normalized_data[[i]] <- as.numeric(ndata[i, ]) / normalization_scalars
        nn_normalized <- c(nn_normalized, "NN normalized")
        next
      }
    }

    if (smooth_method == "line") {
      smooth_x <- pool_injections[unique(pools_to_use_indices)]
      smooth_y <- pool_values[unique(pools_to_use_indices)]
    } else if (smooth_method == "lowess") {
      lowess_res <- stats::lowess(pool_values[unique(pools_to_use_indices)], pool_injections[unique(pools_to_use_indices)], f = 0.2, delta = 0.0, iter = 5)
      smooth_x <- lowess_res[, 1]
      smooth_y <- lowess_res[, 2]
    } else if (smooth_method == "spline") {
      smooth_x <- pool_injections
      smooth_y <- splines::spline(pool_injections[unique(pools_to_use_indices)], pool_values[unique(pools_to_use_indices)], xout = smooth_x)$y
    } else if (smooth_method == "gaussian") {
      smooth_x <- pool_injections[unique(pools_to_use_indices)]
      kernel <- 10 * stats::rbfdot(10)
      gp <- kernlab::gausspr(kernel, smooth_x, pool_values[unique(pools_to_use_indices)])
      smooth_y <- kernlab::predict(gp, matrix(smooth_x, ncol = 1))
    }

    smooth_y_dropna <- smooth_y[!is.na(smooth_y)]
    smooth_x_dropna <- smooth_x[!is.na(smooth_y)]

    m <- min(smooth_x_dropna)
    argmin <- which(smooth_x_dropna == m)
    smooth_min <- smooth_y_dropna[argmin[1]]

    m <- max(smooth_x_dropna)
    argmax <- which(smooth_x_dropna == m)
    smooth_max <- smooth_y_dropna[argmax[1]]

    f <- stats::approxfun(smooth_x_dropna, smooth_y_dropna, method = "linear", yleft = smooth_min, yright = smooth_max)

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
nn_normalize <- function(ndata, sample_information, pref_to_use, prefs_to_remove, pool_missing_p) {
  skipped <- 0
  for (pool in prefs_to_remove) {
    message("Removing pool reference: ", pool)
    sample_information$Broad_name[grepl(pool, sample_information$Broad_name)] <- "do not use"
  }
  prefs_information <- sample_information[grepl(pref_to_use, sample_information$Collaborator_ID), ]
  prefs_raw_data <- ndata[, prefs_information[, 1]]
  prefs_present <- !is.na(prefs_raw_data)
  ndata <- as.matrix(ndata)
  all_p_names <- prefs_information$Collaborator_ID
  all_p_injections <- prefs_information$Injection_order
  sample_injection_order <- sample_information$Injection_order
  ref_to_use <- sample_information$Ref_to_use
  normalized_data <- list()
  nn_normalized <- c()
  for (i in seq_len(nrow(ndata))) {
    key <- prefs_present[i, ]
    pool_names <- all_p_names[key]
    pool_injections <- all_p_injections[key]
    pool_values <- as.numeric(prefs_raw_data[i, key])

    if ((1 - (length(pool_names) / length(all_p_names))) * 100 >= pool_missing_p) {
      normalization_scalars <- NULL
      skipped <- skipped + 1
    } else {
      pools_to_use_indices <- get_normalization_indices(
        sample_injection_order,
        pool_injections, pool_names, ref_to_use
      )
      normalization_scalars <- pool_values[pools_to_use_indices] / stats::median(pool_values[unique(pools_to_use_indices)])
    }

    tryCatch(
      {
        normalized_data[[i]] <- as.numeric(ndata[i, ]) / normalization_scalars
        nn_normalized <- c(nn_normalized, "NN normalized")
      },
      error = function(e) {
        normalized_data[[i]] <- as.numeric(ndata[i, ])
        nn_normalized <- c(nn_normalized, "Not NN normalized")
      }
    )
  }
  return(list(
    normalized_data = do.call(rbind, normalized_data),
    nn_normalized = nn_normalized
  ))
}
