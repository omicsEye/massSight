z_score <- function(x) {
  return((x - mean(x)) / stats::sd(x))
}

cvs <- function(x) {
  x <- x[!is.na(x)]
  temp_mean <- mean(x)
  temp_sd <- sd(x)
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
