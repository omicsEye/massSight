z_score <- function(x) {
  return((x - mean(x)) / stats::sd(x))
}

cvs <- function(x) {
  if (mean(x) == 0 & stats::sd(x) == 0) {
    return(0)
  } else {
    return(stats::sd(x) / mean(x))
  }
}

calc_cvs <- function(data, pref_to_use) {
  if (pref_to_use == "all") {
    pref_columns <- rep(T, ncol(data))
  } else {
    pref_columns <- ncol(data) %in% pref_to_use
  }

  pref_data <- data[, pref_columns] |>
    apply(2, cvs)
  return(pref_data)
}

find_indx <- function(df, word = "Metabolite") {
  if (word %in% colnames(df)) {
    return(c(0, which(colnames(df) == word)))
  } else {
    return(which(df == word, arr.ind = T))
  }
}

verify_df <- function(is_to_use, normalization, data, sample_information, pool_missing_p) {
  error <- character()
  error_flag <- FALSE

  if (!(0.0 <= pool_missing_p && pool_missing_p <= 100.0)) {
    error <- c(error, sprintf("Your pool missingness value must range from 0 to 100. You selected %s\n", pool_missing_p))
  }

  if ("IS" %in% normalization && is_to_use == "") {
    error <- c(error, "You have selected IS normalization but have not specified any internal standards. Either specify a standard or uncheck the IS box in the previous window.\n")
  }

  # Check column labels
  required_columns <- c("Compound_ID", "MZ", "RT", "Metabolite")
  for (col in required_columns) {
    if (!(col %in% colnames(data))) {
      error <- c(error, sprintf("There must be '%s' in your columns.\n", col))
      error_flag <- TRUE
    }
  }

  required_sample_info_columns <- c("Broad_name", "Collaborator_ID", "Injection_order", "Ref_to_use")
  for (col in required_sample_info_columns) {
    if (!(col %in% colnames(sample_information))) {
      error <- c(error, sprintf("There must be '%s' in your sample info columns.\n", col))
      error_flag <- TRUE
    }
  }

  if (error_flag) {
    return(paste(error, collapse = ""))
  }

  # Check sample information compared to sample labels
  metabolite_col <- which(colnames(data) == "Metabolite")
  for (i in seq_along(data[, (metabolite_col + 1):ncol(data)])) {
    column_name <- colnames(data)[metabolite_col + i]
    sample_info <- sample_information$Broad_name[i]
    if (column_name != sample_info) {
      error <- c(error, sprintf("There is at least one mismatch in your sample setup or sample name, starting with '%s' and '%s'.\n", column_name, sample_info))
      break
    }
  }

  # Check if internal standards are present
  for (IS in is_to_use) {
    if (!(IS %in% data$Metabolite)) {
      error <- c(error, sprintf("Your internal standard '%s' was not found in the 'Metabolite' column.\n", IS))
    }
  }

  # Check if there are prefs
  prefa_count <- sum(grepl("PREFA", sample_information$Collaborator_ID))
  prefb_count <- sum(grepl("PREFB", sample_information$Collaborator_ID))
  if (prefa_count == 0 && prefb_count == 0) {
    error <- c(error, sprintf("PREFs are missing. I counted %s PREFAs and %s PREFBs. There must be at least one pool reference.\n", prefa_count, prefb_count))
  }

  # Check if 'Ref_to_use' pools are actually present in the Broad_name name
  for (pref in unique(sample_information$Ref_to_use)) {
    if (!is.na(pref) && !(pref %in% sample_information$Collaborator_ID)) {
      error <- c(error, sprintf("Your 'Ref_to_use' '%s' was not found in the 'Collaborator_ID' column.\n", pref))
    }
  }

  return(paste(error, collapse = ""))
}

get_normalization_indices <- function(sample_injection, pool_injection, pool_names, ref_to_use) {
  pool_to_use_indices <- integer()

  pool_names_set <- unique(pool_names)

  for (i in seq_along(sample_injection)) {
    injection <- sample_injection[i]
    ref_to_use_single <- ref_to_use[i]

    if (!is.na(ref_to_use_single) && ref_to_use_single %in% pool_names_set) {
      pool_to_use_indices <- c(pool_to_use_indices, which(pool_names == ref_to_use_single)[1])
    } else {
      tryCatch(
        {
          closest_pool_index <- which.min(abs(pool_injection - injection))
          pool_to_use_indices <- c(pool_to_use_indices, closest_pool_index)
        },
        error = function(e) {
          return(NULL)
        }
      )
    }
  }

  return(pool_to_use_indices)
}

check_prefs <- function(sample_info,
                        pref_to_use,
                        prefs_to_remove) {
  sample_info |>
    dplyr::mutate(Platform_name = dplyr::case_when(
      Platfrom_name %in% prefs_to_remove ~ "do not use",
      T ~ Platfrom_name
    ))
  prefs_information <- sample_info |>
    dplyr::filter(stringr::str_detect(.data$Collaborator_ID, pref_to_use))
  # TODO add error
  return(prefs_information)
}

numeric_dataframe <- function(input) {
  input[, c(1:ncol(input))] <-
    sapply(sapply(input[, c(1:ncol(input))], as.character), as.numeric)
  return(input)
}
