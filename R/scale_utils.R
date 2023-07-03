z_score <- function(x) {
  return((x - mean(x)) / sd(x))
}

cvs <- function(x) {
  if (mean(x) == 0 & sd(x) == 0) {
    return(0)
  } else {
    return(sd(x) / mean(x))
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

verify_df <-
  function(data,
           sample_info,
           is_to_use,
           norm,
           pool_missing_p) {
    error <- c()
    error_flag <- FALSE
    if (!(0 <= pool_missing_p &
      pool_missing_p <= 100)) {
      error <- error %>%
        append(
          paste(
            "Your pool missingness value must range from 0 to 100. You selected",
            pool_missing_p
          )
        )
    }
    if ("IS" %in% norm & is_to_use == "") {
      error <- error %>%
        append(
          "You have selected IS normalization but have not specified any internal standards. Either specify a standard or uncheck the IS box in the previous window."
        )
    }
    needed_cols <-
      c(
        "Compound_ID",
        "MZ",
        "RT",
        "Metabolite",
        "Platfrom_name",
        "Collaborator_ID",
        "Injection_order",
        "Ref_to_use"
      )
    if (!(c(needed_cols %in% names(data)))) {
      for (col in needed_cols) {
        if (!(col %in% names(data))) {
          error <- error %>%
            append(paste("There must be", col, "in your columns."))
        }
      }
      error_flag <- TRUE
    }
    if (error_flag) {
      message(error)
    }
    # TODO py line 92

    is_not_found <-
      is_not_found[which(is_to_use %in% data$Metabolite)]
    if (length(is_not_found) > 0) {
      purrr::walk(is_not_found, function(x) {
        error <- error %>%
          append(paste(
            "Your internal standard",
            x,
            "was not found in the 'Metabolite' column"
          )) # TODO
      })
    }

    prefa_sum <- sample_info$Collaborator_ID %>%
      str_detect("PREFA") %>%
      sum()
    prefb_sum <- sample_info$Collaborator_ID %>%
      str_detect("PREFB") %>%
      sum()

    if (prefa_sum == 0 &
      prefb_sum == 0) {
      error <- error %>%
        append(
          paste(
            "PREFs are missing. I counted",
            pref_a,
            "PREFAs and",
            preb_b,
            "PREFBs. There must be at least one pool reference"
          )
        )
    }
    # check if 'Ref_to_use' pools are actually present in the Platfrom_name name
    bad_refs <-
      samples_info$Ref_to_use[which(!samples_info$Ref_to_use %in% sample_info$Collaborator_ID)]
    if (length(bad_refs) > 0) {
      error <- error %>%
        append(
          paste0(
            "Your 'Ref_to_use':,",
            bad_refs,
            "was not found in the 'Collaborator_ID' column"
          )
        )
    }
    return(error)
  }

get_norm_indices <-
  function(sample_inj, pool_inj, pool_names, ref_use) {
    pool_use_indx <- c()
    pool_names_unique <- pool_names %>% unique()
    sample_ref_tibble <- dplyr::tibble(
      sample_inj,
      ref_use
    )

    purrr::walk(sample_ref_tibble, function(row) {
      if (!(is.null(row$ref_use)) & row$ref_use %in% pool_names_unique) {
        pool_use_indx <- pool_use_indx %>%
          append(which(pool_names == row$ref_use)[1])
      } else {
        pool_use_indx <- pool_use_indx %>%
          append() # TODO
      }
    })
  }

check_prefs <- function(sample_info,
                        pref_to_use,
                        prefs_to_remove) {
  sample_information |>
    dplyr::mutate(Platfrom_name = case_when(
      Platfrom_name %in% prefs_to_remove ~ "do not use",
      T ~ Platfrom_name
    ))
  prefs_information <- sample_information |>
    filter(stringr::str_detect(Collaborator_ID, pref_to_use))
  # TODO add error
  return(prefs_information)
}

numeric_dataframe <- function(input) {
  input[, c(1:ncol(input))] <-
    sapply(sapply(input[, c(1:ncol(input))], as.character), as.numeric)
  return(input)
}
