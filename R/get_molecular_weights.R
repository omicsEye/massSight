mz_adjust <- function(ms_object,
                      parallel = TRUE,
                      n_cores = NULL) {
  metabolites <- raw_df(ms_object) |>
    dplyr::filter(!is.na(.data$Metabolite)) |>
    dplyr::pull(.data$Metabolite) |>
    unique()

  molecular_weights <- get_molecular_weights(metabolites, parallel = parallel, n_cores = n_cores)
  df_to_merge <- data.frame(Metabolite = metabolites, Molecular_Weight = molecular_weights)
  raw_df(ms_object) <- raw_df(ms_object) |>
    dplyr::left_join(df_to_merge, by = "Metabolite")
  return(ms_object)
}



get_molecular_weights <- function(metabolites,
                                  parallel = TRUE,
                                  n_cores = NULL) {
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 1 # Use one less than the total number of cores to avoid overloading the system
  }

  molecular_weights <- pbapply::pblapply(metabolites, function(name) {
    query_chem(name)
  }, cl = n_cores)

  return(molecular_weights)
}

query_chem <- function(name) {
  compound_info <- webchem::get_cid(name, from = "name")

  if (!is.null(compound_info$cid)) {
    cid <- as.numeric(compound_info$cid[1])
    compound_details <- webchem::pc_prop(cid, properties = c("MolecularWeight"))

    if (is.list(compound_details) &&
        "MolecularWeight" %in% names(compound_details)) {
      return(compound_details$MolecularWeight)
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
}
