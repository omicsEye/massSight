validate_parameters <-
  function(iso_method,
           match_method,
           smooth_method,
           minimum_intensity) {
    checkmate::assert_choice(smooth_method, c("loess", "gam", "gaussian"), .var.name = "smooth_method")
    checkmate::assert_choice(iso_method, c("manual", "dbscan"), .var.name = "iso_method")
    checkmate::assert_choice(match_method, c("supervised", "unsupervised"), .var.name = "match_method")
    checkmate::assert_numeric(minimum_intensity)
  }

log_parameters <- function(log, log_params, log_r, log_date, align_obj, ms1, ms2, time_start) {
  final_results <- all_matched(align_obj)
  compound_id_1 <- paste("Compound_ID", name(ms1), sep = "_")
  compound_id_2 <- paste("Compound_ID", name(ms2), sep = "_")
  n1 <- final_results |>
    dplyr::filter(!is.na(!!rlang::sym(compound_id_1)) &
      is.na(!!rlang::sym(compound_id_2))) |>
    nrow()
  n2 <- final_results |>
    dplyr::filter(is.na(!!rlang::sym(compound_id_1)) &
      !is.na(!!rlang::sym(compound_id_2))) |>
    nrow()
  n12 <- final_results |>
    dplyr::filter(!is.na(!!rlang::sym(compound_id_1)) &
      !is.na(!!rlang::sym(compound_id_2))) |>
    nrow()

  metabolite_1 <- paste("Metabolite", name(ms1), sep = "_")
  metabolite_2 <- paste("Metabolite", name(ms2), sep = "_")
  m_correct <- final_results |>
    dplyr::filter(!!rlang::sym(metabolite_1) == !!rlang::sym(metabolite_2)) |>
    nrow()
  m_total <- final_results |>
    dplyr::filter(!is.na(!!rlang::sym(metabolite_1)) &
      !is.na(!!rlang::sym(metabolite_2))) |>
    nrow()

  log_results <- list(
    "Dataset 1 Unmatched" = n1,
    "Dataset 2 Unmatched" = n2,
    "Matched" = n12,
    "Correct Matched Annotated Metabolites" = m_correct,
    "Total Matched Annotated Metabolites" = m_total,
    "Percentage Correct Matched Annotated Metabolites" = m_correct / m_total
  )

  time_end <- Sys.time()

  log_file <- list(
    date = log_date,
    r_version = log_r,
    parameters = log_params,
    results = log_results,
    runtime = as.numeric(time_end - time_start)
  )
  log_file <-
    jsonlite::toJSON(log_file, auto_unbox = TRUE, pretty = TRUE)

  writeLines(log_file, log)
}
