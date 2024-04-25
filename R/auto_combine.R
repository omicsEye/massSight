#' @export
#' @title Auto Combine
#' @description Combines two `massSight` objects, resulting in a single
#' `MergedMSObject`.
#'
#' @param ms1 A `massSight` object representing the results of a preprocessed
#' LC-MS experiment.
#' @param ms2 A `massSight` object representing the results of a second
#' preprocessed LC-MS experiment.
#' @param rt_lower A numeric indicating the lower bound of the RT
#' range to be considered for aligning two metabolites.
#' @param rt_upper A numeric indicating the upper bound of the RT
#' range to be considered for aligning two metabolites.
#' @param mz_lower A numeric indicating the lower bound of the m/z
#' range to be considered for aligning two metabolites.
#' @param mz_upper A numeric indicating the upper bound of the m/z
#' range to be considered for aligning two metabolites.
#' @param minimum_intensity A numeric indicating the minimum intensity
#' to be considered for alignment.
#' @param iso_method The isolation method used before modeling drift. Can
#' either be "manual" or "dbscan".
#' @param eps Epsilon value for dbscan algorithm. Only used if iso_method =
#' "dbscan"
#' @param rt_iso_threshold A numeric indicating the isolation
#' parameter for RT.
#' @param mz_iso_threshold A numeric indicating the isolation
#' parameter for m/z.
#' @param match_method A character indicating the initial matching method to
#' be used to detect inter-batch variability. Options are "unsupervised" and
#' "supervised".
#' @param smooth_method A character indicating the smoothing method to
#' be used. Options are "lowess", "spline", and "gaussian".
#' @param weights A numeric vector indicating the weights to be used for
#' the alignment.
#' @param keep_features A logical vector indicating whether or not to
#' @param log A character indicating the name of the log file.
#' @param output A character indicating the directory to save the output. If NULL,
#' the output will be saved in the current working directory.
#'
#' @return A `MergedMSObject` containing the combined data.
setGeneric("auto_combine", function(ms1,
                                    ms2,
                                    rt_lower = -.5,
                                    rt_upper = .5,
                                    mz_lower = -15,
                                    mz_upper = 15,
                                    minimum_intensity = 10,
                                    iso_method = "manual",
                                    eps = .1,
                                    rt_iso_threshold = .1,
                                    mz_iso_threshold = 5,
                                    match_method = "unsupervised",
                                    smooth_method = "gam",
                                    weights = c(1, 1, 1),
                                    keep_features = c(F, F),
                                    log = NULL,
                                    output = NULL) {
  standardGeneric("auto_combine")
})

setMethod(
  "auto_combine",
  signature("MSObject", "MSObject"),
  function(ms1,
           ms2,
           rt_lower = -.5,
           rt_upper = .5,
           mz_lower = -15,
           mz_upper = 15,
           minimum_intensity = 10,
           iso_method = "manual",
           eps = .1,
           rt_iso_threshold = .1,
           mz_iso_threshold = 5,
           match_method = "unsupervised",
           smooth_method = "gam",
           weights = c(1, 1, 1),
           keep_features = c(FALSE, FALSE),
           log = "log.json",
           output = NULL) {
    if (!is.null(log)) {
      time_start <- Sys.time()
      check_jsonlite()
      log_params <- match.call(expand.dots = TRUE) |>
        modify_call() |>
        as.list() |>
        (\(x) x[-1])() |>
        lapply(\(x) {
          if (is.language(x)) {
            return(deparse(x))
          } else {
            return(x)
          }
        })
      log_r <- R.Version()
      log_date <- Sys.time()
    }

    validate_parameters(iso_method, match_method, smooth_method, minimum_intensity)

    if (match_method == "unsupervised") {
      if (iso_method == "manual") {
        ref_iso <- getVectors(raw_df(ms1),
          rt_sim = rt_iso_threshold,
          mz_sim = mz_iso_threshold
        )
        query_iso <- getVectors(raw_df(ms2),
          rt_sim = rt_iso_threshold,
          mz_sim = mz_iso_threshold
        )
        isolated(ms1) <- raw_df(ms1) |>
          dplyr::filter(.data$Compound_ID %in% ref_iso)
        isolated(ms2) <- raw_df(ms2) |>
          dplyr::filter(.data$Compound_ID %in% query_iso)
      } else if (iso_method == "dbscan") {
        isolated(ms1) <- iso_dbscan(raw_df(ms1), eps)
        isolated(ms2) <- iso_dbscan(raw_df(ms2), eps)
      }
    } else if (match_method == "supervised") {
      isolated(ms1) <- raw_df(ms1) |>
        dplyr::filter(.data$Metabolite != "")
      isolated(ms2) <- raw_df(ms2) |>
        dplyr::filter(.data$Metabolite != "")
    }

    align_obj <- methods::new("MergedMSObject")
    ms1(align_obj) <- ms1
    ms2(align_obj) <- ms2
    align_obj <- align_obj |>
      # align_pre_isolated_compounds(
      #   rt_minus = rt_lower,
      #   rt_plus = rt_upper,
      #   mz_minus = mz_lower,
      #   mz_plus = mz_upper
      # ) |>
      align_isolated_compounds(
        match_method = match_method,
        rt_minus = rt_lower,
        rt_plus = rt_upper,
        mz_minus = mz_lower,
        mz_plus = mz_upper
      ) |>
      smooth_drift(
        smooth_method = smooth_method,
        minimum_int = minimum_intensity
      ) |>
      final_results(
        keep_features = keep_features,
        weights = weights
      )

    if (!is.null(log)) {
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
      time <- time_end - time_start
      date$runtime <- time

      log_file <- list(
        date = log_date,
        r_version = log_r,
        parameters = log_params,
        results = log_results
      )
      log_file <-
        jsonlite::toJSON(log_file, auto_unbox = TRUE, pretty = TRUE)

      writeLines(log_file, log)

      return(align_obj)
    }
    return(align_obj)
  }
)

setMethod(
  "auto_combine",
  signature("MergedMSObject", "MSObject"),
  function(ms1, ms2) {
    validate_parameters(
      iso_method,
      match_method,
      smooth_method,
      minimum_intensity
    )

    if (iso_method == "manual") {
      ref_iso <- getVectors(raw_df(ms1),
        rt_sim = rt_iso_threshold,
        mz_sim = mz_iso_threshold
      )
      query_iso <- getVectors(raw_df(ms2),
        rt_sim = rt_iso_threshold,
        mz_sim = mz_iso_threshold
      )
      isolated(ms1) <- raw_df(ms1) |>
        dplyr::filter(.data$Compound_ID %in% ref_iso)
      isolated(ms2) <- raw_df(ms2) |>
        dplyr::filter(.data$Compound_ID %in% query_iso)
    } else if (iso_method == "dbscan") {
      isolated(ms1) <- iso_dbscan(raw_df(ms1), eps)
      isolated(ms2) <- iso_dbscan(raw_df(ms2), eps)
    }
  }
)
