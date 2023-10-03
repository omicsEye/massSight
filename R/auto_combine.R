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
#' @param threshold Something
#' @return A `MergedMSObject` containing the combined data.
auto_combine <-
  function(ms1,
           ms2,
           rt_lower = -.1,
           rt_upper = .1,
           mz_lower = -10,
           mz_upper = 10,
           minimum_intensity = 10,
           iso_method = "manual",
           eps = .1,
           rt_iso_threshold = .5,
           mz_iso_threshold = 5,
           threshold = "manual",
           match_method = "unsupervised",
           smooth_method = "loess",
           weights = c(1, 1, 1),
           keep_features = c(F, F)) {
    # stop conditions ---------------------------------------------------------
    stopifnot("`smooth_method` must be either 'loess' or 'gam'" =
                smooth_method %in% c("loess", "gam"))
    # raw_df(ms1) <- raw_df(ms1) |>
    #   dplyr::mutate(
    #     MZ = round(.data$MZ, 4),
    #     RT = round(.data$RT, 3)
    #   )
    #
    # raw_df(ms2) <- raw_df(ms2) |>
    #   dplyr::mutate(
    #     MZ = round(MZ, 4),
    #     RT = round(RT, 3)
    #   )

    logr::log_open(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M"),
                          ".log"))

    options("logr.notes" = FALSE)
    options("logr.traceback" = FALSE)

    string1 <- glue::glue(
      "auto_combine(ms1,
                             ms2,
                             rt_lower = {rt_lower},
                             rt_upper = {rt_upper},
                             mz_lower = {mz_lower},
                             mz_upper = {mz_upper},
                             minimum_intensity = {minimum_intensity},
                             iso_method = {iso_method},
                             eps = {eps},
                             rt_iso_threshold = {rt_iso_threshold},
                             mz_iso_threshold = {mz_iso_threshold},
                             threshold = {threshold},
                             match_method = {match_method},
                             smooth_method = {smooth_method})"
    )

    logr::log_print("--- massSight Run and Parameters ---",
                    console = F,
                    hide_notes = T)
    logr::log_print(string1, console = F, hide_notes = T)

    if (match_method == "unsupervised") {
      if (iso_method == "manual") {
        ref_iso <- get_vectors(raw_df(ms1),
                               rt_sim = rt_iso_threshold,
                               mz_sim = mz_iso_threshold)
        query_iso <- get_vectors(raw_df(ms2),
                                 rt_sim = rt_iso_threshold,
                                 mz_sim = mz_iso_threshold)
        isolated(ms1) <- raw_df(ms1) |>
          dplyr::filter(.data$Compound_ID %in% ref_iso)
        isolated(ms2) <- raw_df(ms2) |>
          dplyr::filter(.data$Compound_ID %in% query_iso)
      } else if (iso_method == "dbscan") {
        isolated(ms1) <- iso_dbscan(raw_df(ms1), eps)
        isolated(ms2) <- iso_dbscan(raw_df(ms2), eps)
      } else {
        stop("`iso_method` must be either 'manual' or 'dbscan'")
      }
    } else if (match_method == "supervised") {
      isolated(ms1) <- raw_df(ms1) |>
        dplyr::filter(.data$Metabolite != "")
      isolated(ms2) <- raw_df(ms2) |>
        dplyr::filter(.data$Metabolite != "")
    } else {
      stop("`match_method` must be either 'unsupervised' or 'supervised'.")
    }

    align_obj <- methods::new("MergedMSObject")
    ms1(align_obj) <- ms1
    ms2(align_obj) <- ms2
    smooth(align_obj) <- smooth_method
    align_obj <- align_obj |>
      align_pre_isolated_compounds(
        rt_minus = rt_lower,
        rt_plus = rt_upper,
        mz_minus = mz_lower,
        mz_plus = mz_upper
      ) |>
      align_isolated_compounds(
        match_method = match_method,
        rt_minus = rt_lower,
        rt_plus = rt_upper,
        mz_minus = mz_lower,
        mz_plus = mz_upper
      ) |>
      smooth_drift(smooth_method = smooth_method,
                   minimum_int = minimum_intensity) |>
      final_results(keep_features = keep_features,
                    weights = weights)

    logr::log_print(paste0("Numbers of matched/kept features: ",
                   nrow(all_matched(align_obj))), console = T)

    logr::log_close()
    return(align_obj)
  }
