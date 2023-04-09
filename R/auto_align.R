#' @export
#' @title Auto Align
#' @description This function will automatically align your data based on the
#' smoothing method you choose.
#' @param df1 A data frame representing the results of a preprocessed
#' LC-MS experiment.
#' @param df2 A data frame of representing the results of a second
#' preprocessed LC-MS experiment.
#' @param rt_lower A numeric indicating the lower bound of the RT
#' range to be considered for alignment.
#' @param rt_upper A numeric indicating the upper bound of the RT
#' range to be considered for alignment.
#' @param mz_lower A numeric indicating the lower bound of the m/z
#' range to be considered for alignment.
#' @param mz_upper A numeric indicating the upper bound of the m/z
#' range to be considered for alignment.
#' @param rt_smooth A numeric indicating the smoothing parameter for
#' RT.
#' @param mz_smooth A numeric indicating the smoothing parameter for
#' m/z.
#' @param minimum_intensity A numeric indicating the minimum intensity
#' to be considered for alignment.
#' @param rt_iso_threshold A numeric indicating the simplification
#' parameter for RT.
#' @param mz_iso_threshold A numeric indicating the simplification
#' parameter for m/z.
#' @param match_method A character indicating the initial matching method to
#' be used to detect inter-batch variability. Options are "unsupervised" and
#' "supervised".
#' @param smooth_method A character indicating the smoothing method to
#' be used. Options are "lowess", "spline", and "gaussian".
#' @param multipliers A numeric vector indicating the multipliers to be
#' used for the alignment.
#' @param weights A numeric vector indicating the weights to be used for
#' the alignment.
#' @param keep_features A logical vector indicating whether or not to
#' keep features that are not matched.
#' @return A data frame containing the aligned data.
auto_align <-
  function(ref,
           query,
           rt_lower = -.5,
           rt_upper = .5,
           mz_lower = -15,
           mz_upper = 15,
           rt_smooth = .2,
           mz_smooth = .2,
           minimum_intensity = 1000,
           rt_iso_threshold = .5,
           mz_iso_threshold = 5,
           threshold = "manual",
           match_method = "unsupervised",
           smooth_method = "lowess",
           multipliers = c(6, 6, 6),
           weights = c(1, 1, 1),
           keep_features = c(F, F)) {
    raw_df(ref) <- raw_df(ref) |>
      dplyr::mutate(MZ = round(MZ, 4),
                    RT = round(RT, 2))

    raw_df(query) <- raw_df(query) |>
      dplyr::mutate(MZ = round(MZ, 4),
                    RT = round(RT, 2))

    if (match_method == "unsupervised") {
      ref_iso <- get_vectors(raw_df(ref),
                             rt_sim = rt_iso_threshold,
                             mz_sim = mz_iso_threshold)
      query_iso <- get_vectors(raw_df(query),
                               rt_sim = rt_iso_threshold,
                               mz_sim = mz_iso_threshold)
      isolated(ref) <- raw_df(ref) |>
        dplyr::filter(Compound_ID %in% ref_iso)
      isolated(query) <- raw_df(query) |>
        dplyr::filter(Compound_ID %in% query_iso)
    } else if (match_method == "supervised") {
      isolated(ref) <- raw_df(ref) |>
        filter(Metabolite != "")
      isolated(query) <- raw_df(query) |>
        filter(Metabolite != "")
    } else {
      stop("`match_method` must be either 'unsupervised' or 'supervised'.")
    }

    align_obj <- new("MergedMSObject")
    ms1(align_obj) <- ref
    ms2(align_obj) <- query
    align_obj <- align_obj |>
      align_isolated_compounds(match_method = match_method) |>
      smooth_drift(smooth_method = smooth_method,
                   minimum_int = minimum_intensity) |>
      final_results(keep_features = keep_features,
                    multipliers = multipliers,
                    weights = weights)

    message(paste0("Numbers of matched/kept features: ",
                   nrow(all_matched(align_obj))))

    return(align_obj)
  }
