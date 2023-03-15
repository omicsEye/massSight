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
  function(df1,
           df2,
           rt_lower = -.5,
           rt_upper = .5,
           mz_lower = -15,
           mz_upper = 15,
           rt_smooth = .2,
           mz_smooth = .2,
           minimum_intensity = 1000,
           rt_iso_threshold = .5,
           mz_iso_threshold = 50,
           threshold = "manual",
           match_method = "unsupervised",
           smooth_method = "lowess",
           multipliers = c(6, 6, 6),
           weights = c(1, 1, 1),
           keep_features = c(F, F)) {
    df1 <- df1 |>
      dplyr::mutate(MZ = round(MZ, 4),
                    RT = round(RT, 2))

    df2 <- df2 |>
      dplyr::mutate(MZ = round(MZ, 4),
                    RT = round(RT, 2))

    results_list <-
      find_isolated_compounds(
        df1,
        df2,
        rt_lower = rt_lower,
        rt_upper = rt_upper,
        mz_lower = mz_lower,
        mz_upper = mz_upper,
        rt_smooth = rt_smooth,
        mz_smooth = mz_smooth,
        minimum_intensity = minimum_intensity,
        rt_iso_threshold = rt_iso_threshold,
        mz_iso_threshold = mz_iso_threshold,
        threshold = threshold,
        match_method = match_method,
        smooth_method = smooth_method
      )

    results <- results_list[[1]]
    scaled_values <- results_list[[2]]
    cutoffs <- results_list[[3]]

    final_results_list <-
      final_results(
        df1,
        df2,
        scaled_values,
        cutoffs,
        keep_features = keep_features,
        multipliers = multipliers,
        weights = weights
      )

    results_df_complete <- final_results_list[[1]]
    adjusted_df <- final_results_list[[2]]
    message(paste0(
      "Numbers of matched/kept features: ",
      nrow(results_df_complete)
    ))
    return(
      list(
        "results_df_complete" = results_df_complete,
        "adjusted_df" = adjusted_df,
        "smooth_for_plot" = results
      )
    )
  }
