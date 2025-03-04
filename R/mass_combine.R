#' @importFrom ParamHelpers makeParamSet makeNumericParam
#' @importFrom mlrMBO makeMBOControl setMBOControlTermination mbo
#' @export
#' @title Mass Combine
#' @description Combines two `massSight` objects by aligning their features and
#' correcting for systematic differences in retention time and mass-to-charge ratios.
#'
#' @param ms1 A `massSight` object representing the results of a preprocessed
#'   LC-MS experiment that will serve as the reference dataset.
#' @param ms2 A `massSight` object representing the results of a second
#'   preprocessed LC-MS experiment that will be aligned to ms1.
#' @param optimize Logical: whether to optimize alignment parameters using known metabolites.
#'   Default is TRUE
#' @param rt_delta Numeric: the retention time window (+/-) in minutes to consider
#'   when aligning features. Default is 0.5.
#' @param mz_delta Numeric: the mass-to-charge ratio window (+/-) in ppm to consider
#'   when aligning features. Default is 15.
#' @param minimum_intensity Numeric: the minimum intensity threshold for features
#'   to be considered in the alignment. Default is 10.
#' @param iso_method Character: the method for isolating high-quality features before
#'   modeling drift. Options are "manual" or "dbscan". Default is "manual".
#' @param eps Numeric: epsilon parameter for DBSCAN clustering when iso_method = "dbscan".
#'   Default is 0.1.
#' @param rt_iso_threshold Numeric: the retention time similarity threshold for manual
#'   isolation. Default is 0.01.
#' @param mz_iso_threshold Numeric: the m/z similarity threshold for manual isolation.
#'   Default is 2.
#' @param match_method Character: the method for initial feature matching. Options are
#'   "unsupervised" (uses all features) or "supervised" (uses only known metabolites).
#'   Default is "unsupervised".
#' @param smooth_method Character: the method for smoothing systematic drift. Options are
#'   "gam", "bayesian_gam", or "gp". Default is "gam".
#' @param weights Numeric vector: weights for RT, m/z, and intensity in the alignment
#'   scoring. Default is c(1, 1, 1).
#' @param log Character: path to save the log file. Set to NULL to disable logging.
#'   Default is NULL.
#' @param output Character: directory to save output files. Default is NULL (current directory).
#' @param n_iter Integer: number of optimization iterations if optimize=TRUE. Default is 50.
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A `MergedMSObject` containing:
#'   \itemize{
#'     \item Original ms1 and ms2 objects
#'     \item Matched features between datasets
#'     \item Drift correction models
#'     \item Alignment quality metrics
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' combined <- mass_combine(ms1_data, ms2_data)
#'
#' # Use optimization with known metabolites
#' combined <- mass_combine(ms1_data, ms2_data,
#'                         optimize = TRUE,
#'                         n_iter = 50)
#'
#' # Custom parameters for stricter matching
#' combined <- mass_combine(ms1_data, ms2_data,
#'                         rt_delta = 0.3,
#'                         mz_delta = 10,
#'                         minimum_intensity = 100)
#' }
#'
#' @export
mass_combine <- function(ms1,
                         ms2,
                         optimize = TRUE,
                         rt_delta = 0.5,
                         mz_delta = 15,
                         minimum_intensity = 10,
                         iso_method = "manual",
                         eps = .1,
                         rt_iso_threshold = .01,
                         mz_iso_threshold = 2,
                         match_method = "unsupervised",
                         smooth_method = "gam",
                         weights = c(1, 1, 1),
                         log = NULL,
                         output = NULL,
                         n_iter = 50,
                         ...) {
  if (!is.null(log)) {
    time_start <- Sys.time()
    check_jsonlite()
    log_params <- match.call(expand.dots = TRUE) %>%
      modify_call() %>%
      as.list() %>%
      (\(x) x[-1])() %>%
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

  # check if ms1 and ms2 are massSight objects
  if (!is(ms1, "MSObject") || !is(ms2, "MSObject")) {
    stop("ms1 and ms2 must be massSight objects")
  }

  validate_parameters(iso_method, match_method, smooth_method, minimum_intensity)

  # Either optimize parameters or use provided/default parameters
  if (optimize) {
    message("Optimizing parameters using Bayesian optimization...")
    opt_result <- optimize_parameters(ms1, ms2, ...)
    params <- opt_result$parameters
    align_obj <- opt_result$best_object # Use the stored best object

    message(sprintf(
      "Optimization complete. Final score: %.3f",
      opt_result$final_score
    ))
    message("\nOptimal parameters:")
    message(sprintf("  RT delta: %.3f", params[["rt_delta"]]))
    message(sprintf("  MZ delta: %.3f", params[["mz_delta"]]))
    message(sprintf("  RT isolation threshold: %.3f", params[["rt_iso_threshold"]]))
    message(sprintf("  MZ isolation threshold: %.3f", params[["mz_iso_threshold"]]))
    message(sprintf("  Alpha rank: %.3f", params[["alpha_rank"]]))
    message(sprintf("  Alpha RT: %.3f", params[["alpha_rt"]]))
    message(sprintf("  Alpha MZ: %.3f", params[["alpha_mz"]]))
  } else {
    params <- list(
      rt_delta = rt_delta,
      mz_delta = mz_delta,
      rt_iso_threshold = rt_iso_threshold,
      mz_iso_threshold = mz_iso_threshold,
      alpha_rank = .5,
      # default values when not optimizing
      alpha_rt = .5,
      alpha_mz = .5
    )

    if (match_method == "unsupervised") {
      if (iso_method == "manual") {
        ref_iso <- getVectors(raw_df(ms1), rt_sim = params[["rt_iso_threshold"]], # Use optimized parameter
                              mz_sim = params[["mz_iso_threshold"]]) # Use optimized parameter
        query_iso <- getVectors(raw_df(ms2), rt_sim = params[["rt_iso_threshold"]], # Use optimized parameter
                                mz_sim = params[["mz_iso_threshold"]]) # Use optimized parameter
        isolated(ms1) <- raw_df(ms1) %>%
          dplyr::filter(.data$Compound_ID %in% ref_iso)
        isolated(ms2) <- raw_df(ms2) %>%
          dplyr::filter(.data$Compound_ID %in% query_iso)
      } else if (iso_method == "dbscan") {
        isolated(ms1) <- iso_dbscan(raw_df(ms1), eps)
        isolated(ms2) <- iso_dbscan(raw_df(ms2), eps)
      }
    } else if (match_method == "supervised") {
      isolated(ms1) <- raw_df(ms1) %>%
        dplyr::filter(.data$Metabolite != "")
      isolated(ms2) <- raw_df(ms2) %>%
        dplyr::filter(.data$Metabolite != "")
    }

    align_obj <- methods::new("MergedMSObject")
    ms1(align_obj) <- ms1
    ms2(align_obj) <- ms2
    align_obj <- align_obj %>%
      align_isolated_compounds(
        match_method = match_method,
        rt_delta = params[["rt_delta"]],
        mz_delta = params[["mz_delta"]]
      ) %>%
      smooth_drift(smooth_method = smooth_method, minimum_int = minimum_intensity) %>%
      final_results(
        rt_threshold = params[["rt_delta"]],
        mz_threshold = params[["mz_delta"]],
        alpha_rank = params[["alpha_rank"]],
        alpha_rt = params[["alpha_rt"]],
        alpha_mz = params[["alpha_mz"]]
      )
  }




  if (!is.null(log)) {
    log_parameters(log,
                   log_params,
                   log_r,
                   log_date,
                   align_obj,
                   ms1,
                   ms2,
                   time_start)
  }
  if (optimize) {
    # Store optimization history in the object attributes
    attr(align_obj, "optimization") <- list(
      parameters = opt_result$parameters,
      final_score = opt_result$final_score,
      history = opt_result$optimization_history
    )
  }
  return(align_obj)
}

mad_based_outlier <- function(points, thresh = 5) {
  diff <- sqrt((points - stats::median(points, na.rm = TRUE)) ^ 2)
  med_abs_deviation <- stats::median(diff, na.rm = TRUE)
  if (med_abs_deviation == 0) {
    mod_z_score <- rep(0, length(diff))
  } else {
    mod_z_score <- .6745 * diff / med_abs_deviation
  }
  return(mod_z_score > thresh)
}

scale_intensity <- function(data, intensity) {
  return(data / intensity)
}

scale_intensity_parameters <-
  function(data_int1, data_int2, min_int = 0) {
    min_int <- log10(min_int)
    # Filter both vectors using the minimum intensity threshold
    valid_idx <- data_int1 > min_int & data_int2 > min_int &
      !is.na(data_int1) & !is.na(data_int2)

    # Use only valid indices for both vectors
    data_int1_filtered <- data_int1[valid_idx]
    data_int2_filtered <- data_int2[valid_idx]

    # Calculate mean ratio only for valid pairs
    fit <- mean(data_int2_filtered / data_int1_filtered, na.rm = TRUE)
    return(fit)
  }

scale_smooth <- function(query_values, smooth_x, smooth_y) {
  suppressWarnings(f <- stats::approx(
    x = smooth_x,
    y = smooth_y,
    xout = query_values,
    rule = 2
  ))

  query_out <- query_values - f$y
  return(query_out)
}

align_isolated_compounds <-
  function(align_ms_obj,
           match_method,
           rt_delta = 0.5,
           mz_delta = 15,
           keep = FALSE) {
    df1 <- isolated(ms1(align_ms_obj))
    df2 <- isolated(ms2(align_ms_obj))
    if (match_method == "unsupervised") {
      df1 <- df1 %>%
        dplyr::rename(dplyr::any_of(
          c(
            "RT_1" = "RT",
            "MZ_1" = "MZ",
            "Intensity_1" = "Intensity",
            "Compound_ID_1" = "Compound_ID",
            "Metabolite_1" = "Metabolite"
          )
        ))
      df2 <- df2 %>% dplyr::rename(dplyr::any_of(
        c(
          "RT_2" = "RT",
          "MZ_2" = "MZ",
          "Intensity_2" = "Intensity",
          "Compound_ID_2" = "Compound_ID",
          "Metabolite_2" = "Metabolite"
        )
      ))
      results <- df1 %>%
        dplyr::mutate(
          rt_upper = RT_1 + rt_delta,
          rt_lower = RT_1 - rt_delta,
          mz_upper = MZ_1 + (MZ_1 * mz_delta / 1e6),
          mz_lower = MZ_1 - (MZ_1 * mz_delta / 1e6)
        ) %>%
        dplyr::inner_join(df2,
                          by = dplyr::join_by(
                            rt_lower < RT_2,
                            rt_upper > RT_2,
                            mz_lower < MZ_2,
                            mz_upper > MZ_2
                          )) %>%
        dplyr::select(-rt_upper, -rt_lower, -mz_upper, -mz_lower)
    } else if (match_method == "supervised") {
      stopifnot("Metabolite" %in% colnames(df1) &
                  "Metabolite" %in% colnames(df2))
      vec_1 <- df1 %>%
        dplyr::rename(df1 = .data$Compound_ID) %>%
        dplyr::filter(.data$Metabolite != "")
      vec_2 <- df2 %>%
        dplyr::rename(
          RT_2 = .data$RT_1,
          MZ_2 = .data$MZ_1,
          Intensity_2 = .data$Intensity,
          df2 = .data$Compound_ID
        ) %>%
        dplyr::filter(.data$Metabolite != "")
      results <- vec_1 %>%
        dplyr::inner_join(vec_2, by = c("Metabolite"))
    }
    iso_matched(align_ms_obj) <- results
    return(align_ms_obj)
  }

final_results <- function(align_ms_obj,
                          rt_threshold,
                          mz_threshold,
                          alpha_rank,
                          alpha_rt,
                          alpha_mz) {
  study1_name <- name(ms1(align_ms_obj))
  study2_name <- name(ms2(align_ms_obj))
  df1 <- raw_df(ms1(align_ms_obj))
  df2 <- raw_df(ms2(align_ms_obj))
  scaled_df <- scaled_values(align_ms_obj)
  df2$RT_adj_2 <- scaled_df$RT
  df2$MZ_adj_2 <- scaled_df$MZ
  scaled_df$Intensity_2 <- df2$Intensity
  scaled_df$Metabolite_2 <- df2$Metabolite
  stds <- cutoffs(align_ms_obj)

  # Create a joined matrix of potential matches
  message("Creating potential final matches")
  potential_matches <- find_all_matches(
    df1,
    df2,
    rt_threshold = rt_threshold,
    mz_threshold = mz_threshold,
    alpha_rank = alpha_rank,
    alpha_rt = alpha_rt,
    alpha_mz = alpha_mz
  )
  # Calculate match scores for all potential matches
  message("Calculating match scores")
  potential_matches <- potential_matches %>%
    dplyr::mutate(delta_RT = RT_adj_2 - RT_1,
                  delta_MZ = (MZ_adj_2 - MZ_1) / MZ_1 * 1e6,
                  # Convert to ppm)

                  best_matches <- potential_matches
                  # Prepare the final results dataframe
                  results <- best_matches %>%
                    dplyr::rename(dplyr::any_of(
                      c(
                        # Rename columns for study1
                        rlang::set_names("Compound_ID_1", paste0("Compound_ID_", study1_name)),
                        rlang::set_names("RT_1", paste0("RT_", study1_name)),
                        rlang::set_names("MZ_1", paste0("MZ_", study1_name)),
                        rlang::set_names("Intensity_1", paste0("Intensity_", study1_name)),
                        rlang::set_names("Metabolite_1", paste0("Metabolite_", study1_name)),

                        # Rename columns for study2
                        rlang::set_names("Compound_ID_2", paste0("Compound_ID_", study2_name)),
                        rlang::set_names("RT_2", paste0("RT_", study2_name)),
                        rlang::set_names("RT_adj_2", paste0("RT_adj_", study2_name)),
                        rlang::set_names("MZ_2", paste0("MZ_", study2_name)),
                        rlang::set_names("MZ_adj_2", paste0("MZ_adj_", study2_name)),
                        rlang::set_names("Intensity_2", paste0("Intensity_", study2_name)),
                        rlang::set_names("Metabolite_2", paste0("Metabolite_", study2_name))
                      )
                    ))

                  # Add representative columns

                  unmatched_study1 <- df1 %>%
                    dplyr::anti_join(results, by = setNames(paste0("Compound_ID_", study1_name), "Compound_ID")) %>%
                    dplyr::rename_with( ~ paste0(.x, "_", study1_name), .cols = dplyr::everything()) %>%
                    dplyr::mutate(
                      !!rlang::sym(paste0("Compound_ID_", study2_name)) := NA,!!rlang::sym(paste0("RT_", study2_name)) := NA,!!rlang::sym(paste0("MZ_", study2_name)) := NA,!!rlang::sym(paste0("Intensity_", study2_name)) := NA,!!rlang::sym(paste0("Metabolite_", study2_name)) := NA,
                      rep_Compound_ID = !!rlang::sym(paste0("Compound_ID_", study1_name)),
                      rep_RT = !!rlang::sym(paste0("RT_", study1_name)),
                      rep_MZ = !!rlang::sym(paste0("MZ_", study1_name)),
                      rep_Intensity = if ("Intensity" %in% names(df1))
                        !!rlang::sym(paste0("Intensity_", study1_name))
                      else
                        NULL,
                      rep_Metabolite = if ("Metabolite" %in% names(df1))
                        !!rlang::sym(paste0("Metabolite_", study1_name))
                      else
                        NULL
                    )

                  # Get unmatched metabolites from study2
                  unmatched_study2 <- df2 %>%
                    dplyr::select(-c(RT_adj_2, MZ_adj_2)) %>%
                    dplyr::anti_join(results, by = setNames(paste0("Compound_ID_", study2_name), "Compound_ID")) %>%
                    dplyr::rename_with( ~ paste0(.x, "_", study2_name), .cols = dplyr::everything()) %>%
                    dplyr::mutate(
                      !!rlang::sym(paste0("Compound_ID_", study1_name)) := NA,!!rlang::sym(paste0("RT_", study1_name)) := NA,!!rlang::sym(paste0("MZ_", study1_name)) := NA,!!rlang::sym(paste0("Intensity_", study1_name)) := NA,!!rlang::sym(paste0("Metabolite_", study1_name)) := NA,
                      rep_Compound_ID = !!rlang::sym(paste0("Compound_ID_", study2_name)),
                      rep_RT = !!rlang::sym(paste0("RT_", study2_name)),
                      rep_MZ = !!rlang::sym(paste0("MZ_", study2_name)),
                      rep_Intensity = if ("Intensity" %in% names(df2))
                        !!rlang::sym(paste0("Intensity_", study2_name))
                      else
                        NULL,
                      rep_Metabolite = if ("Metabolite" %in% names(df2))
                        !!rlang::sym(paste0("Metabolite_", study2_name))
                      else
                        NULL
                    )

                  # Combine matched and unmatched results
                  all_results <- list(results, unmatched_study1, unmatched_study2) %>%
                    purrr::keep( ~ nrow(.) > 0) %>%
                    dplyr::bind_rows()
                  all_results <- all_results %>%
                    dplyr::mutate(
                      rep_Compound_ID = ifelse(
                        is.na(!!rlang::sym(paste0(
                          "Compound_ID_", study1_name
                        ))),!!rlang::sym(paste0("Compound_ID_", study2_name)),!!rlang::sym(paste0("Compound_ID_", study1_name))
                      ),
                      rep_RT = ifelse(
                        is.na(!!rlang::sym(paste0(
                          "RT_", study1_name
                        ))),!!rlang::sym(paste0("RT_", study2_name)),!!rlang::sym(paste0("RT_", study1_name))
                      ),
                      rep_MZ = ifelse(
                        is.na(!!rlang::sym(paste0(
                          "MZ_", study1_name
                        ))),!!rlang::sym(paste0("MZ_", study2_name)),!!rlang::sym(paste0("MZ_", study1_name))
                      ),
                      rep_Intensity = if ("Intensity" %in% names(df1)) {
                        ifelse(is.na(!!rlang::sym(paste0(
                          "Intensity_", study1_name
                        ))),!!rlang::sym(paste0("Intensity_", study2_name)),!!rlang::sym(paste0("Intensity_", study1_name)))
                      } else {
                        NULL
                      },
                      rep_Metabolite = if ("Metabolite" %in% names(df1)) {
                        ifelse(is.na(!!rlang::sym(paste0(
                          "Metabolite_", study1_name
                        ))),!!rlang::sym(paste0("Metabolite_", study2_name)),!!rlang::sym(paste0("Metabolite_", study1_name)))
                      } else {
                        NULL
                      }
                    )
                  # Reorder columns
                  all_results <- all_results %>%
                    dplyr::mutate(matched = ifelse(is.na(!!rlang::sym(
                      paste0("Compound_ID_", study1_name)
                    )) |
                      is.na(!!rlang::sym(
                        paste0("Compound_ID_", study2_name)
                      )), FALSE, TRUE)) %>%
                    dplyr::select(dplyr::starts_with("rep_"), matched, dplyr::everything())

                  # match metadata from ms1 and ms2 to final results
                  metadata2 <- metadata(ms2(align_ms_obj))
                  metadata1 <- metadata(ms1(align_ms_obj))
                  if (nrow(metadata1) > 0) {
                    all_results <- all_results %>%
                      dplyr::left_join(metadata1, by = structure("Compound_ID", names = paste0("Compound_ID_", study1_name)))
                  }

                  if (nrow(metadata2) > 0) {
                    all_results <- all_results %>%
                      dplyr::left_join(metadata2, by = structure("Compound_ID", names = paste0("Compound_ID_", study2_name)))
                  }

                  all_matched(align_ms_obj) <- all_results
                  adjusted_df(align_ms_obj) <- scaled_df

                  return(align_ms_obj)
}

find_all_matches <- function(ref,
                             query,
                             rt_threshold,
                             mz_threshold,
                             alpha_rank,
                             alpha_rt,
                             alpha_mz) {
  denom <- exp(alpha_rank) + exp(alpha_rt) + exp(alpha_mz)
  rank_weight <- exp(alpha_rank) / denom
  rt_weight <- exp(alpha_rt) / denom
  mz_weight <- exp(alpha_mz) / denom

  ref <- ref %>%
    dplyr::rename(dplyr::any_of(
      c(
        "RT_1" = "RT",
        "MZ_1" = "MZ",
        "Compound_ID_1" = "Compound_ID",
        "Intensity_1" = "Intensity",
        "Metabolite_1" = "Metabolite"
      )
    ))
  query <- query %>%
    dplyr::rename(dplyr::any_of(
      c(
        "RT_2" = "RT",
        "MZ_2" = "MZ",
        "Compound_ID_2" = "Compound_ID",
        "Intensity_2" = "Intensity",
        "Metabolite_2" = "Metabolite"
      )
    ))

  # Add normalized RT rank information (0 to 1 scale)
  ref <- ref %>%
    dplyr::arrange(RT_1) %>%
    dplyr::mutate(RT_rank_1 = (dplyr::row_number() - 1) / (dplyr::n() - 1))

  query <- query %>%
    dplyr::arrange(RT_adj_2) %>%
    dplyr::mutate(RT_rank_2 = (dplyr::row_number() - 1) / (dplyr::n() - 1))

  matches <- ref %>%
    dplyr::mutate(
      mz_upper = MZ_1 + mz_threshold * MZ_1 / 1e6,
      mz_lower = MZ_1 - mz_threshold * MZ_1 / 1e6,
      rt_upper = RT_1 + rt_threshold,
      rt_lower = RT_1 - rt_threshold
    ) %>%
    dplyr::left_join(
      query,
      dplyr::join_by(
        mz_lower < MZ_adj_2,
        mz_upper > MZ_adj_2,
        rt_lower < RT_adj_2,
        rt_upper > RT_adj_2
      )
    ) %>%
    dplyr::select(-mz_upper, -mz_lower, -rt_upper, -rt_lower) %>%
    # Calculate normalized rank difference
    dplyr::mutate(
      # Calculate differences
      rank_diff = abs(RT_rank_2 - RT_rank_1),
      rt_diff = abs(RT_adj_2 - RT_1),
      mz_diff = abs((MZ_adj_2 - MZ_1) / MZ_1 * 1e6),

      # Min-max scaling with NA removal and inversion (1 - score)
      # This way, smaller differences = scores closer to 1
      rank_score = (rank_diff - min(rank_diff, na.rm = TRUE)) /
        (max(rank_diff, na.rm = TRUE) - min(rank_diff, na.rm = TRUE)),
      rt_score = (rt_diff - min(rt_diff, na.rm = TRUE)) /
        (max(rt_diff, na.rm = TRUE) - min(rt_diff, na.rm = TRUE)),
      mz_score = (mz_diff - min(mz_diff, na.rm = TRUE)) /
        (max(mz_diff, na.rm = TRUE) - min(mz_diff, na.rm = TRUE)),
      score = rank_weight * rank_score + rt_weight * rt_score + mz_weight * mz_score
    )
  return(matches)
}

get_cutoffs <- function(df1, df2, has_int = TRUE) {
  data_rt <- df2$RT - df1$RT_1
  data_mz <- (df2$MZ - df1$MZ_1) / df1$MZ_1 * 1e6
  if (has_int) {
    data_int <- log10(df2$Intensity) - log10(df1$Intensity_1)
    not_outliers <- !mad_based_outlier(data_rt) &
      !mad_based_outlier(data_mz) &
      !mad_based_outlier(data_int)
    data_int <- replace(data_int, is.infinite(data_int), NA)
    cutoffs <- c(sd(data_rt[not_outliers], na.rm = TRUE),
                 sd(data_mz[not_outliers], na.rm = TRUE),
                 sd(data_int[not_outliers], na.rm = TRUE))
  } else {
    not_outliers <- !mad_based_outlier(data_rt) &
      !mad_based_outlier(data_mz)
    cutoffs <- c(sd(data_rt[not_outliers], na.rm = TRUE), sd(data_mz[not_outliers], na.rm = TRUE))
  }
  # Properly identify outliers
  rt_outliers <- df1$Compound_ID_1[mad_based_outlier(data_rt)]
  mz_outliers <- df1$Compound_ID_1[mad_based_outlier(data_mz)]
  if (has_int) {
    int_outliers <- df1$Compound_ID_1[mad_based_outlier(data_int)]
    outliers <- unique(c(rt_outliers, mz_outliers, int_outliers))
  } else {
    outliers <- unique(c(rt_outliers, mz_outliers))
  }

  return(list("cutoffs" = cutoffs, "outliers" = outliers))
}

smooth_drift <- function(align_ms_obj,
                         smooth_method,
                         minimum_int,
                         mz_bin_size = 100) {
  df1 <- align_ms_obj %>%
    ms1() %>%
    raw_df()
  df2 <- align_ms_obj %>%
    ms2() %>%
    raw_df()
  results <- iso_matched(align_ms_obj)

  # RT drift correction remains the same
  results <- results %>%
    dplyr::arrange(.data$RT_1) %>%
    dplyr::mutate(delta_RT = .data$RT_2 - .data$RT_1)

  # RT smoothing methods remain unchanged
  if (smooth_method == "bayesian_gam") {
    message("Starting bayesian GAM smoothing for RT")
    gam_data <- data.frame(RT_1 = results$RT_1, delta_RT = results$delta_RT)

    # Define the model formula
    formula <- brms::bf(delta_RT ~ s(RT_1))

    # Fit the Bayesian GAM
    brms_model <- brms::brm(
      formula = formula,
      data = gam_data,
      family = stats::gaussian(),
      cores = 4,
      chains = 4,
      iter = 2000,
      warmup = 1000,
      # refresh = 0,
      # silent = 2,
      # open_progress = FALSE,
      control = list(adapt_delta = 0.95),
    )

    # Generate predictions
    smooth_x_rt <- results$RT_1
    smooth_y_rt <- brms::posterior_predict(brms_model, newdata = data.frame(RT_1 = smooth_x_rt))
    smooth_y_rt <- colMeans(smooth_y_rt) # Use posterior mean as point estimate

    # Store the smoothing results
    smooth_method(align_ms_obj)[["rt_x"]] <- smooth_x_rt
    smooth_method(align_ms_obj)[["rt_y"]] <- smooth_y_rt
  } else if (smooth_method == "gam") {
    message("GAM smoothing for RT drift")
    gam_fit <- mgcv::gam(
      delta_RT ~ s(RT_1),
      data = results,
      family = mgcv::scat(),
      method = "REML"
    )

    smooth_x_rt <- results$RT_1
    smooth_y_rt <- predict(gam_fit, newdata = data.frame(RT_1 = smooth_x_rt))
    smooth_method(align_ms_obj)[["rt_x"]] <- smooth_x_rt
    smooth_method(align_ms_obj)[["rt_y"]] <- smooth_y_rt
  } else if (smooth_method == "lm") {
    message("Linear smoothing for RT drift")
    lm_fit <- MASS::rlm(delta_RT ~ RT_1, data = results)
    smooth_x_rt <- results$RT_1
    smooth_y_rt <- predict(lm_fit, newdata = data.frame(RT_1 = smooth_x_rt))
    smooth_method(align_ms_obj)[["rt_x"]] <- smooth_x_rt
    smooth_method(align_ms_obj)[["rt_y"]] <- smooth_y_rt
  } else if (smooth_method == "gp") {
    message("Starting gaussian smoothing for RT")

    gp_fit <- brms::brm(delta_RT ~ gp(RT, cov = "matern52", scale = TRUE),
                        results,
                        algorithm = "meanfield")

    smooth_x_rt <- results$RT_1
    smooth_y_rt <- brms::posterior_predict(gp_fit, newdata = data.frame(RT_1 = smooth_x_rt))
    smooth_y_rt <- colMeans(smooth_y_rt)

    smooth_method(align_ms_obj)[["rt_x"]] <- smooth_x_rt
    smooth_method(align_ms_obj)[["rt_y"]] <- smooth_y_rt
  }

  # Apply RT corrections as before
  scaled_rts <- scale_smooth(df2$RT, smooth_x_rt + smooth_y_rt, smooth_y_rt)
  scaled_rts_res <- scale_smooth(results$RT_1, smooth_x_rt + smooth_y_rt, smooth_y_rt)

  results <- results %>%
    dplyr::mutate(smooth_rt = smooth_y_rt, srt = scaled_rts_res)

  ## MZ correction using smoothing approach
  message("Starting mass error correction")
  results <- results %>%
    dplyr::mutate(mass_error_ppm = (MZ_2 - MZ_1) / MZ_1 * 1e6)

  # MZ smoothing using the selected method
  if (smooth_method == "bayesian_gam") {
    message("Starting bayesian GAM smoothing for mass error")
    gam_data <- data.frame(MZ_1 = results$MZ_1,
                           mass_error_ppm = results$mass_error_ppm)

    # Define the model formula for mass error
    formula <- brms::bf(mass_error_ppm ~ s(MZ_1))

    # Fit the Bayesian GAM for mass error
    brms_model_mz <- brms::brm(
      formula = formula,
      data = gam_data,
      family = stats::gaussian(),
      cores = 4,
      chains = 4,
      iter = 2000,
      warmup = 1000,
      control = list(adapt_delta = 0.95)
    )

    # Generate predictions for mass error
    smooth_x_mz <- results$MZ_1
    smooth_y_mz <- brms::posterior_predict(brms_model_mz, newdata = data.frame(MZ_1 = smooth_x_mz))
    smooth_y_mz <- colMeans(smooth_y_mz) # Use posterior mean as point estimate

  } else if (smooth_method == "gam") {
    message("GAM smoothing for mass error")
    gam_fit_mz <- mgcv::gam(
      mass_error_ppm ~ s(MZ_1),
      data = results,
      family = mgcv::scat(),
      method = "REML"
    )

    smooth_x_mz <- results$MZ_1
    smooth_y_mz <- predict(gam_fit_mz, newdata = data.frame(MZ_1 = smooth_x_mz))

  } else if (smooth_method == "lm") {
    message("Linear smoothing for mass error")
    lm_fit_mz <- MASS::rlm(mass_error_ppm ~ MZ_1, data = results)
    smooth_x_mz <- results$MZ_1
    smooth_y_mz <- predict(lm_fit_mz, newdata = data.frame(MZ_1 = smooth_x_mz))

  } else if (smooth_method == "gp") {
    message("Starting gaussian process smoothing for mass error")
    gp_fit_mz <- brms::brm(mass_error_ppm ~ gp(MZ_1, cov = "matern52", scale = TRUE),
                           results,
                           algorithm = "meanfield")

    smooth_x_mz <- results$MZ_1
    smooth_y_mz <- brms::posterior_predict(gp_fit_mz, newdata = data.frame(MZ_1 = smooth_x_mz))
    smooth_y_mz <- colMeans(smooth_y_mz)
  }

  # Store smoothing results for MZ
  smooth_method(align_ms_obj)[["mz_x"]] <- smooth_x_mz
  smooth_method(align_ms_obj)[["mz_y"]] <- smooth_y_mz

  # Apply MZ corrections using smoothed values
  # For all points in df2, interpolate the correction
  suppressWarnings({
    mz_corrections <- stats::approx(
      x = smooth_x_mz,
      y = smooth_y_mz,
      xout = df2$MZ,
      rule = 2
    )$y
  })

  # Apply corrections to all masses
  scaled_mzs <- df2$MZ / (1 + mz_corrections / 1e6)
  scaled_mzs_res <- results$MZ_1 / (1 + smooth_y_mz / 1e6)

  results <- results %>%
    dplyr::mutate(smooth_mz = smooth_y_mz, smz = scaled_mzs_res)

  # Intensity scaling remains unchanged
  if ("Intensity_1" %in% names(results)) {
    temp_df1_int <- log10(results$Intensity_1)
    temp_df2_int <- log10(results$Intensity_2)

    intensity_parameters <- scale_intensity_parameters(temp_df1_int, temp_df2_int, min_int = minimum_int)

    scaled_vector_intensity <- scale_intensity(temp_df2_int, intensity_parameters)
    scaled_vector_intensity <- 10 ^ scaled_vector_intensity
    results$sintensity <- scaled_vector_intensity

    log_df2 <- log10(df2$Intensity)
    scaled_intensity <- scale_intensity(log_df2, intensity_parameters)
    scaled_intensity <- 10 ^ scaled_intensity

    scaled_df <- data.frame(
      "RT" = results$srt,
      "MZ" = results$smz,
      "Intensity" = results$sintensity
    )
    scaled_values <- data.frame(
      "Compound_ID" = df2$Compound_ID,
      "RT" = scaled_rts,
      "MZ" = scaled_mzs,
      "Intensity" = scaled_intensity
    )
  } else {
    scaled_df <- data.frame("RT" = results$srt, "MZ" = results$smz)
    scaled_values <- data.frame(
      "Compound_ID" = df2$Compound_ID,
      "RT" = scaled_rts,
      "MZ" = scaled_mzs
    )
  }

  # Calculate final deviations
  dev_out <- get_cutoffs(
    df1 = results %>%
      dplyr::select(dplyr::any_of(
        c("Compound_ID_1", "RT_1", "MZ_1", "Intensity_1")
      )),
    df2 = scaled_df,
    has_int = ("Intensity_1" %in% names(results))
  )

  # Update object
  iso_matched(align_ms_obj) <- results
  scaled_values(align_ms_obj) <- scaled_values
  cutoffs(align_ms_obj) <- dev_out$cutoffs

  return(align_ms_obj)
}

optimize_parameters <- function(ms1,
                                ms2,
                                n_iter = 50,
                                match_method = "unsupervised",
                                iso_method = "manual",
                                smooth_method = "gam",
                                minimum_intensity = 10) {
  message("Initializing optimization...")

  # Define parameter space
  param_set <- ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam("rt_delta", lower = 0.1, upper = 1.0),
    ParamHelpers::makeNumericParam("mz_delta", lower = 1, upper = 20),
    ParamHelpers::makeNumericParam("rt_iso_threshold", lower = 0.01, upper = 0.1),
    ParamHelpers::makeNumericParam("mz_iso_threshold", lower = 1, upper = 5),
    ParamHelpers::makeNumericParam("alpha_rank", lower = -2, upper = 2),
    ParamHelpers::makeNumericParam("alpha_rt", lower = -2, upper = 2),
    ParamHelpers::makeNumericParam("alpha_mz", lower = -2, upper = 2)
  )

  # Create progress bar with more informative format
  total_iters <- n_iter + 20  # Initial design + iterations
  pb <- progress::progress_bar$new(
    format = "  Optimization [:bar] :percent | Iter :current/:total | Best score: :best_score | Elapsed: :elapsed",
    total = total_iters,
    clear = FALSE,
    width = 100
  )

  # Create an environment to store optimization state
  opt_state <- new.env(parent = emptyenv())
  opt_state$best_score <- -Inf
  opt_state$best_params <- NULL
  opt_state$best_object <- NULL  # Add storage for best MergedMSObject
  opt_state$target_achieved <- FALSE

  # Define objective function with progress updates
  obj_fun <- smoof::makeSingleObjectiveFunction(
    name = "alignment_score",
    fn = function(x) {
      # Check if we already found a perfect score
      if (opt_state$target_achieved) {
        return(1)
      }

      score <- suppressMessages(try({
        if (match_method == "unsupervised") {
          if (iso_method == "manual") {
            ref_iso <- getVectors(raw_df(ms1),
                                  rt_sim = x[["rt_iso_threshold"]],
                                  mz_sim = x[["mz_iso_threshold"]])
            query_iso <- getVectors(raw_df(ms2),
                                    rt_sim = x[["rt_iso_threshold"]],
                                    mz_sim = x[["mz_iso_threshold"]])
            isolated(ms1) <- raw_df(ms1) %>%
              dplyr::filter(.data$Compound_ID %in% ref_iso)
            isolated(ms2) <- raw_df(ms2) %>%
              dplyr::filter(.data$Compound_ID %in% query_iso)
          }
        } else if (match_method == "supervised") {
          isolated(ms1) <- raw_df(ms1) %>%
            dplyr::filter(.data$Metabolite != "")
          isolated(ms2) <- raw_df(ms2) %>%
            dplyr::filter(.data$Metabolite != "")
        }

        align_obj <- methods::new("MergedMSObject")
        ms1(align_obj) <- ms1
        ms2(align_obj) <- ms2
        align_obj <- align_obj %>%
          align_isolated_compounds(
            match_method = match_method,
            rt_delta = x[["rt_delta"]],
            mz_delta = x[["mz_delta"]]
          ) %>%
          smooth_drift(smooth_method = smooth_method, minimum_int = minimum_intensity) %>%
          final_results(
            rt_threshold = x[["rt_delta"]],
            mz_threshold = x[["mz_delta"]],
            alpha_rank = x[["alpha_rank"]],
            alpha_rt = x[["alpha_rt"]],
            alpha_mz = x[["alpha_mz"]]
          )

        score <- evaluate_matches(align_obj)

        # Update best object if score is better
        if (score > opt_state$best_score) {
          opt_state$best_object <- align_obj
        }

        score
      }, silent = TRUE))

      # Update best score and parameters if current score is better
      if (!inherits(score, "try-error") &&
          !is.infinite(score) && !is.na(score)) {
        if (score > opt_state$best_score && !opt_state$target_achieved) {
          opt_state$best_score <- score
          opt_state$best_params <- x
          # Check if we've found a perfect score
          if (score >= .99 && !opt_state$target_achieved) {
            opt_state$target_achieved <- TRUE
            message("\nTarget score achieved! Stopping optimization.")
            return(1)
          }
        }
      }

      # Tick progress bar with updated best score
      pb$tick(tokens = list(best_score = sprintf("%.4f", opt_state$best_score)))

      if (inherits(score, "try-error") ||
          is.infinite(score) || is.na(score)) {
        return(0)
      }
      return(score)
    },
    par.set = param_set,
    minimize = FALSE,
  )

  # Configure MBO control
  control <- mlrMBO::makeMBOControl() %>%
    mlrMBO::setMBOControlTermination(iters = n_iter, target.fun.value = .99)

  # Run optimization
  suppressWarnings({
    design <- ParamHelpers::generateDesign(n = 20, par.set = param_set) %>%
      as.data.frame() %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x + runif(dplyr::n(), -1e-6, 1e-6)))

    result <- mlrMBO::mbo(
      fun = obj_fun,
      design = design,
      control = control,
      show.info = FALSE
    )
  })

  # Clear progress bar
  pb$terminate()

  return(
    list(
      parameters = opt_state$best_params,
      final_score = opt_state$best_score,
      optimization_history = as.data.frame(result$opt.path),
      best_object = opt_state$best_object  # Return the best MergedMSObject
    )
  )
}


#' Get pairs of known metabolites between datasets
#' @param ms1 First MSObject
#' @param ms2 Second MSObject
#' @return Data frame of matching compound IDs for known metabolites
get_known_pairs <- function(merged_ms) {
  # Get data frames from MSObjects
  df1 <- raw_df(ms1(merged_ms))
  df2 <- raw_df(ms2(merged_ms))

  metabs1 <- df1$Metabolite[!is.na(df1$Metabolite)]
  metabs2 <- df2$Metabolite[!is.na(df2$Metabolite)]

  # Find shared metabolites
  shared_metabolites <- intersect(df1$Metabolite[!is.na(df1$Metabolite)], df2$Metabolite[!is.na(df2$Metabolite)])

  # Create pairs dataframe
  pairs <- data.frame(
    metabolite = shared_metabolites,
    ms1_id = df1$Compound_ID[match(shared_metabolites, df1$Metabolite)],
    ms2_id = df2$Compound_ID[match(shared_metabolites, df2$Metabolite)]
  ) %>%
    dplyr::distinct()


  return(pairs)
}

#' Evaluate matching quality using known metabolite pairs
#' @param result Result from mass_combine
#' @return Numeric score between 0 and 1
evaluate_matches <- function(result) {
  # Get raw dataframes
  df1 <- raw_df(ms1(result))
  df2 <- raw_df(ms2(result))

  # Get all metabolites with multiplicities from original data
  true_metabs1 <- df1$Metabolite[!is.na(df1$Metabolite)]
  true_metabs2 <- df2$Metabolite[!is.na(df2$Metabolite)]

  # Count occurrences in each set
  count1 <- table(true_metabs1)
  count2 <- table(true_metabs2)

  # Find common elements and take the minimum count
  common_elements <- intersect(names(count1), names(count2))
  min_counts <- pmin(count1[common_elements], count2[common_elements])

  # Create the resulting set Z
  result_set <- rep(names(min_counts), min_counts)

  # Get the total number of elements in the resulting set
  total_number <- length(result_set)

  # Get all matches from the result
  matches <- result |>
    get_unique_matches()

  # Pre-compute column names
  ms1_metab_col <- paste0("Metabolite_", name(result@ms1))
  ms2_metab_col <- paste0("Metabolite_", name(result@ms2))

  # Create a named vector for renaming
  rename_cols <- c(metabolite_name_1 = ms1_metab_col,
                   metabolite_name_2 = ms2_metab_col)

  # Calculate valid matches (where metabolites are same)
  n_matches <- matches |>
    dplyr::rename(dplyr::any_of(rename_cols)) |>
    dplyr::filter(metabolite_name_1 == metabolite_name_2) |>
    nrow()
  # Calculate possible matches per metabolite
  return(n_matches / total_number)
}

#' Get unique 1-1 matches from mass_combine output
#'
#' @param ms_object MergedMSObject containing match results
#' @param pref Logical: whether to ensure every metabolite from dataset 1 gets a match (TRUE)
#'   or to optimize for overall match quality (FALSE). Default is FALSE.
#' @return Data frame containing only unique 1-1 matches, where each feature appears only once
#' @export
get_unique_matches <- function(ms_object, pref = FALSE) {
  # Get study names from ms_object
  study1_name <- name(ms_object@ms1)
  study2_name <- name(ms_object@ms2)

  # Create column names dynamically based on study names
  id_col1 <- paste0("Compound_ID_", study1_name)
  id_col2 <- paste0("Compound_ID_", study2_name)

  # Get matched data
  matched_data <- all_matched(ms_object)

  # Verify required columns exist
  required_cols <- c(id_col1, id_col2)
  missing_cols <- setdiff(required_cols, names(matched_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  # Check if score column exists, if not create simple ranking
  if (!"score" %in% names(matched_data)) {
    warning("No score column found. Using row order as score.")
    matched_data$score <- seq_len(nrow(matched_data))
  }

  if (pref) {
    # Preference-based matching: ensure every metabolite from dataset 1 gets a match
    result <- matched_data %>%
      # Remove rows where either ID is NA
      dplyr::filter(!is.na(.data[[id_col1]]), !is.na(.data[[id_col2]])) %>%
      # Sort by score (best scores first)
      dplyr::arrange(score) %>%
      # First, get best match for each ID in dataset 1
      dplyr::group_by(.data[[id_col1]]) %>%
      dplyr::slice_min(order_by = score,
                       n = 1,
                       with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      # Then handle any duplicates in dataset 2 by keeping best score
      dplyr::group_by(.data[[id_col2]]) %>%
      dplyr::slice_min(order_by = score,
                       n = 1,
                       with_ties = FALSE) %>%
      dplyr::ungroup()
  } else {
    # Original behavior: optimize for overall match quality
    result <- matched_data %>%
      # Remove rows where either ID is NA
      dplyr::filter(!is.na(.data[[id_col1]]), !is.na(.data[[id_col2]])) %>%
      # Sort by score (best scores first)
      dplyr::arrange(score) %>%
      # Keep only first occurrence of each ID in either column
      dplyr::filter(!duplicated(.data[[id_col1]]) &
                      !duplicated(.data[[id_col2]]))
  }

  # Add informative attributes
  attr(result, "n_input_rows") <- nrow(matched_data)
  attr(result, "n_output_rows") <- nrow(result)

  if (nrow(result) == 0) {
    warning("No unique matches found. Input had ",
            nrow(matched_data),
            " rows.")
  }

  return(result)
}
