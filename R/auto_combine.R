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

    # check if ms1 and ms2 are massSight objects
    if (!is(ms1, "MSObject") || !is(ms2, "MSObject")) {
      stop("ms1 and ms2 must be massSight objects")
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
        weights = weights
      )

    if (!is.null(log)) {
      log_parameters(log, log_params, log_r, log_date, align_obj, ms1, ms2, time_start)
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

dedup <- function(cols, item) {
  item_locations <- which(cols == item)
  for (i in seq_along(item_locations)) {
    if (i != 1) {
      cols[item_locations[i]] <- paste0(item, "_", i)
    }
  }
  return(cols)
}

mad_based_outlier <- function(points, thresh = 5) {
  diff <- sqrt((points - stats::median(points, na.rm = TRUE))^2)
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
  function(data_int1,
           data_int2,
           int_col = "Intensity",
           min_int = 0) {
    min_int <- log10(min_int)
    int_key <- (data_int1 > min_int) & (data_int2 > min_int)
    fit <- mean(data_int2[int_key] / data_int1[int_key])
    return(fit)
  }

calculate_mean <- function(values, w, method = "arit") {
  if (sum(values) == 0) {
    return(0)
  }
  if (method == "whm") {
    num <- replace(w, values < 0, 0)
    denom <- replace(w, values < 0, 0)
    denom <- replace(denom, values > 0, denom / values)
    w_h_mean <- sum(num) / sum(denom)
    nonzeroes <- sum(values != 0)
    mean_value <- w_h_mean / (nonzeroes + 1)
  } else if (method == "arit") {
    mean_value <- sum(w * values) / sum(w)
  }
  return(mean_value)
}

rms <- function(a, b, std) {
  rt_difference <- a$RT - b$RT
  if (std[1] == 0) {
    rt_contribution <- abs(rt_difference)
  } else {
    rt_contribution <- (rt_difference / std[1])**2
  }

  # convert to mzs to ppm to score
  ppm_difference <- (a$MZ - b$MZ) * 1e6 / mean(c(a$MZ, b$MZ))
  if (std[2] == 0) {
    mz_contribution <- abs(ppm_difference)
  } else {
    mz_contribution <- (ppm_difference / std[2])^2
  }

  score <- mean(c(rt_contribution, mz_contribution))
  return(score)
}

scale_smooth <- function(query_values, smooth_x, smooth_y) {
  suppressWarnings(f <- stats::approx(
    x = smooth_x,
    y = smooth_y,
    xout = query_values,
    rule = 2
  ))

  query_out <- query_values - f$y
  # smooth_min <- smooth_y[which.min(smooth_x)]
  # smooth_max <- smooth_y[which.max(smooth_x)]
  #
  # query_values <- na.omit(query_values)
  #
  # query_values[query_values <= min(smooth_x)] <-
  #   query_values[query_values <= min(smooth_x)] - smooth_min
  # query_values[query_values >= max(smooth_x)] <-
  #   query_values[query_values >= max(smooth_x)] - smooth_max
  #
  # smooth_df <- data.frame("x" = smooth_x,
  #                         "y" = smooth_y) |>
  #   dplyr::arrange(x)
  #
  # query_values[query_values > min(smooth_x) &
  #                query_values < max(smooth_x)] <-
  #   pracma::interp1(smooth_df$x,
  #                   smooth_df$y,
  #                   query_values[query_values > min(smooth_x) &
  #                                  query_values < max(smooth_x)])

  return(query_out)
}

align_isolated_compounds <-
  function(align_ms_obj,
           match_method,
           rt_minus = -.5,
           rt_plus = .5,
           mz_minus = -15,
           mz_plus = 15,
           keep = FALSE) {
    df1 <- isolated(ms1(align_ms_obj))
    df2 <- isolated(ms2(align_ms_obj))
    if (match_method == "unsupervised") {
      if ("Intensity" %in% names(df1)) {
        pb <-
          progress::progress_bar$new(
            format = "Matching isolated features from datasets [:bar] :percent :eta",
            total = nrow(df1),
            clear = F
          )
        results <- data.frame(
          "df1" = character(),
          "RT" = numeric(),
          "MZ" = numeric(),
          "Intensity" = numeric(),
          "df2" = character(),
          "RT_2" = numeric(),
          "MZ_2" = numeric(),
          "Intensity_2" = numeric()
        )
        for (row in 1:nrow(df1)) {
          pb$tick()
          df2_filter <- df2 |>
            dplyr::filter(
              .data$RT > df1$RT[row] + rt_minus &
                .data$RT < df1$RT[row] + rt_plus &
                .data$MZ > df1$MZ[row] + mz_minus * df1$MZ[row] / 1e6 &
                .data$MZ < df1$MZ[row] + mz_plus * df1$MZ[row] / 1e6
            )
          if (nrow(df2_filter) > 0) {
            for (row_2 in 1:nrow(df2_filter)) {
              res_add <- data.frame(
                df1 = df1$Compound_ID[row],
                RT = df1$RT[row],
                MZ = df1$MZ[row],
                Intensity = df1$Intensity[row],
                df2 = df2_filter$Compound_ID[row_2],
                RT_2 = df2_filter$RT[row_2],
                MZ_2 = df2_filter$MZ[row_2],
                Intensity_2 = df2_filter$Intensity[row_2]
              )
              results <- results |>
                rbind(res_add)
            }
          }
        }
      } else {
        pb <-
          progress::progress_bar$new(
            format = "Matching all features from datasets [:bar] :percent :eta",
            total = nrow(df1),
            clear = F
          )
        results <- data.frame(
          "df1" = character(),
          "RT" = numeric(),
          "MZ" = numeric(),
          "df2" = character(),
          "RT_2" = numeric(),
          "MZ_2" = numeric()
        )

        for (row in 1:nrow(df1)) {
          pb$tick()
          df2_filter <- df2 |>
            dplyr::filter(
              .data$RT > (df1[row, "RT"] + rt_minus),
              .data$RT < (df1[row, "RT"] + rt_plus),
              .data$MZ > (df1[row, "MZ"] + mz_minus / 1e6),
              .data$MZ < (df1[row, "MZ"] + mz_plus / 1e6)
            )

          if (nrow(df2_filter) > 0) {
            for (row_2 in 1:nrow(df2_filter)) {
              results <- results |>
                dplyr::bind_rows(
                  data.frame(
                    "df1" = df1[row, "Compound_ID"],
                    "RT" = df1[row, "RT"],
                    "MZ" = df1[row, "MZ"],
                    "df2" = df2_filter[row_2, "Compound_ID"],
                    "RT_2" = df2_filter[row_2, "RT"],
                    "MZ_2" = df2_filter[row_2, "MZ"]
                  )
                )
            }
          }
        }
      }
    } else if (match_method == "supervised") {
      stopifnot("Metabolite" %in% colnames(df1) &
        "Metabolite" %in% colnames(df2))
      vec_1 <- df1 |>
        dplyr::rename(df1 = .data$Compound_ID) |>
        dplyr::filter(.data$Metabolite != "")
      vec_2 <- df2 |>
        dplyr::rename(
          RT_2 = .data$RT,
          MZ_2 = .data$MZ,
          Intensity_2 = .data$Intensity,
          df2 = .data$Compound_ID
        ) |>
        dplyr::filter(.data$Metabolite != "")
      results <- vec_1 |>
        dplyr::inner_join(vec_2, by = c("Metabolite"))
    }
    iso_matched(align_ms_obj) <- results
    return(align_ms_obj)
  }

final_results <-
  function(align_ms_obj,
           weights = c(1, 1, 1)) {
    study1_name <- name(ms1(align_ms_obj))
    study2_name <- name(ms2(align_ms_obj))
    Compound_ID_1 <- paste("Compound_ID", study1_name, sep = "_")
    Compound_ID_2 <- paste("Compound_ID", study2_name, sep = "_")
    Metabolite_1 <- paste("Metabolite", study1_name, sep = "_")
    Metabolite_2 <- paste("Metabolite", study2_name, sep = "_")
    RT_1 <- paste("RT", study1_name, sep = "_")
    RT_2 <- paste("RT", study2_name, sep = "_")
    MZ_1 <- paste("MZ", study1_name, sep = "_")
    MZ_2 <- paste("MZ", study2_name, sep = "_")
    Intensity_1 <- paste("Intensity", study1_name, sep = "_")
    Intensity_2 <- paste("Intensity", study2_name, sep = "_")
    df1 <- raw_df(ms1(align_ms_obj))
    df2 <- raw_df(ms2(align_ms_obj))
    scaled_df <- scaled_values(align_ms_obj)
    stds <- cutoffs(align_ms_obj)
    df2$RT_2_adj <- scaled_df$RT
    df2$MZ_2_adj <- scaled_df$MZ
    df2_adj <- df2
    df2_adj$RT <- scaled_df$RT
    df2_adj$MZ <- scaled_df$MZ
    df2_adj$Intensity <- scaled_df$Intensity
    df1_for_align <- df1 |>
      dplyr::select(dplyr::any_of(c("Compound_ID", "RT", "MZ", "Intensity")))
    df2_for_align <- df2_adj |>
      dplyr::select(dplyr::any_of(c("Compound_ID", "RT", "MZ", "Intensity")))

    best_hits_df1 <- c()
    best_hits_found <- c()
    features_not_aligned <- c()
    pb <-
      progress::progress_bar$new(
        format = "Aligning datasets [:bar] :percent :eta",
        total = nrow(df1_for_align),
        clear = FALSE
      )
    for (i in seq_len(nrow(df1_for_align))) {
      best_match <-
        find_closest_match(df1_for_align[i, ], df2_for_align, stds)
      if (!is.null(best_match)) {
        pb$tick()
        best_reverse_match <-
          find_closest_match(
            df2_for_align |>
              dplyr::filter(.data$Compound_ID == best_match),
            df1_for_align,
            stds
          )
      } else {
        features_not_aligned <-
          c(features_not_aligned, df1_for_align[i, "Compound_ID"])
        pb$tick()
        next
      }

      if (df1_for_align[i, "Compound_ID"] %in% best_reverse_match) {
        best_hits_df1 <- c(best_hits_df1, best_match)
        best_hits_found <-
          c(best_hits_found, rep(df1_for_align[[i, "Compound_ID"]], length(best_match)))
      }
    }

    match_df <- data.frame(df1 = best_hits_found, df2 = best_hits_df1)
    df2_raw <- df2

    if (nrow(metadata(ms1(align_ms_obj))) > 0) {
      df1 <- df1 |>
        merge(metadata(ms1(align_ms_obj)), by = "Compound_ID")
    }
    if (nrow(metadata(ms2(align_ms_obj))) > 0) {
      df2 <- df2 |>
        merge(metadata(ms2(align_ms_obj)), by = "Compound_ID")
    }
    df <-
      merge(df1,
        match_df,
        by.x = "Compound_ID",
        by.y = "df1",
        all = TRUE
      ) |>
      merge(df2,
        by.x = "df2",
        by.y = "Compound_ID",
        all = TRUE
      )

    df <- df |>
      dplyr::rename_with(
        ~ ifelse(
          .x %in% names(df),
          c(
            paste("Compound_ID", study1_name, sep = "_"),
            paste("Compound_ID", study2_name, sep = "_"),
            paste("Metabolite", study1_name, sep = "_"),
            paste("Metabolite", study2_name, sep = "_"),
            paste("RT", study1_name, sep = "_"),
            paste("RT", study2_name, sep = "_"),
            paste("MZ", study1_name, sep = "_"),
            paste("MZ", study2_name, sep = "_"),
            paste("Intensity", study1_name, sep = "_"),
            paste("Intensity", study2_name, sep = "_")
          ),
          .x
        ),
        .cols = c(
          "Compound_ID",
          "df2",
          "Metabolite.x",
          "Metabolite.y",
          "RT.x",
          "RT.y",
          "MZ.x",
          "MZ.y",
          "Intensity.x",
          "Intensity.y"
        )
      ) |>
      dplyr::mutate(
        rep_Compound_ID = dplyr::case_when(
          !is.na(get(Compound_ID_1)) ~
            get(Compound_ID_1),
          is.na(get(Compound_ID_1)) &
            !is.na(get(Compound_ID_2)) ~ get(Compound_ID_2),
          TRUE ~ NA
        ),
        rep_RT = dplyr::case_when(
          !is.na(get(RT_1)) ~ get(RT_1),
          is.na(get(RT_1)) & !is.na(get(RT_2)) ~ get(RT_2),
          TRUE ~ NA
        ),
        rep_MZ = dplyr::case_when(
          !is.na(get(MZ_1)) ~ get(MZ_1),
          is.na(get(MZ_1)) & !is.na(get(MZ_2)) ~ get(MZ_2),
          TRUE ~ NA
        ),
        rep_Intensity = dplyr::case_when(
          !is.na(get(Intensity_1)) ~ get(Intensity_1),
          is.na(get(Intensity_1)) &
            !is.na(get(Intensity_2)) ~ get(Intensity_2),
          TRUE ~ NA
        ),
        rep_Metabolite = dplyr::case_when(
          !is.na(get(Metabolite_1)) ~ get(Metabolite_1),
          is.na(get(Metabolite_1)) &
            !is.na(get(Metabolite_2)) ~ get(Metabolite_2),
          TRUE ~ NA
        )
      ) |>
      dplyr::select(
        c(
          "rep_Compound_ID",
          "rep_RT",
          "rep_MZ",
          "rep_Intensity",
          "rep_Metabolite",
          !!dplyr::sym(Compound_ID_1),
          !!dplyr::sym(Compound_ID_2),
          !!dplyr::sym(Metabolite_1),
          !!dplyr::sym(Metabolite_2),
          !!dplyr::sym(RT_1),
          !!dplyr::sym(RT_2),
          !!dplyr::sym(MZ_1),
          !!dplyr::sym(MZ_2),
          !!dplyr::sym(Intensity_1),
          !!dplyr::sym(Intensity_2),
          dplyr::everything(),
          -dplyr::contains("_adj")
        )
      )

    df <- df %>%
      dplyr::group_by(rep_Compound_ID) |>
      dplyr::mutate(dup_count = dplyr::row_number()) |>
      dplyr::ungroup() |>
      dplyr::mutate(rep_Compound_ID = ifelse(
        dup_count > 1,
        paste0(rep_Compound_ID, "_", dup_count),
        rep_Compound_ID
      )) |>
      dplyr::select(-dup_count)

    all_matched(align_ms_obj) <- df
    adjusted_df(align_ms_obj) <- df2_adj

    return(align_ms_obj)
  }

find_all_matches <- function(ref, query, rt_threshold, mz_threshold) {
  matches <- data.frame()
  ref <- ref |>
    dplyr::rename(
      RT_1 = .data$RT,
      MZ_1 = .data$MZ,
      Compound_ID_1 = .data$Compound_ID
    )
  query <- query |>
    dplyr::rename(
      RT_2 = .data$RT,
      MZ_2 = .data$MZ,
      Compound_ID_2 = .data$Compound_ID
    )
  for (i in 1:nrow(ref)) {
    rt <- ref$RT_1[i]
    mz <- ref$MZ_1[i]
    rt_lower <- rt - rt_threshold
    rt_upper <- rt + rt_threshold
    mz_lower <- mz - mz_threshold
    mz_upper <- mz + mz_threshold
    query_matches <- query |>
      dplyr::filter(.data$RT_2 >= rt_lower &
        .data$RT_2 <= rt_upper &
        .data$MZ_2 >= mz_lower &
        .data$MZ_2 <= mz_upper)
    if (nrow(query_matches) > 0) {
      matches <- matches |>
        dplyr::bind_rows(dplyr::bind_cols(ref[i, ], query_matches))
    }
  }
  return(matches)
}

find_closest_match <-
  function(query,
           ref,
           stds) {
    ref_index <- ref$Compound_ID

    rt_hits <- (ref$RT <= (query$RT + .5)) &
      (ref$RT >= (query$RT - .5))

    mz_hits <- (ref$MZ < (query$MZ + .1)) &
      (ref$MZ > (query$MZ - .1))

    combined_hits <- rt_hits & mz_hits

    if (!(TRUE %in% (combined_hits))) {
      return(NULL)
    }

    hits <- ref |>
      dplyr::select(dplyr::any_of(c("Compound_ID", "RT", "MZ", "Intensity"))) |>
      dplyr::filter(rt_hits & mz_hits)
    hits_index <- ref_index[rt_hits & mz_hits]
    hits_results <- c()
    for (i in seq_len(nrow(hits))) {
      score <- rms(
        query,
        hits[i, ],
        stds
      )
      hits_results <- c(hits_results, score)
    }
    return(hits_index[hits_results == min(hits_results)])
  }

get_cutoffs <-
  function(df1,
           df2,
           has_int = TRUE) {
    data_rt <- df2$RT - df1$RT
    data_mz <- (df2$MZ - df1$MZ) / df1$MZ * 1e6

    if (has_int) {
      data_int <- log10(df2$Intensity) - log10(df1$Intensity)
      not_outliers <- !mad_based_outlier(data_rt) &
        !mad_based_outlier(data_mz) &
        !mad_based_outlier(data_int)
      data_int <- replace(data_int, data_int %in% c(Inf, -Inf), NA)
      cutoffs <- c(
        stats::sd(data_rt[not_outliers]),
        stats::sd(data_mz[not_outliers]),
        stats::sd(data_int[not_outliers])
      )
    } else {
      not_outliers <- !mad_based_outlier(data_rt) &
        !mad_based_outlier(data_mz)
      cutoffs <- c(
        stats::sd(data_rt[not_outliers]),
        stats::sd(data_mz[not_outliers])
      )
    }

    # TODO Fix below
    # rt_outliers <- df1[mad_based_outlier(data_rt)] |>
    #   rownames()
    # mz_outliers <- df1[mad_based_outlier(data_mz)] |>
    #   rownames()
    # int_outliers <- df1[mad_based_outlier(data_int)] |>
    #   rownames()
    # outliers <- c(rt_outliers, mz_outliers, int_outliers)
    outliers <- "tmp"
    return(list(
      "cutoffs" = cutoffs,
      "outliers" = outliers
    ))
  }

smooth_drift <- function(align_ms_obj,
                         smooth_method,
                         minimum_int) {
  df1 <- align_ms_obj |>
    ms1() |>
    raw_df()
  df2 <- align_ms_obj |>
    ms2() |>
    raw_df()
  results <- iso_matched(align_ms_obj)
  results <- results |>
    dplyr::arrange(.data$RT) |>
    dplyr::mutate(delta_RT = .data$RT_2 - .data$RT)

  if (smooth_method == "gam") {
    spline_func <- mgcv::gam(delta_RT ~ s(RT),
      method = "REML",
      data = results
    )
    smooth_x_rt <- results$RT
    smooth_y_rt <-
      stats::predict(spline_func, data.frame(RT = smooth_x_rt)) |>
      as.vector()
    smooth_method(align_ms_obj)[["rt_x"]] <- smooth_x_rt
    smooth_method(align_ms_obj)[["rt_y"]] <- smooth_y_rt
  } else if (smooth_method == "gaussian") {
    smooth_x_rt <- results$RT
    # TODO check for RBF Kernel
    message("Starting gaussian smoothing")
    gp <-
      GauPro::gpkm(
        smooth_x_rt,
        results$RT_2 - results$RT,
        kernel = "matern52",
        parallel = FALSE,
        normalize = TRUE,
        verbose = 0
      )
    smooth_y_rt <- gp$pred(smooth_x_rt)
    message("Finished gaussian smoothing")
    smooth_method(align_ms_obj)[["rt_x"]] <- smooth_x_rt
    smooth_method(align_ms_obj)[["rt_y"]] <- smooth_y_rt
  }
  smooth_x_rt_dropna <- smooth_x_rt |> stats::na.omit()
  smooth_y_rt_dropna <- smooth_y_rt |> stats::na.omit()
  if (length(smooth_x_rt_dropna) == 0) {
    stop(
      "There were not enough matches found to generate a predicted smoothing curve for your RT. Try increasing the range of your rt cutoffs (default -0.5 to +0.5), or increasing the 'smooth rt' value (default 0.1). Or if all else fails, try a manual/custom scaling (rt_custom)."
    )
  }

  suppressWarnings(
    f <- stats::approx(
      x = smooth_x_rt_dropna,
      y = smooth_y_rt_dropna,
      xout = smooth_x_rt,
      rule = 2
    )
  )
  smooth_x_rt <- results$RT
  smooth_y_rt <- f$y
  scaled_rts <-
    scale_smooth(df2$RT, smooth_x_rt + smooth_y_rt, smooth_y_rt)
  scaled_rts_res <-
    scale_smooth(
      results$RT,
      smooth_x_rt + smooth_y_rt,
      smooth_y_rt
    )

  results <- results |>
    dplyr::mutate(
      smooth_rt = smooth_y_rt,
      srt = scaled_rts_res
    )


  ## scale mzs ---------------------------------------------------------------
  results <- results |>
    dplyr::arrange(.data$MZ) |>
    dplyr::mutate(delta_MZ = .data$MZ_2 - .data$MZ)

  if (smooth_method == "gam") {
    mz_gam <- mgcv::gam(delta_MZ ~ s(MZ),
      method = "REML",
      data = results
    )
    smooth_x_mz <- results$MZ
    smooth_y_mz <-
      stats::predict(mz_gam, data.frame(MZ = smooth_x_mz)) |>
      as.vector()
    smooth_method(align_ms_obj)[["mz_x"]] <- smooth_x_mz
    smooth_method(align_ms_obj)[["mz_y"]] <- smooth_y_mz
  } else if (smooth_method == "gaussian") {
    smooth_x_mz <- results$MZ
    message("Starting gaussian smoothing")
    gp <-
      GauPro::gpkm(
        smooth_x_mz,
        results$MZ_2 - results$MZ,
        kernel = "matern52",
        parallel = FALSE,
        normalize = TRUE,
        verbose = 0
      )
    smooth_y_mz <- gp$predict(smooth_x_mz)
    message("Finished gaussian smoothing")
    smooth_method(align_ms_obj)[["mz_x"]] <- smooth_x_mz
    smooth_method(align_ms_obj)[["mz_y"]] <- smooth_y_mz
  }
  smooth_x_mz_dropna <- smooth_x_mz |> stats::na.omit()
  smooth_y_mz_dropna <- smooth_y_mz |> stats::na.omit()
  if (length(smooth_x_mz_dropna) == 0) {
    stop(
      "There were not enough matches found to generate a predicted smoothing curve for your RT. Try increasing the range of your rt cutoffs (default -0.5 to +0.5), or increasing the 'smooth rt' value (default 0.1). Or if all else fails, try a manual/custom scaling (rt_custom)."
    )
  }

  suppressWarnings(
    f <- stats::approx(
      x = smooth_x_mz_dropna,
      y = smooth_y_mz_dropna,
      rule = 2,
      xout = smooth_x_mz
    )
  )

  smooth_x_mz <- results$MZ
  smooth_y_mz <- f$y
  scaled_mzs <-
    scale_smooth(df2$MZ, smooth_x_mz + smooth_y_mz, smooth_y_mz)
  scaled_mzs_res <- scale_smooth(
    results$MZ,
    smooth_x_mz + smooth_y_mz,
    smooth_y_mz
  )

  results <- results |>
    dplyr::mutate(
      smooth_mz = smooth_y_mz,
      smz = scaled_mzs_res
    )

  # scale intensities -------------------------------------------------------
  if ("Intensity" %in% names(results)) {
    temp_df1_int <- log10(results$Intensity)
    temp_df2_int <- log10(results$Intensity_2)

    ## find slope for linear adjustment of log-intensity parameters
    intensity_parameters <- scale_intensity_parameters(temp_df1_int,
      temp_df2_int,
      min_int = minimum_int
    )

    ## scale potential matches
    scaled_vector_intensity <-
      scale_intensity(temp_df2_int, intensity_parameters)
    scaled_vector_intensity <- 10^scaled_vector_intensity
    results$sintensity <- scaled_vector_intensity

    # scale full results
    log_df2 <- log10(df2$Intensity)
    scaled_intensity <-
      scale_intensity(log_df2, intensity_parameters)
    scaled_intensity <- 10^scaled_intensity
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
    scaled_df <- data.frame(
      "RT" = results$srt,
      "MZ" = results$smz
    )
    scaled_values <- data.frame(
      "RT" = scaled_rts,
      "MZ" = scaled_mzs
    )
  }

  dev_out <- get_cutoffs(
    df1 = results |>
      dplyr::select(dplyr::any_of(c(
        "RT", "MZ", "Intensity"
      ))),
    df2 = scaled_df,
    has_int = ("Intensity" %in% names(results))
  )

  deviations <- dev_out$cutoffs
  outliers <- dev_out$outliers

  iso_matched(align_ms_obj) <- results
  scaled_values(align_ms_obj) <- scaled_values
  cutoffs(align_ms_obj) <- deviations
  return(align_ms_obj)
}
