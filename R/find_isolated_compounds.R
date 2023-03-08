find_isolated_compounds <-
  function(df1,
           df2,
           rt_lower,
           rt_upper,
           mz_lower,
           mz_upper,
           rt_smooth,
           mz_smooth,
           minimum_intensity,
           threshold,
           per,
           rt_iso_threshold,
           mz_iso_threshold,
           match_method,
           smooth_method,
           custom_rt = "None") {
    disallowed_names <-
      "" # ['RT_1', 'RT_2', 'MZ_1', 'MZ_2', 'Intensity_1', 'Intensity_2']

    if (!is.null(disallowed_names)) {
      if (disallowed_names %in% colnames(df1)) {
        names <- dplyr::intersect(disallowed_names, df1)
        return(
          paste(
            "You have a column labeled",
            names,
            "in file 1. Please change this to something else."
          )
        )
      }
    }
    if (disallowed_names %in% colnames(df2)) {
      names <- dplyr::intersect(disallowed_names, df2)
      return(
        paste(
          "You have a column labeled",
          names,
          "in file 2. Please change this to something else."
        )
      )
    }

    if (match_method == "unsupervised") {
      if (threshold == "manual") {
        vec_1 <-
          get_vectors_manual(df1, rt_iso_threshold, mz_iso_threshold)
        vec_2 <-
          get_vectors_manual(df2, rt_iso_threshold, mz_iso_threshold)
      } else if (threshold == "auto") {
        vec_1 <-
          get_vectors_auto(df1, per)
      }
      vec_1 <- df1 |>
        dplyr::filter(Compound_ID %in% vec_1) |>
        dplyr::select("RT", "MZ", "Intensity", "Metabolite", "Compound_ID")
      vec_2 <- df2 |>
        dplyr::filter(Compound_ID %in% vec_2) |>
        dplyr::select("RT", "MZ", "Intensity", "Metabolite", "Compound_ID")

      results <-
        align_isolated_compounds(vec_1, vec_2, rt_lower, rt_upper, mz_lower, mz_upper)
    } else if (match_method == "supervised") {
      stopifnot("Metabolite" %in% colnames(df1) &
                  "Metabolite" %in% colnames(df2))
      vec_1 <- df1 |>
        dplyr::rename(df1 = Compound_ID) |>
        dplyr::filter(Metabolite != "")
      vec_2 <- df2 |>
        dplyr::rename(
          RT_2 = RT,
          MZ_2 = MZ,
          Intensity_2 = Intensity,
          df2 = Compound_ID
        ) |>
        dplyr::filter(Metabolite != "")
      results <- vec_1 |>
        dplyr::inner_join(vec_2, by = c("Metabolite"))
    } else {
      stop("Valid arguments for `match_method` are 'supervised' or 'unsupervised'")
    }
    if (results$RT_2 |>
        na.omit() |>
        length() == 0) {
      stop(
        "Couldn't find any potential matches between the datasets. Check your RT and mz values in your input files. They may mislabeled or reversed."
      )
    }
    results <- results |>
      dplyr::arrange(RT)

    if (typeof(custom_rt) == "list") {
      smooth_x_rt <- custom_rt[1]
      smooth_y_rt <- custom_rt[2]
    } else if (smooth_method == "lowess") {
      res_low <- lowess(x = results$RT,
                        y = results$RT_2 - results$RT,
                        f = rt_smooth)
      smooth_x_rt <- res_low$x
      smooth_y_rt <- res_low$y
    } else if (smooth_method == "spline") {
      results <- results |>
        dplyr::distinct(RT, .keep_all = TRUE)
      smooth_x_rt <- results$RT
      smooth_y_rt <-
        splinefun(smooth_x_rt, results$RT_2 - results$RT)
      smooth_y_rt <- smooth_y_rt(smooth_x_rt)
    } else if (smooth_method == "gaussian") {
      smooth_x_rt <- results$RT
      # TODO check for RBF Kernel
      gp <- GauPro::GauPro(smooth_x_rt, results$RT_2 - results$RT)
      smooth_y_rt <- gp$predict(smooth_x_rt)
    }
    smooth_x_rt_dropna <- smooth_x_rt |> na.omit()
    smooth_y_rt_dropna <- smooth_y_rt |> na.omit()
    if (length(smooth_x_rt_dropna) == 0) {
      return(
        "There were not enough matches found to generate a predicted smoothing curve for your RT. Try increasing the range of your rt cutoffs (default -0.5 to +0.5), or increasing the 'smooth rt' value (default 0.1). Or if all else fails, try a manual/custom scaling (rt_custom)."
      )
    }

    suppressWarnings(f <- approx(
      x = smooth_x_rt_dropna,
      y = smooth_y_rt_dropna,
      xout = smooth_x_rt,
      rule = 2
    ))
    smooth_x_rt <- results$RT
    smooth_y_rt <- f$y
    scaled_rts <-
      scale_smooth(df2$RT, smooth_x_rt + smooth_y_rt, smooth_y_rt)
    scaled_rts_res <-
      scale_smooth(results$RT,
                   smooth_x_rt + smooth_y_rt,
                   smooth_y_rt)

    results <- results |>
      dplyr::mutate(smooth_rt = smooth_y_rt,
                    srt = scaled_rts_res)
    # results <- results |>
    #   dplyr::mutate(smooth_rt = smooth_y_rt,
    #                 srt = RT_2 - smooth_rt)
    #
    # # message("Scaling rts for final match")
    # scaled_rts <-
    #   scale_smooth(df2$RT, smooth_x_rt + smooth_y_rt, smooth_y_rt)
    #
    # scaled_vector_rts <- data.frame("RT" = scaled_rts)
    # scaled_vector_rts$Metabolite <- df2$Metabolite
    # scaled_vector_rts <- scaled_vector_rts |>
    #   filter("Metabolite" %in% results$df2)
    # # TODO ask ali about m2align line 520
    #
    # results$metabolite <- scaled_vector_rts$Metabolite
    # results$srt <- scaled_vector_rts$RT


    ## scale mzs ---------------------------------------------------------------
    results <- results |>
      dplyr::arrange("MZ")

    mz_df <- results |>
      dplyr::select(MZ, MZ_2)

    if (smooth_method == "lowess") {
      mz_low <- lowess(x = mz_df$MZ,
                       y = mz_df$MZ_2 - mz_df$MZ,
                       f = mz_smooth)
      smooth_x_mz <- mz_low$x
      smooth_y_mz <- mz_low$y
    } else if (smooth_method == "spline") {
      mz_df <- mz_df |>
        dplyr::distinct(MZ, .keep_all = TRUE)
      results <- results |>
        dplyr::distinct(MZ, .keep_all = TRUE)
      smooth_x_mz <- mz_df$MZ
      smooth_y_mz <-
        splinefun(smooth_x_mz, mz_df$MZ_2 - mz_df$MZ)
      smooth_y_mz <- smooth_y_mz(smooth_x_mz)
    } else if (smooth_method == "gaussian") {
      smooth_x_mz <- mz_df$MZ
      # TODO check for RBF Kernel
      gp <- GauPro::GauPro(smooth_x_mz, mz_df$MZ_2 - mz_df$MZ)
      smooth_y_mz <- gp$predict(smooth_x_mz)
    }
    smooth_x_mz_dropna <- smooth_x_mz |> na.omit()
    smooth_y_mz_dropna <- smooth_y_mz |> na.omit()
    if (length(smooth_x_mz_dropna) == 0) {
      return(
        "There were not enough matches found to generate a predicted smoothing curve for your RT. Try increasing the range of your rt cutoffs (default -0.5 to +0.5), or increasing the 'smooth rt' value (default 0.1). Or if all else fails, try a manual/custom scaling (rt_custom)."
      )
    }

    suppressWarnings(f <- approx(
      x = smooth_x_mz_dropna,
      y = smooth_y_mz_dropna,
      rule = 2,
      xout = smooth_x_mz
    ))

    smooth_x_mz <- results$MZ
    smooth_y_mz <- f$y
    scaled_mzs <-
      scale_smooth(df2$MZ, smooth_x_mz + smooth_y_mz, smooth_y_mz)
    scaled_mzs_res <- scale_smooth(results$MZ,
                                   smooth_x_mz + smooth_y_mz,
                                   smooth_y_mz)

    results <- results |>
      dplyr::mutate(smooth_mz = smooth_y_mz,
                    smz = scaled_mzs_res)
    # results <- results |>
    #   dplyr::arrange(MZ) |>
    #   dplyr::mutate(smooth_mz = smooth_y_mz,
    #                 smz = MZ_2 - smooth_mz)

    # message("Scaling mzs for final match")
    # scaled_mzs <-
    #   scale_smooth(df2$MZ, smooth_x_mz + smooth_y_mz, smooth_y_mz)
    # scaled_vector_mzs <- data.frame("mz" = scaled_mzs)
    # scaled_vector_mzs <- scaled_vector_mzs[results$df2,]
    # results <- results |>
    #   dplyr::arrange(mz)
    # scaled_mzs_df <- data.frame("smz" = scaled_mzs,
    #                          "df2" = df2$Compound_ID)
    # results <- dplyr::left_join(results, scaled_mzs_df, by = "df2")
    # results$smz <- scaled_vector_mzs$mz


    # scale intensities -------------------------------------------------------
    temp_df1_int <- log10(results$Intensity)
    temp_df2_int <- log10(results$Intensity_2)

    ## find slope for linear adjustment of log-intensity parameters
    intensity_parameters <- scale_intensity_parameters(temp_df1_int,
                                                       temp_df2_int,
                                                       min_int = minimum_intensity)

    ## scale potential matches
    scaled_vector_intensity <-
      scale_intensity(temp_df2_int, intensity_parameters)
    scaled_vector_intensity <- 10 ^ scaled_vector_intensity
    results$sintensity <- scaled_vector_intensity

    # scale full results
    log_df2 <- log10(df2$Intensity)
    scaled_intensity <-
      scale_intensity(log_df2, intensity_parameters)
    scaled_intensity <- 10 ^ scaled_intensity

    dev_out <- get_cutoffs(
      results |>
        dplyr::select("RT", "MZ", "Intensity"),
      data.frame(
        "RT" = results$srt,
        "MZ" = results$smz,
        "Intensity" = results$sintensity
      ),
      c("RT", "MZ", "Intensity")
    )
    deviations <- dev_out$cutoffs
    outliers <- dev_out$outliers
    scaled_values <- data.frame("RT" = scaled_rts,
                                "MZ" = scaled_mzs,
                                "Intensity" = scaled_intensity)
    return(
      list(
        "results" = results,
        "scaled_values" = scaled_values,
        "deviations" = deviations
      )
    )
  }
