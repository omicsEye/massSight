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

  if (smooth_method == "loess") {
    res_low <- stats::loess(delta_RT ~ RT, data = results)
    smooth_x_rt <- results$RT
    smooth_y_rt <- stats::predict(res_low, smooth_x_rt)
  } else if (smooth_method == "gam") {
    spline_func <- mgcv::gam(delta_RT ~ s(RT), data = results)
    smooth_x_rt <- results$RT
    smooth_y_rt <- stats::predict(spline_func, data.frame(RT = smooth_x_rt)) |>
      as.vector()
  } else if (smooth_method == "gaussian") {
    smooth_x_rt <- results$RT
    # TODO check for RBF Kernel
    message("Starting gaussian smoothing")
    gp <-
      GauPro::GauPro(smooth_x_rt,
        results$RT_2 - results$RT,
        parallel = F
      )
    smooth_y_rt <- gp$predict(smooth_x_rt)
    message("Finished gaussian smoothing")
  }
  smooth_x_rt_dropna <- smooth_x_rt |> stats::na.omit()
  smooth_y_rt_dropna <- smooth_y_rt |> stats::na.omit()
  if (length(smooth_x_rt_dropna) == 0) {
    stop(
      "There were not enough matches found to generate a predicted smoothing curve for your RT. Try increasing the range of your rt cutoffs (default -0.5 to +0.5), or increasing the 'smooth rt' value (default 0.1). Or if all else fails, try a manual/custom scaling (rt_custom)."
    )
  }

  suppressWarnings(f <- stats::approx(
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

  if (smooth_method == "loess") {
    mz_low <- stats::loess(delta_MZ ~ MZ, data = results)
    smooth_x_mz <- results$MZ
    smooth_y_mz <- stats::predict(mz_low, smooth_x_mz)
  } else if (smooth_method == "gam") {
    mz_gam <- mgcv::gam(delta_MZ ~ s(MZ), data = results)
    smooth_x_mz <- results$MZ
    smooth_y_mz <- stats::predict(mz_gam, data.frame(MZ = smooth_x_mz)) |>
      as.vector()
  } else if (smooth_method == "gaussian") {
    smooth_x_mz <- results$MZ
    message("Starting gaussian smoothing")
    gp <-
      GauPro::GauPro(smooth_x_mz,
        results$MZ_2 - results$MZ,
        parallel = FALSE
      )
    smooth_y_mz <- gp$predict(smooth_x_mz)
    message("Finished gaussian smoothing")
  } else {
    stop("Invalid smooth method")
  }
  smooth_x_mz_dropna <- smooth_x_mz |> stats::na.omit()
  smooth_y_mz_dropna <- smooth_y_mz |> stats::na.omit()
  if (length(smooth_x_mz_dropna) == 0) {
    stop(
      "There were not enough matches found to generate a predicted smoothing curve for your RT. Try increasing the range of your rt cutoffs (default -0.5 to +0.5), or increasing the 'smooth rt' value (default 0.1). Or if all else fails, try a manual/custom scaling (rt_custom)."
    )
  }

  suppressWarnings(f <- stats::approx(
    x = smooth_x_mz_dropna,
    y = smooth_y_mz_dropna,
    rule = 2,
    xout = smooth_x_mz
  ))

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
  smooth_method(align_ms_obj) <- smooth_method
  return(align_ms_obj)
}
