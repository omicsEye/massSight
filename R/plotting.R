qc_plots <-
  function(results,
           smooth = NULL,
           rt_lim = c(-.5, .5),
           mz_lim = c(-15, 15)) {
    plot_res <- dplyr::tibble(
      rt_x = results$RT,
      raw_rts = results$RT_2 - results$RT,
      scaled_rts = results$srt - results$RT,
      mz_x = results$MZ,
      raw_mz = ((results$MZ_2 - results$MZ) / ((
        results$MZ + results$MZ_2
      ) / 2) * 1e6),
      scaled_mz = ((results$smz - results$MZ) / ((
        results$MZ + results$smz
      ) / 2) * 1e6),
      raw_int1 = results$Intensity,
      raw_int2 = results$Intensity_2,
      scaled_int1 = log10(results$Intensity),
      scaled_int2 = log10(results$Intensity_2),
      int_scaled_diff = log10(results$sintensity) - log10(results$Intensity),
      int_diff = results$Intensity - results$Intensity_2
    )

    plot_scale <- dplyr::tibble(
      smooth_rt_x = smooth$RT,
      smooth_rt_y = smooth$smooth_rt,
      smooth_mz_x = smooth$MZ,
      smooth_mz_y = smooth$smooth_mz / smooth$MZ * 1e6,
    )

    rt_plot_raw <-
      ggplot2::ggplot() +
      ggplot2::geom_point(
        data = plot_res,
        ggplot2::aes(x = rt_x, y = raw_rts),
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(1.5),
        stroke = 0.05
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        data = plot_scale,
        ggplot2::aes(x = smooth_rt_x, y = smooth_rt_y),
        alpha = .8,
        col = "#AA9868"
      ) +
      ggplot2::labs(title = "Raw",
                    x = "RT 1",
                    y = expression(Delta * "RT")) +
      ggplot2::theme_classic(base_size = 1.54) +
      theme_omicsEye()

    rt_plot_scaled <-
      ggplot2::ggplot() +
      ggplot2::geom_point(
        data = plot_res,
        ggplot2::aes(x = rt_x, y = scaled_rts),
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(2),
        stroke = 0.05
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        data = plot_res,
        ggplot2::aes(x = rt_x, y = scaled_rts),
        alpha = .8,
        col = "#AA9868"
      ) +
      ggplot2::labs(title = "Scaled",
                    x = "RT 1",
                    y = expression(Delta * "RT")) +
      ggplot2::theme_classic(base_size = 1.54) +
      theme_omicsEye()

    mz_plot_raw <-
      ggplot2::ggplot() +
      ggplot2::geom_point(
        data = plot_res,
        ggplot2::aes(x = mz_x, y = raw_mz),
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(2),
        stroke = 0.05
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        data = plot_scale,
        ggplot2::aes(x = smooth_mz_x, y = smooth_mz_y),
        alpha = 1,
        color = "#AA9868"
      ) +
      ggplot2::labs(title = "Raw",
                    x = "MZ 1",
                    y = expression(Delta * "MZ")) +
      ggplot2::theme_classic(base_size = 1.54) +
      theme_omicsEye()

    mz_plot_scaled <-
      ggplot2::ggplot() +
      ggplot2::geom_point(
        data = plot_res,
        ggplot2::aes(x = mz_x, y = scaled_mz),
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(2),
        stroke = 0.05
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        data = plot_res,
        ggplot2::aes(x = mz_x, y = scaled_mz),
        alpha = .8,
        col = "#AA9868"
      ) +
      ggplot2::labs(title = "Scaled",
                    x = "MZ 1",
                    y = expression(Delta * "MZ")) +
      ggplot2::theme_classic(base_size = 1.54) +
      theme_omicsEye()

    int_plot_raw <-
      ggplot2::ggplot() +
      ggplot2::geom_point(
        data = plot_res,
        ggplot2::aes(x = raw_int1, y = int_diff),
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(2),
        stroke = 0.05
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        data = plot_res,
        ggplot2::aes(x = raw_int1, y = int_diff),
        alpha = .8,
        col = "#AA9868"
      ) +
      ggplot2::labs(title = "Raw",
                    x = "Int 1",
                    y = expression(Delta * "RT")) +
      ggplot2::theme_classic(base_size = 1.54) +
      theme_omicsEye()

    int_plot_scaled <-
      ggplot2::ggplot() +
      ggplot2::geom_point(
        data = plot_res,
        ggplot2::aes(x = scaled_int1, y = int_scaled_diff),
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(2),
        stroke = 0.05
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        data = plot_res,
        ggplot2::aes(x = scaled_int1, y = int_scaled_diff),
        alpha = .8,
        col = "#AA9868"
      ) +
      ggplot2::labs(title = "Scaled",
                    x = "log10(Int 1)",
                    y = expression(Delta * "log10(RT)")) +
      ggplot2::theme_classic(base_size = 1.54) +
      theme_omicsEye()


    plots <- cowplot::plot_grid(
      rt_plot_raw,
      rt_plot_scaled,
      mz_plot_raw,
      mz_plot_scaled,
      int_plot_raw,
      int_plot_scaled,
      ncol = 2
    )
    return(plots)
  }

#' @export
final_plots <-
  function(merged_ms_obj,
           rt_lim = c(-.5, .5),
           mz_lim = c(-15, 15)) {
    results <- all_matched(merged_ms_obj)
    smooth <- iso_matched(merged_ms_obj)
    scaled_results <- adjusted_df(merged_ms_obj)
      both_found_key <- !is.na(results$RT) & !is.na(results$RT_2)
    results_df <-
      results[both_found_key, c("RT", "MZ", "Intensity", "RT_2", "MZ_2", "Intensity_2")]
    results_df$srt <- scaled_results$rt_2_adj
    results_df$smz <- scaled_results$mz_2_adj
    results_df$sintensity <- scaled_results$int_2_adj
    return(qc_plots(results_df, smooth, rt_lim, mz_lim))
  }
