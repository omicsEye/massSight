#' Final Plots
#'
#' @param merged_ms_obj A Merged MS Object, typically the result of `auto_align`
#' @param rt_lim RT bounds
#' @param mz_lim MZ bounds
#'
#' @return cowplot grid of plots
#' @export
#' @examples
#' data("hp1_ms")
#' data("hp2_ms")
#' merged_ms_obj <- auto_align(hp1_ms, hp2_ms)
#' final_plots(merged_ms_obj)
final_plots <-
  function(merged_ms_obj,
           rt_lim = c(-.5, .5),
           mz_lim = c(-15, 15)) {
    rt_pre_iso <- pre_iso_matched(merged_ms_obj) |>
      ggplot2::ggplot(ggplot2::aes(x = RT, y = RT_2 - RT)) +
      ggplot2::geom_point(
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(1.5),
        stroke = 0.05
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        alpha = .3,
        col = "#AA9868"
      ) +
      ggplot2::labs(
        title = "All Matches",
        x = "RT 1",
        y = expression(Delta * "RT")
      ) +
      ggplot2::theme_classic() +
      theme_omicsEye()

    rt_iso <- iso_matched(merged_ms_obj) |>
      ggplot2::ggplot(ggplot2::aes(
        x = RT,
        y = RT_2 - RT
      )) +
      ggplot2::geom_point(
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(1.5),
        stroke = 0.05
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        alpha = .3,
        col = "#AA9868"
      ) +
      ggplot2::labs(
        title = "Isolated Matches",
        x = "RT 1",
        y = expression(Delta * "RT")
      ) +
      ggplot2::theme_classic() +
      theme_omicsEye()

    rt_all <- all_matched(merged_ms_obj) |>
      dplyr::mutate(scaled_rts = adjusted_df(merged_ms_obj)$rt_2_adj -
        all_matched(merged_ms_obj)$RT) |>
      ggplot2::ggplot(ggplot2::aes(x = RT, y = scaled_rts)) +
      ggplot2::geom_point(
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(1.5),
        stroke = 0.05
      ) +
      ggplot2::ylim(rt_lim) +
      ggplot2::geom_smooth(
        method = "loess",
        alpha = .3,
        col = "#AA9868"
      ) +
      ggplot2::labs(
        title = "Scaled Matches",
        x = "RT 1",
        y = expression(Delta * "RT")
      ) +
      ggplot2::theme_classic() +
      theme_omicsEye()

    mz_pre_iso <- pre_iso_matched(merged_ms_obj) |>
      dplyr::mutate(delta_MZ = ((MZ_2 - MZ) / ((MZ + MZ_2) / 2) * 1e6)) |>
      ggplot2::ggplot(ggplot2::aes(x = MZ, y = delta_MZ)) +
      ggplot2::geom_point(
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(1.5),
        stroke = 0.05
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        alpha = .3,
        col = "#AA9868"
      ) +
      ggplot2::labs(
        title = "All Matches",
        x = "MZ 1",
        y = expression(Delta * "MZ")
      ) +
      ggplot2::theme_classic() +
      theme_omicsEye()

    mz_iso <- iso_matched(merged_ms_obj) |>
      dplyr::mutate(delta_MZ = ((MZ_2 - MZ) / ((MZ + MZ_2) / 2) * 1e6)) |>
      ggplot2::ggplot(ggplot2::aes(x = MZ, y = delta_MZ)) +
      ggplot2::geom_point(
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(1.5),
        stroke = 0.05
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        alpha = .3,
        col = "#AA9868"
      ) +
      ggplot2::labs(
        title = "Isolated Matches",
        x = "MZ 1",
        y = expression(Delta * "MZ")
      ) +
      ggplot2::theme_classic() +
      theme_omicsEye()

    mz_all <- all_matched(merged_ms_obj) |>
      dplyr::mutate(scaled_mz = (
        adjusted_df(merged_ms_obj)$mz_2_adj - all_matched(merged_ms_obj)$MZ
      ) / ((
        adjusted_df(merged_ms_obj)$mz_2_adj + all_matched(merged_ms_obj)$MZ
      ) / 2) * 1e6) |>
      ggplot2::ggplot(ggplot2::aes(x = MZ, y = scaled_mz)) +
      ggplot2::geom_point(
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(1.5),
        stroke = 0.05
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        alpha = .3,
        col = "#AA9868"
      ) +
      ggplot2::labs(
        title = "Scaled Matches",
        x = "MZ 1",
        y = expression(Delta * "MZ")
      ) +
      ggplot2::theme_classic() +
      ggplot2::ylim(mz_lim) +
      theme_omicsEye()

    out <- cowplot::plot_grid(rt_pre_iso, rt_iso, rt_all,
      mz_pre_iso, mz_iso, mz_all,
      nrow = 2
    )

    return(out)
  }