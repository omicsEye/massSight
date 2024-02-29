#' Final Plots
#'
#' @param merged_ms_obj A Merged MS Object, typically the result of `auto_combine`
#' @param rt_lim RT bounds
#' @param mz_lim MZ bounds
#'
#' @return `cowplot` grid of plots
#' @export
final_plots <-
  function(merged_ms_obj,
           rt_lim = c(-.5, .5),
           mz_lim = c(-15, 15)) {
    smooth_list <- merged_ms_obj@smooth_method
    pairs_base <- "Number of pairs:"
    pre_iso_pairs <-
      paste0(pairs_base, nrow(pre_iso_matched(merged_ms_obj)))
    iso_pairs <-
      paste0(pairs_base, nrow(iso_matched(merged_ms_obj)))
    all_pairs <-
      paste0(pairs_base, nrow(
        all_matched(merged_ms_obj) |> dplyr::filter(!is.na(Compound_ID_1) &
          !is.na(Compound_ID_2))
      ))

    rt_pre_iso <- pre_iso_matched(merged_ms_obj) |>
      ggplot2::ggplot(ggplot2::aes(x = .data$RT, y = .data$RT_2 - .data$RT)) +
      ggplot2::geom_point(
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(1.5),
        stroke = 0.05
      ) +
      ggplot2::labs(
        title = "All Matches",
        x = "RT 1",
        y = expression(Delta * "RT")
      ) +
      theme_omicsEye() +
      ggplot2::annotate(
        "label",
        x = Inf,
        y = -Inf,
        hjust = 1,
        vjust = 0,
        label = pre_iso_pairs,
        size = 8 / ggplot2::.pt
      )

    rt_pre_iso <-
      rt_pre_iso |> ggExtra::ggMarginal(
        type = "histogram",
        xparams = list(fill = "light gray", size = .25),
        yparams = list(fill = "light gray", size = .25)
      )

    rt_iso <- iso_matched(merged_ms_obj) |>
      ggplot2::ggplot(ggplot2::aes(
        x = .data$RT,
        y = .data$RT_2 - .data$RT
      )) +
      ggplot2::geom_point(
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(1.5),
        stroke = 0.05
      ) +
      ggplot2::geom_line(
        data = data.frame(
          "x" = smooth_list$rt_x,
          "y" = smooth_list$rt_y
        ),
        ggplot2::aes(x, y),
        col = "#AA9868"
      ) +
      ggplot2::annotate(
        "label",
        x = Inf,
        y = -Inf,
        hjust = 1,
        vjust = 0,
        label = iso_pairs,
        size = 8 / ggplot2::.pt
      ) +
      ggplot2::labs(
        title = "Isolated Matches",
        x = "RT 1",
        y = expression(Delta * "RT")
      ) +
      theme_omicsEye()

    rt_iso <-
      rt_iso |> ggExtra::ggMarginal(
        type = "histogram",
        xparams = list(fill = "light gray", size = .25),
        yparams = list(fill = "light gray", size = .25)
      )

    rt_all <- all_matched(merged_ms_obj) |>
      dplyr::filter(!is.na(Compound_ID_1) &
        !is.na(Compound_ID_2)) |>
      ggplot2::ggplot(ggplot2::aes(x = .data$RT_1, y = .data$RT_2_adj - .data$RT_1)) +
      ggplot2::geom_point(
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(1.5),
        stroke = 0.05
      ) +
      ggplot2::geom_hline(
        yintercept = 0,
        col = "#AA9868"
      ) +
      ggplot2::annotate(
        "label",
        x = Inf,
        y = -Inf,
        hjust = 1,
        vjust = 0,
        label = all_pairs,
        size = 8 / ggplot2::.pt
      ) +
      ggplot2::labs(
        title = "Scaled Matches",
        x = "RT 1",
        y = expression(Delta * "RT")
      ) +
      ggplot2::ylim(rt_lim) +
      theme_omicsEye()

    rt_all <-
      rt_all |> ggExtra::ggMarginal(
        type = "histogram",
        xparams = list(fill = "light gray", size = .25),
        yparams = list(fill = "light gray", size = .25)
      )

    mz_pre_iso <- pre_iso_matched(merged_ms_obj) |>
      dplyr::mutate(delta_MZ = ((.data$MZ_2 - .data$MZ) / ((
        .data$MZ + .data$MZ_2
      ) / 2) * 1e6)) |>
      ggplot2::ggplot(ggplot2::aes(x = .data$MZ, y = .data$delta_MZ)) +
      ggplot2::geom_point(
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(1.5),
        stroke = 0.05
      ) +
      ggplot2::labs(
        title = "All Matches",
        x = "MZ 1",
        y = expression(Delta * "MZ")
      ) +
      theme_omicsEye()

    mz_pre_iso <-
      mz_pre_iso |> ggExtra::ggMarginal(
        type = "histogram",
        xparams = list(fill = "light gray", size = .25),
        yparams = list(fill = "light gray", size = .25)
      )

    mz_iso <- iso_matched(merged_ms_obj) |>
      dplyr::mutate(delta_MZ = ((.data$MZ_2 - .data$MZ) / ((
        .data$MZ + .data$MZ_2
      ) / 2) * 1e6)) |>
      ggplot2::ggplot(ggplot2::aes(x = .data$MZ, y = .data$delta_MZ)) +
      ggplot2::geom_point(
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(1.5),
        stroke = 0.05
      ) +
      ggplot2::geom_line(
        data = data.frame(
          "x" = smooth_list$mz_x,
          "y" = smooth_list$mz_y
        ),
        ggplot2::aes(x, y),
        col = "#AA9868"
      ) +
      ggplot2::labs(
        title = "Isolated Matches",
        x = "MZ 1",
        y = expression(Delta * "MZ")
      ) +
      theme_omicsEye()

    mz_iso <-
      mz_iso |> ggExtra::ggMarginal(
        type = "histogram",
        xparams = list(fill = "light gray", size = .25),
        yparams = list(fill = "light gray", size = .25)
      )

    mz_all <- all_matched(merged_ms_obj) |>
      dplyr::filter(!is.na(Compound_ID_1) &
        !is.na(Compound_ID_2)) |>
      ggplot2::ggplot(ggplot2::aes(
        x = .data$MZ_1,
        y = (.data$MZ_2_adj - .data$MZ_1) / (.data$MZ_2_adj + MZ_1) / 2
          * 1e6
      )) +
      ggplot2::geom_point(
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(1.5),
        stroke = 0.05
      ) +
      ggplot2::geom_hline(
        yintercept = 0,
        col = "#AA9868"
      ) +
      ggplot2::labs(
        title = "Scaled Matches",
        x = "MZ 1",
        y = expression(Delta * "MZ")
      ) +
      ggplot2::ylim(mz_lim) +
      theme_omicsEye()

    mz_all <-
      mz_all |> ggExtra::ggMarginal(
        type = "histogram",
        xparams = list(fill = "light gray", size = .25),
        yparams = list(fill = "light gray", size = .25)
      )

    out <- cowplot::plot_grid(rt_pre_iso, rt_iso, rt_all,
      mz_pre_iso, mz_iso, mz_all,
      nrow = 2
    ) +
      ggplot2::theme(
        plot.background =
          ggplot2::element_rect(fill = "white")
      )

    return(out)
  }
