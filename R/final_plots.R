#' Final Plots
#'
#' @param merged_ms_obj A Merged MS Object, typically the result of `mass_combine`
#' @param rt_lim RT bounds
#' @param mz_lim MZ bounds
#'
#' @return `cowplot` grid of plots
#' @export
final_plots <-
  function(merged_ms_obj,
           rt_lim = c(-.5, .5),
           mz_lim = c(-15, 15)) {
    ci_name1 <- paste0("Compound_ID_", ms1(merged_ms_obj)@name)
    ci_name2 <- paste0("Compound_ID_", ms2(merged_ms_obj)@name)
    rt_name1 <- paste0("RT_", ms1(merged_ms_obj)@name)
    rt_name2 <- paste0("RT_", ms2(merged_ms_obj)@name)
    rt_adj_name2 <- paste0("RT_adj_", ms2(merged_ms_obj)@name)
    mz_name1 <- paste0("MZ_", ms1(merged_ms_obj)@name)
    mz_name2 <- paste0("MZ_", ms2(merged_ms_obj)@name)
    mz_adj_name2 <- paste0("MZ_adj_", ms2(merged_ms_obj)@name)
    smooth_list <- merged_ms_obj@smooth_method
    pairs_base <- "Number of pairs:"
    pre_iso_pairs <-
      paste0(pairs_base, nrow(pre_iso_matched(merged_ms_obj)))
    iso_pairs <-
      paste0(pairs_base, nrow(iso_matched(merged_ms_obj)))

    all_matched <- all_matched(merged_ms_obj) |>
      dplyr::filter(!is.na(.data[[ci_name1]]) &
        !is.na(.data[[ci_name2]])) |>
      dplyr::inner_join(
        adjusted_df(merged_ms_obj),
        by = stats::setNames("Compound_ID", ci_name2)
      )

    all_pairs <- paste0(pairs_base, nrow(all_matched))


    rt_iso <- iso_matched(merged_ms_obj) |>
      ggplot2::ggplot(ggplot2::aes(x = .data$RT_1, y = .data$RT_2 - .data$RT_1)) +
      ggplot2::geom_point(
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(1.5),
        stroke = 0.05
      ) +
      ggplot2::geom_line(
        data = data.frame("x" = smooth_list$rt_x, "y" = smooth_list$rt_y),
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

    rt_all <- all_matched |>
      ggplot2::ggplot(ggplot2::aes(x = .data[[rt_name1]], y = .data[[rt_adj_name2]] - .data[[rt_name1]])) +
      ggplot2::geom_point(
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(1.5),
        stroke = 0.05
      ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
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
      theme_omicsEye()

    rt_all <-
      rt_all |> ggExtra::ggMarginal(
        type = "histogram",
        xparams = list(fill = "light gray", size = .25),
        yparams = list(fill = "light gray", size = .25)
      )

    mz_iso <- iso_matched(merged_ms_obj) |>
      dplyr::mutate(
        delta_MZ = ((.data$MZ_2 - .data$MZ_1) / .data$MZ_1 * 1e6),
        mz_bin = factor(floor(MZ_1 / 100) * 100)
      ) |>
      ggplot2::ggplot(ggplot2::aes(x = .data$mz_bin, y = .data$delta_MZ)) +
      ggplot2::geom_violin(
        fill = "#033C5A",
        alpha = 0.5
      ) +
      ggplot2::geom_boxplot(
        width = 0.2,
        fill = "white",
        alpha = 0.5,
        outlier.shape = NA
      ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#AA9868") +
      ggplot2::labs(
        title = "Mass Error Distribution by m/z Bin",
        x = "m/z Bin",
        y = expression(Delta * "MZ (ppm)")
      ) +
      ggplot2::ylim(-5, 5) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      theme_omicsEye()

    mz_all <- all_matched |>
      ggplot2::ggplot(ggplot2::aes(
        x = .data[[mz_name1]], 
        y = (.data[[mz_adj_name2]] - .data[[mz_name1]]) / .data[[mz_name1]] * 1e6
      )) +
      ggplot2::geom_point(
        alpha = I(0.25),
        shape = 21,
        colour = "black",
        fill = "#033C5A",
        size = I(1.5),
        stroke = 0.05
      ) +
      ggplot2::geom_hline(yintercept = 0, col = "#AA9868") +
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

    out <- cowplot::plot_grid(rt_iso, rt_all, mz_iso, mz_all, nrow = 2) +
      ggplot2::theme(
        plot.background =
          ggplot2::element_rect(fill = "white")
      )

    return(out)
  }
