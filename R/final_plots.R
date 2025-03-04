#' Create publication-ready visualization of mass spectrometry alignment results
#'
#' @param merged_ms_obj A Merged MS Object from `mass_combine`
#' @param rt_lim RT bounds for plotting, default c(-.5, .5)
#' @param mz_lim MZ bounds for plotting, default c(-15, 15)
#' @param point_params List of visual parameters for points (alpha, size, etc.)
#' @param smooth_color Color for smoothing lines, default "#AA9868"
#' @param point_color Color for points, default "#033C5A"
#' @param label_size Text size for annotations, default 8
#'
#' @return A `cowplot` grid of alignment visualization plots
#' @export
final_plots <- function(merged_ms_obj,
                        rt_lim = c(-.5, .5),
                        mz_lim = c(-15, 15),
                        point_params = list(alpha = 0.25,
                                            size = 1.5,
                                            stroke = 0.05),
                        smooth_color = "#AA9868",
                        point_color = "#033C5A",
                        label_size = 8) {
  # Input validation
  if (!inherits(merged_ms_obj, "MergedMSObject")) {
    stop("merged_ms_obj must be a MergedMSObject")
  }

  # Helper function to create column names
  get_names <- function(ms_obj) {
    list(
      ci = paste0("Compound_ID_", ms_obj@name),
      rt = paste0("RT_", ms_obj@name),
      rt_adj = paste0("RT_adj_", ms_obj@name),
      mz = paste0("MZ_", ms_obj@name),
      mz_adj = paste0("MZ_adj_", ms_obj@name)
    )
  }

  # Get column names
  names1 <- get_names(ms1(merged_ms_obj))
  names2 <- get_names(ms2(merged_ms_obj))

  # Create base plotting function
  create_scatter <- function(data, x, y, title, xlabel, ylabel) {
    ggplot2::ggplot(data, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
      ggplot2::geom_point(
        alpha = I(point_params$alpha),
        shape = 21,
        colour = "black",
        fill = point_color,
        size = I(point_params$size),
        stroke = point_params$stroke
      ) +
      ggplot2::labs(title = title, x = xlabel, y = ylabel) +
      theme_omicsEye()
  }

  ci_name1 <- names1$ci
  ci_name2 <- names2$ci
  rt_name1 <- names1$rt
  rt_name2 <- names2$rt
  rt_adj_name2 <- names2$rt_adj
  mz_name1 <- names1$mz
  mz_name2 <- names2$mz
  mz_adj_name2 <- names2$mz_adj
  smooth_list <- merged_ms_obj@smooth_method
  pairs_base <- "Number of pairs:"
  pre_iso_pairs <-
    paste0(pairs_base, nrow(pre_iso_matched(merged_ms_obj)))
  iso_pairs <-
    paste0(pairs_base, nrow(iso_matched(merged_ms_obj)))

  all_matched <- merged_ms_obj |>
    get_unique_matches()

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
    ggplot2::labs(title = "Isolated Matches",
                  x = "RT 1",
                  y = expression(Delta * "RT")) +
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
    ggplot2::labs(title = "Scaled Matches",
                  x = "RT 1",
                  y = expression(Delta * "RT")) +
    theme_omicsEye()

  rt_all <-
    rt_all |> ggExtra::ggMarginal(
      type = "histogram",
      xparams = list(fill = "light gray", size = .25),
      yparams = list(fill = "light gray", size = .25)
    )

  mz_iso <- iso_matched(merged_ms_obj) |>
    dplyr::mutate(delta_MZ = ((.data$MZ_2 - .data$MZ_1) / .data$MZ_1 * 1e6)) |>
    ggplot2::ggplot(ggplot2::aes(x = .data$MZ_1, y = .data$delta_MZ)) +
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
        "x" = merged_ms_obj@smooth_method$mz_x,
        "y" = merged_ms_obj@smooth_method$mz_y
      ),
      ggplot2::aes(x, y),
      col = "#AA9868"
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(title = "Mass Error Distribution",
                  x = "m/z",
                  y = expression(Delta * "MZ (ppm)")) +
    ggplot2::ylim(mz_lim) +
    theme_omicsEye()

  mz_all <- all_matched |>
    ggplot2::ggplot(ggplot2::aes(x = .data[[mz_name1]], y = (.data[[mz_adj_name2]] - .data[[mz_name1]]) / .data[[mz_name1]] * 1e6)) +
    ggplot2::geom_point(
      alpha = I(0.25),
      shape = 21,
      colour = "black",
      fill = "#033C5A",
      size = I(1.5),
      stroke = 0.05
    ) +
    ggplot2::geom_hline(yintercept = 0, col = "#AA9868") +
    ggplot2::labs(title = "Scaled Matches",
                  x = "MZ 1",
                  y = expression(Delta * "MZ")) +
    ggplot2::ylim(mz_lim) +
    theme_omicsEye()

  mz_all <-
    mz_all |> ggExtra::ggMarginal(
      type = "histogram",
      xparams = list(fill = "light gray", size = .25),
      yparams = list(fill = "light gray", size = .25)
    )

  out <- cowplot::plot_grid(rt_iso, rt_all, mz_iso, mz_all, nrow = 2) +
    ggplot2::theme(plot.background =
                     ggplot2::element_rect(fill = "white"))

  return(out)
}
