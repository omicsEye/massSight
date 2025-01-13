#' @title Distribution Plot
#' @description Plot distributions of retention times and mass to charge ratios
#' of individual MS experiments
#'
#' @param ms_obj Either a `MSObject`, or `MergedMSObject`
#' @param subset Whether to plot all metabolites or only isolated metabolites.
#' Can either be "all" or "isolated"
#'
#' @return a scatter plot of the retention times vs mass to charge ratios with
#' marginal histograms
#' @export

distribution_plot <- function(ms_obj, subset = "all") {
  if (methods::is(ms_obj, "MSObject")) {
    if (subset == "all") {
      ms_df <- raw_df(ms_obj)
    } else if (subset == "isolated") {
      ms_df <- ms_obj@isolated
      if (nrow(ms_df) == 0) {
        stop("No isolated subset detected.")
      } else {
        stop("subset must be one of 'all' or 'isolated'")
      }
    }
    p_center <- ms_df |>
      ggplot2::ggplot(ggplot2::aes(.data$RT, .data$MZ)) +
      ggplot2::geom_point(alpha = .4) +
      ggplot2::theme_classic(base_size = 1.54)

    # Create top histogram (x density)
    p_top <- ms_df |>
      ggplot2::ggplot(ggplot2::aes(x = .data$RT)) +
      ggplot2::geom_histogram(fill = "grey70") +
      ggplot2::theme_classic(base_size = 1.54) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank()
      )

    # Create right histogram (y density)
    p_right <- ms_df |>
      ggplot2::ggplot(ggplot2::aes(x = .data$MZ)) +
      ggplot2::geom_histogram(fill = "grey70") +
      ggplot2::coord_flip() +
      ggplot2::theme_classic(base_size = 1.54) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank()
      )

    design <- "
AAAAAAAAAAAAAAAAAA####
BBBBBBBBBBBBBBBBBBCCCC
BBBBBBBBBBBBBBBBBBCCCC
BBBBBBBBBBBBBBBBBBCCCC
BBBBBBBBBBBBBBBBBBCCCC
"

    p_combined <- p_top + p_center + p_right +
      patchwork::plot_layout(design = design)

    return(p_combined)
  } else if (methods::is(ms_obj, "MergedMSObject")) {
    if (subset == "all") {
      ms_df1 <- raw_df(ms1(ms_obj))
      ms_df2 <- raw_df(ms2(ms_obj))
    } else if (subset == "isolated") {
      ms_df1 <- isolated(ms1(ms_obj))
      ms_df2 <- isolated(ms2(ms_obj))
    } else {
      stop("subset must be one of 'all' or 'isolated'")
    }
    p1 <- ms_df1 |>
      ggplot2::ggplot(ggplot2::aes(.data$RT, .data$MZ)) +
      ggplot2::geom_point(alpha = .4) +
      ggplot2::theme_classic(base_size = 1.54) +
      theme_omicsEye()
    p2 <- ms_df2 |>
      ggplot2::ggplot(ggplot2::aes(.data$RT, .data$MZ)) +
      ggplot2::geom_point(alpha = .4) +
      ggplot2::theme_classic(base_size = 1.54) +
      theme_omicsEye()
    p1 <-
      ggExtra::ggMarginal(p1, type = "histogram")
    p2 <-
      ggExtra::ggMarginal(p2, type = "histogram")
    return(cowplot::plot_grid(
      p1,
      p2,
      nrow = 1,
      labels = c("Dataset 1", "Dataset 2"),
      label_x = .5
    ))
  }
}
