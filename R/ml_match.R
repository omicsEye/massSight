#' ML Match
#'
#' @description Trains, fits, and uses a machine learning model on known
#' metabolite data to predict unknown metabolite pair matches.
#'
#' @param ms1 A `massSight` object representing the results of a preprocessed
#' LC-MS experiment
#' @param ms2 A `massSight` object representing the results of a second
#' preprcoessed LC-MS experiment
#' @param mz_thresh `numeric` Mass to Charge threshold. Used to limit potential matches
#' between metabolites.
#' @param rt_thresh `numeric` Retention Time threshold. Used to limit potential matches
#' between metabolites
#' @param seed Seed value for reproducibility
#'
#' @note This function requires semi-annotated data (some metabolites must be
#' named)
#'
#' @return A dataframe consisting of predicted metabolite pairs
#' @export
#'
#' @examples
#' \dontrun{
#' ml_match(ms1, ms2, mz_thresh = 15, rt_thresh = .5, seed = 2)
#' }
ml_match <-
  function(ms1,
           ms2,
           mz_thresh = 15,
           rt_thresh = 1,
           prob_thresh = .5,
           seed = 72) {
    known_metabolites <- get_shared_metabolites(ms1, ms2)
    training_data <-
      get_training_data(known_metabolites, ms2, mz_thresh, rt_thresh)
    fitted_model <- fit_model(training_data, seed)
    pred_fmt_data <-
      create_pred_data(ms1, ms2, mz_thresh, rt_thresh)
    pred_fmt_data$prob_matched <-
      stats::predict(fitted_model, pred_fmt_data, type = "prob")$matched
    pred_fmt_data <- pred_fmt_data |>
      dplyr::filter(.data$prob_matched > prob_thresh)
    ml_match_list <- list(
      "matched_data" = pred_fmt_data,
      "model" = fitted_model,
      "training_data" = training_data
    )
    return(ml_match_list)
  }

get_shared_metabolites <- function(ms1, ms2) {
  ms1_df <- raw_df(ms1)
  ms2_df <- raw_df(ms2)
  ms1_known <- ms1_df |>
    dplyr::filter(.data$Metabolite != "" &
                    .data$Compound_ID != .data$Metabolite)
  ms1_known_scaled <- ms1_df |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), scale)) |>
    dplyr::filter(.data$Metabolite != "" &
                    .data$Compound_ID != .data$Metabolite)
  ms2_known <- ms2_df |>
    dplyr::filter(.data$Metabolite != "" &
                    .data$Compound_ID != .data$Metabolite)
  known <-
    dplyr::inner_join(ms1_known, ms2_known, by = "Metabolite") |>
    dplyr::select(.data$MZ.x,
                  .data$MZ.y,
                  .data$RT.x,
                  .data$RT.y,
                  .data$Metabolite) |>
    dplyr::rename(
      MZ_1 = .data$MZ.x,
      MZ_2 = .data$MZ.y,
      RT_1 = .data$RT.x,
      RT_2 = .data$RT.y,
      Metabolite_1 = .data$Metabolite
    ) |>
    dplyr::mutate(Class = "matched", Metabolite_2 = .data$Metabolite_1)
  out <- list(
    "known" = known,
    "ms1_known" = ms1_known,
    "ms2_known" = ms2_known,
    "ms1_known_scaled" = ms1_known_scaled
  )
  return(out)
}

get_training_data <-
  function(shared,
           ms2,
           mz_thresh = 30,
           rt_thresh = 2) {
    known <- shared$known
    ms1_known <- shared$ms1_known
    ms2_known <- shared$ms2_known
    ms1_known_scaled <- shared$ms1_known_scaled
    for (i in 1:nrow(known)) {
      ms1_sample <- ms1_known[i, ]
      ms1_sample_scaled <- ms1_known_scaled[i, ]
      ms2_sample <- ms2@raw_df |>
        dplyr::filter(.data$Metabolite != ms1_sample_scaled[1, "Metabolite"])
      ms2_sample_scaled <- ms2@raw_df |>
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric), scale)) |>
        dplyr::filter(.data$Metabolite != ms1_sample_scaled[1, "Metabolite"])
      distance <- sqrt((ms1_sample$MZ - ms2_sample$MZ) ^ 2 + (ms1_sample$RT - ms2_sample$RT) ^
                         2)
      ms2_sample <- ms2_sample[which.min(distance), ]
      row <- data.frame(
        MZ_1 = ms1_sample$MZ,
        MZ_2 = ms2_sample$MZ,
        RT_1 = ms1_sample$RT,
        RT_2 = ms2_sample$RT,
        Metabolite_1 = ms1_sample$Metabolite,
        Metabolite_2 = ms2_sample$Metabolite,
        Class = "unmatched"
      )
      known <- dplyr::bind_rows(known, row)
    }
    known <- known |>
      dplyr::mutate(
        Class = as.factor(.data$Class),
        delta_RT = .data$RT_1 - .data$RT_2,
        delta_MZ = .data$MZ_1 - .data$MZ_2
      )
    return(known)
  }

fit_model <- function(known, seed) {
  set.seed(seed)

  rf <- caret::train(
    Class ~ delta_RT + delta_MZ + RT_1 + RT_2 + MZ_1 + MZ_2,
    data = known,
    method = "rf",
    ntree = 501,
    tuneGrid = data.frame(mtry = 1:6),
    trControl = caret::trainControl(method = "cv", number = 10)
  )

  return(rf)
}

create_pred_data <-
  function(ms1,
           ms2,
           mz_thresh = 15,
           rt_thresh = .5) {
    ms1_df_unknown <- ms1 |>
      raw_df() |>
      dplyr::rename(
        Compound_ID_1 = .data$Compound_ID,
        MZ_1 = .data$MZ,
        RT_1 = .data$RT,
        Metabolite_1 = .data$Metabolite
      ) |>
      dplyr::select(-.data$Intensity)
    ms2_df_unknown <- ms2 |>
      raw_df() |>
      dplyr::rename(
        Compound_ID_2 = .data$Compound_ID,
        MZ_2 = .data$MZ,
        RT_2 = .data$RT,
        Metabolite_2 = .data$Metabolite
      ) |>
      dplyr::select(-.data$Intensity)

    combined <- tidyr::crossing(ms1_df_unknown, ms2_df_unknown) |>
      dplyr::mutate(delta_MZ = .data$MZ_1 - .data$MZ_2,
                    delta_RT = .data$RT_1 - .data$RT_2) |>
      dplyr::filter(abs(.data$delta_MZ) < mz_thresh &
                      abs(.data$delta_RT) < rt_thresh)
    return(combined)
  }

plot_decision_boundary <- function(matched_list) {
  rt1_min <- min(matched_list$matched_data$RT_1)
  rt1_max <- max(matched_list$matched_data$RT_1)
  rt2_min <- min(matched_list$matched_data$RT_2)
  rt2_max <- max(matched_list$matched_data$RT_2)
  mz1_min <- min(matched_list$matched_data$MZ_1)
  mz1_max <- max(matched_list$matched_data$MZ_1)
  mz2_min <- min(matched_list$matched_data$MZ_2)
  mz2_max <- max(matched_list$matched_data$MZ_2)
  delta_mz_min <- min(matched_list$matched_data$delta_MZ)
  delta_mz_max <- max(matched_list$matched_data$delta_MZ)
  delta_rt_min <- min(matched_list$matched_data$delta_RT)
  delta_rt_max <- max(matched_list$matched_data$delta_RT)
  plot_data <- expand.grid(
    "delta_MZ" = seq(delta_mz_min, delta_mz_max, length.out = 20),
    "delta_RT" = seq(delta_rt_min, delta_rt_max, length.out = 20),
    "RT_1" = seq(rt1_min, rt1_max, length.out = 20),
    "RT_2" = seq(rt2_min, rt2_max, length.out = 20),
    "MZ_1" = seq(mz1_min, mz1_max, length.out = 20),
    "MZ_2" = seq(mz2_min, mz2_max, length.out = 20)
  )
  plot_data$pred_prob <-
    stats::predict(matched_list$model, plot_data, type = "prob")$matched
  plot_data |>
    dplyr::mutate(match = dplyr::case_when(pred_prob > .5 ~ "yes", T ~ "no")) |>
    ggplot2::ggplot(ggplot2::aes(x = delta_RT, y = delta_MZ, fill = match)) +
    ggplot2::geom_tile()
}
