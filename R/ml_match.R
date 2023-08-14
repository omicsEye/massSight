#' ML Match
#'
#' @description Trains, fits, and uses a machine learning model on known
#' metabolite data to predict unknown metabolite pair matches.
#'
#' @param ms1 A `massSight` object representing the results of a preprocessed
#' LC-MS experiment
#' @param ms2 A `massSight` object representing the results of a second
#' preprcoessed LC-MS experiment
#' @param mz_thresh Mass to Charge threshold. Used to limit potential matches
#' between metabolites.
#' @param rt_thresh Retention Time threshold. Used to limit potential matches
#' between metabolites
#' @param seed Seed value for reproducibility
#'
#' @return A dataframe consisting of predicted metabolite pairs
#' @export
ml_match <-
  function(ms1,
           ms2,
           mz_thresh = 15,
           rt_thresh = .5,
           seed = 72) {
    known_metabolites <- get_shared_metabolites(ms1, ms2)
    training_data <-
      get_training_data(known_metabolites, mz_thresh, rt_thresh)
    fitted_model <- fit_model(training_data, seed)
    pred_fmt_data <-
      create_pred_data(ms1, ms2, mz_thresh, rt_thresh)
    pred_fmt_data$prob_matched <-
      predict(fitted_model, pred_fmt_data, type = "prob")$matched
    pred_fmt_data <- pred_fmt_data |>
      dplyr::filter(prob_matched > .95)
    return(pred_fmt_data)
  }

get_shared_metabolites <- function(ms1, ms2) {
  ms1_df <- raw_df(ms1)
  ms2_df <- raw_df(ms2)
  ms1_known <- ms1_df |>
    dplyr::filter(Metabolite != "" &
                    Compound_ID != Metabolite)
  ms2_known <- ms2_df |>
    dplyr::filter(Metabolite != "" &
                    Compound_ID != Metabolite)
  known <-
    dplyr::inner_join(ms1_known, ms2_known, by = "Metabolite") |>
    dplyr::select(MZ.x, MZ.y, RT.x, RT.y, Metabolite) |>
    dplyr::rename(
      MZ_1 = MZ.x,
      MZ_2 = MZ.y,
      RT_1 = RT.x,
      RT_2 = RT.y,
      Metabolite_1 = Metabolite
    ) |>
    dplyr::mutate(Class = "matched",
                  Metabolite_2 = Metabolite_1)
  out <- list("known" = known,
              "ms1_known" = ms1_known,
              "ms2_known" = ms2_known)
  return(out)
}

get_training_data <-
  function(shared,
           mz_thresh = 15,
           rt_thresh = .5) {
    known <- shared$known
    ms1_known <- shared$ms1_known
    ms2_known <- shared$ms2_known
    for (i in 1:nrow(known)) {
      ms1_sample <- ms1_known |>
        dplyr::slice_sample(n = 1)
      ms2_sample <-
        ms2_sample <- ms2_known |>
        dplyr::filter(
          Metabolite != ms1_sample[1, "Metabolite"] &
            abs(ms1_sample$MZ - MZ) < mz_thresh &
            abs(ms1_sample$RT - RT) < rt_thresh
        ) |>
        dplyr::slice_sample(n = 1)
      if (nrow(ms2_sample) == 0) {
        i <- i - 1
        next
      }
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
        Class = as.factor(Class),
        delta_RT = RT_1 - RT_2,
        delta_MZ = MZ_1 - MZ_2
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
      dplyr::rename(Compound_ID_1 = Compound_ID,
                    MZ_1 = MZ,
                    RT_1 = RT,
                    Metabolite_1 = Metabolite) |>
      dplyr::select(-Intensity)
    ms2_df_unknown <- ms2 |>
      raw_df() |>
      dplyr::rename(Compound_ID_2 = Compound_ID,
                    MZ_2 = MZ,
                    RT_2 = RT,
                    Metabolite_2 = Metabolite) |>
      dplyr::select(-Intensity)

    combined <- tidyr::crossing(ms1_df_unknown, ms2_df_unknown) |>
      dplyr::mutate(delta_MZ = MZ_1 - MZ_2,
                    delta_RT = RT_1 - RT_2) |>
      dplyr::filter(abs(delta_MZ) < mz_thresh &
                      abs(delta_RT) < rt_thresh)
    return(combined)
  }
