final_results <-
  function(align_ms_obj,
           keep_features = c(FALSE, FALSE),
           weights = c(1, 1, 1)) {
    df1 <- raw_df(ms1(align_ms_obj))
    df2 <- raw_df(ms2(align_ms_obj))
    scaled_df <- scaled_values(align_ms_obj)
    stds <- cutoffs(align_ms_obj)
    df2$RT_2_adj <- scaled_df$RT
    df2$MZ_2_adj <- scaled_df$MZ
    df2_adj <- df2
    df2_adj$RT <- scaled_df$RT
    df2_adj$MZ <- scaled_df$MZ
    df2_adj$Intensity <- scaled_df$Intensity
    df1_for_align <- df1 |>
      dplyr::select(dplyr::any_of(c("Compound_ID", "RT", "MZ", "Intensity")))
    df2_for_align <- df2_adj |>
      dplyr::select(dplyr::any_of(c("Compound_ID", "RT", "MZ", "Intensity")))

    best_hits_df1 <- c()
    best_hits_found <- c()
    features_not_aligned <- c()
    pb <-
      progress::progress_bar$new(
        format = "Aligning datasets [:bar] :percent :eta",
        total = nrow(df1_for_align),
        clear = F
      )
    for (i in 1:nrow(df1_for_align)) {
      best_match <-
        find_closest_match(
          df1_for_align[i, ],
          df2_for_align,
          stds
        )
      if (!is.null(best_match)) {
        pb$tick()
        best_reverse_match <-
          find_closest_match(
            df2_for_align |>
              dplyr::filter(.data$Compound_ID == best_match),
            df1_for_align,
            stds
          )
      } else {
        features_not_aligned <-
          c(features_not_aligned, df1_for_align[i, "Compound_ID"])
        pb$tick()
        next
      }

      if (df1_for_align[i, "Compound_ID"] %in% best_reverse_match) {
        best_hits_df1 <- c(best_hits_df1, best_match)
        best_hits_found <-
          c(best_hits_found, rep(df1_for_align[[i, "Compound_ID"]], length(best_match)))
      }
    }

    match_df <- data.frame(
      df1 = best_hits_found,
      df2 = best_hits_df1
    )
    df2_raw <- df2

    df <-
      merge(df1,
        match_df,
        by.x = "Compound_ID",
        by.y = "df1",
        all = TRUE
      ) |>
      merge(df2,
        by.x = "df2",
        by.y = "Compound_ID",
        all = TRUE
      )

    df <- df %>%
      dplyr::rename_with(
        ~ ifelse(
          .x %in% names(df),
          c(
            "Compound_ID_1",
            "Compound_ID_2",
            "Metabolite_1",
            "Metabolite_2",
            "RT_1",
            "RT_2",
            "MZ_1",
            "MZ_2",
            "Intensity_1",
            "Intensity_2"
          ),
          .x
        ),
        .cols = c(
          "Compound_ID",
          "df2",
          "Metabolite.x",
          "Metabolite.y",
          "RT.x",
          "RT.y",
          "MZ.x",
          "MZ.y",
          "Intensity.x",
          "Intensity.y"
        )
      )

    all_matched(align_ms_obj) <- df
    adjusted_df(align_ms_obj) <- df2_adj
    metadata_1 <- metadata(ms1(align_ms_obj))
    metadata_2 <- metadata(ms2(align_ms_obj))
    if (nrow(metadata_1) > 0) {
      metadata_1 <-
        metadata_1 |>
        dplyr::semi_join(results_df,
          by = dplyr::join_by(Compound_ID == df1)
        )
    }
    if (nrow(metadata_2) > 0) {
      metadata_2 <-
        metadata_2 |> dplyr::semi_join(results_df, by = dplyr::join_by(Compound_ID == df2))
    }
    if (nrow(metadata_1) > 0 && nrow(metadata_2) > 0) {
      metadata(align_ms_obj) <- cbind(metadata_1, metadata_2)
    }
    return(align_ms_obj)
  }
