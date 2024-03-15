final_results <-
  function(align_ms_obj,
           keep_features = c(FALSE, FALSE),
           weights = c(1, 1, 1)) {
    study1_name <- name(ms1(align_ms_obj))
    study2_name <- name(ms2(align_ms_obj))
    Compound_ID_1 <- paste("Compound_ID", study1_name, sep = "_")
    Compound_ID_2 <- paste("Compound_ID", study2_name, sep = "_")
    Metabolite_1 <- paste("Metabolite", study1_name, sep = "_")
    Metabolite_2 <- paste("Metabolite", study2_name, sep = "_")
    RT_1 <- paste("RT", study1_name, sep = "_")
    RT_2 <- paste("RT", study2_name, sep = "_")
    MZ_1 <- paste("MZ", study1_name, sep = "_")
    MZ_2 <- paste("MZ", study2_name, sep = "_")
    Intensity_1 <- paste("Intensity", study1_name, sep = "_")
    Intensity_2 <- paste("Intensity", study2_name, sep = "_")
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

    df1_full <- df1 |>
      merge(metadata(ms1(align_ms_obj)), by = "Compound_ID")
    df2_full <- df2 |>
      merge(metadata(ms2(align_ms_obj)), by = "Compound_ID")
    df <-
      merge(df1_full,
        match_df,
        by.x = "Compound_ID",
        by.y = "df1",
        all = TRUE
      ) |>
      merge(df2_full,
        by.x = "df2",
        by.y = "Compound_ID",
        all = TRUE
      )

    df <- df %>%
      dplyr::rename_with(
        ~ ifelse(
          .x %in% names(df),
          c(
            paste("Compound_ID", study1_name, sep = "_"),
            paste("Compound_ID", study2_name, sep = "_"),
            paste("Metabolite", study1_name, sep = "_"),
            paste("Metabolite", study2_name, sep = "_"),
            paste("RT", study1_name, sep = "_"),
            paste("RT", study2_name, sep = "_"),
            paste("MZ", study1_name, sep = "_"),
            paste("MZ", study2_name, sep = "_"),
            paste("Intensity", study1_name, sep = "_"),
            paste("Intensity", study2_name, sep = "_")
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
      ) |>
      dplyr::mutate(
        rep_Compound_ID = dplyr::case_when(
          !is.na(get(Compound_ID_1)) ~
            get(Compound_ID_1),
          is.na(get(Compound_ID_1)) & !is.na(get(Compound_ID_2)) ~ get(Compound_ID_2),
          TRUE ~ NA
        ),
        rep_RT = dplyr::case_when(
          !is.na(get(RT_1)) ~ get(RT_1),
          is.na(get(RT_1)) & !is.na(get(RT_2)) ~ get(RT_2),
          TRUE ~ NA
        ),
        rep_MZ = dplyr::case_when(
          !is.na(get(MZ_1)) ~ get(MZ_1),
          is.na(get(MZ_1)) & !is.na(get(MZ_2)) ~ get(MZ_2),
          TRUE ~ NA
        ),
        rep_Intensity = dplyr::case_when(
          !is.na(get(Intensity_1)) ~ get(Intensity_1),
          is.na(get(Intensity_1)) & !is.na(get(Intensity_2)) ~ get(Intensity_2),
          TRUE ~ NA
        ),
        rep_Metabolite = dplyr::case_when(
          !is.na(get(Metabolite_1)) ~ get(Metabolite_1),
          is.na(get(Metabolite_1)) & !is.na(get(Metabolite_2)) ~ get(Metabolite_2),
          TRUE ~ NA
        )
      ) |>
      dplyr::select(
        c(
          "rep_Compound_ID", "rep_RT", "rep_MZ", "rep_Intensity",
          "rep_Metabolite", Compound_ID_1,
          Compound_ID_2, Metabolite_1, Metabolite_2, RT_1, RT_2,
          MZ_1, MZ_2, Intensity_1, Intensity_2, dplyr::everything(),
          -dplyr::contains("_adj")
        )
      )

    df <- df %>%
      dplyr::group_by(rep_Compound_ID) |>
      dplyr::mutate(dup_count = dplyr::row_number()) |>
      dplyr::ungroup() |>
      dplyr::mutate(rep_Compound_ID = ifelse(dup_count > 1,
        paste0(rep_Compound_ID, "_", dup_count),
        rep_Compound_ID
      )) |>
      dplyr::select(-dup_count)

    all_matched(align_ms_obj) <- df
    adjusted_df(align_ms_obj) <- df2_adj

    return(align_ms_obj)
  }
