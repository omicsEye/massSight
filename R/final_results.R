final_results <-
  function(align_ms_obj,
           keep_features = c(FALSE, FALSE),
           weights = c(1, 1, 1)) {
    df1 <- raw_df(ms1(align_ms_obj))
    df2 <- raw_df(ms2(align_ms_obj))
    scaled_df <- scaled_values(align_ms_obj)
    stds <- cutoffs(align_ms_obj)
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
              dplyr::filter(Compound_ID == best_match),
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
          c(best_hits_found, rep(df1_for_align[i, "Compound_ID"], length(best_match)))
      } else {
        features_not_aligned <-
          c(features_not_aligned, rownames(df1_for_align[i, "Compound_ID"]))
      }
    }
    df_not_found <- df1 |>
      dplyr::filter(Compound_ID %in% features_not_aligned)
    rownames(df1) <- df1$Compound_ID
    rownames(df2) <- df2$Compound_ID
    rownames(df2_adj) <- df2_adj$Compound_ID
    df1 <- df1[best_hits_found, ]
    df2_raw <- df2
    df2 <- df2[best_hits_df1, ]
    results_df_complete <- cbind(df1, df2)
    df2_adj <- df2_adj[best_hits_df1, ]
    if ("Intensity" %in% colnames(df1)) {
      if ("Metabolite" %in% colnames(df1)) {
        results_df <- data.frame(
          "df1" = df1$Compound_ID,
          "RT" = df1$RT,
          "MZ" = df1$MZ,
          "Intensity" = df1$Intensity,
          "Metabolite" = df1$Metabolite,
          "df2" = df2$Compound_ID,
          "RT_2" = df2$RT,
          "MZ_2" = df2$MZ,
          "Intensity_2" = df2$Intensity,
          "Metabolite_2" = df2$Metabolite
        )
      } else {
        results_df <- data.frame(
          "df1" = df1$Compound_ID,
          "RT" = df1$RT,
          "MZ" = df1$MZ,
          "Intensity" = df1$Intensity,
          "df2" = df2$Compound_ID,
          "RT_2" = df2$RT,
          "MZ_2" = df2$MZ,
          "Intensity_2" = df2$Intensity
        )
      }

      adjusted_df <- data.frame(
        "rt_2_adj" = df2_adj$RT,
        "mz_2_adj" = df2_adj$MZ,
        "int_2_adj" = df2_adj$Intensity
      )
    } else {
      if ("Metabolite" %in% colnames(df1)) {
        results_df <- data.frame(
          "df1" = df1$Compound_ID,
          "df2" = df2$Compound_ID,
          "RT" = df1$RT,
          "RT_2" = df2$RT,
          "MZ" = df1$MZ,
          "MZ_2" = df2$MZ,
          "Metabolite" = df1$Metabolite,
          "Metabolite_2" = df2$Metabolite
        )
      } else {
        results_df <- data.frame(
          "df1" = df1$Compound_ID,
          "df2" = df2$Compound_ID,
          "RT" = df1$RT,
          "RT_2" = df2$RT,
          "MZ" = df1$MZ,
          "MZ_2" = df2$MZ
        )
      }

      adjusted_df <- data.frame(
        "rt_2_adj" = df2_adj$RT,
        "mz_2_adj" = df2_adj$MZ
      )
    }


    if (keep_features[1] == T) {
      message("keeping file 1 features")
      columns <- colnames(results_df_complete)
      # TODO
    }

    if (keep_features[2] == T) {
      # TODO
    }

    # columns <- colnames(results_df_complete) |>
    #   dedup("RT") |>
    #   dedup("MZ") |>
    #   dedup("Intensity") |>
    #   dedup("Compound_ID") |>
    #   dedup("Metabolite")
    #
    # colnames(results_df_complete) <- columns
    all_matched(align_ms_obj) <- results_df
    adjusted_df(align_ms_obj) <- adjusted_df
    return(align_ms_obj)
  }
