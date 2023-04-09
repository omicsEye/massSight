final_results <-
  function(align_ms_obj,
           keep_features = c(FALSE, FALSE),
           multipliers = c(6, 6, 6),
           weights = c(1, 1, 1)) {
    df1 <- ms1(align_ms_obj)
    df2 <- ms2(align_ms_obj)
    scaled_df <- scaled_values(align_ms_obj)
    stds <- cutoffs(align_ms_obj)
    df2_adj <- df2
    df2_adj$RT <- scaled_df$RT
    df2_adj$MZ <- scaled_df$MZ
    df2_adj$Intensity <- scaled_df$Intensity
    df1_for_align <-
      df1[, c("Compound_ID", "RT", "MZ", "Intensity")]
    df2_for_align <-
      df2_adj[, c("Compound_ID", "RT", "MZ", "Intensity")]

    best_hits_df1 <- c()
    best_hits_found <- c()
    features_not_aligned <- c()
    pb <-
      progress::progress_bar$new(format = "Aligning datasets [:bar] :percent :eta",
                                 total = nrow(df1_for_align),
                                 clear = F)
    for (i in 1:nrow(df1_for_align)) {
      best_match <-
        find_closest_match(df1_for_align[i,],
                           df2_for_align,
                           stds,
                           multipliers,
                           weights)
      if (!is.null(best_match)) {
        pb$tick()
        best_reverse_match <-
          find_closest_match(
            df2_for_align |>
              dplyr::filter(Compound_ID == best_match),
            df1_for_align,
            stds,
            multipliers,
            weights
          )
      } else {
        features_not_aligned <-
          c(features_not_aligned, df1_for_align[i, "Compound_ID"])
        pb$tick()
        next
      }

      if (as.vector(df1_for_align[i, "Compound_ID"]) == best_reverse_match) {
        best_hits_df1 <- c(best_hits_df1, best_match)
        best_hits_found <-
          c(best_hits_found, df1_for_align[i, "Compound_ID"])
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
    df1 <- df1[best_hits_found,]
    df2_raw <- df2
    df2 <- df2[best_hits_df1,]
    results_df_complete <- cbind(df1, df2)
    df2_adj <- df2_adj[best_hits_df1,]
    results_df <- data.frame(
      "df1_name" = df1$Compound_ID,
      "df2_name" = df2$Compound_ID,
      "df1_rt" = df1$RT,
      "df2_rt" = df2$RT,
      "df1_mz" = df1$MZ,
      "df2_mz" = df2$MZ,
      "df1_int" = df1$Intensity,
      "df2_int" = df2$Intensity
    )

    adjusted_df <- data.frame(
      "rt_2_adj" = df2_adj$RT,
      "mz_2_adj" = df2_adj$MZ,
      "int_2_adj" = df2_adj$Intensity
    )

    if (keep_features[1] == T) {
      message("keeping file 1 features")
      columns <- colnames(results_df_complete)
      # TODO
    }

    if (keep_features[2] == T) {
      # TODO
    }

    columns <- colnames(results_df_complete) |>
      dedup("RT") |>
      dedup("MZ") |>
      dedup("Intensity") |>
      dedup("Compound_ID") |>
      dedup("Metabolite")

    colnames(results_df_complete) <- columns
    all_matched(align_ms_obj) <- results_df_complete
    adjusted_df(align_ms_obj) <- adjusted_df
    return(align_ms_obj)
  }
