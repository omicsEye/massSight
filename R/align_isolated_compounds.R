align_isolated_compounds <-
  function(align_ms_obj,
           match_method,
           rt_minus = -.5,
           rt_plus = .5,
           mz_minus = -15,
           mz_plus = 15,
           keep = FALSE) {
    df1 <- isolated(ms1(align_ms_obj))
    df2 <- isolated(ms2(align_ms_obj))
    if (match_method == "unsupervised") {
      if ("Intensity" %in% names(df1)) {
        pb <-
          progress::progress_bar$new(
            format = "Matching isolated features from datasets [:bar] :percent :eta",
            total = nrow(df1),
            clear = F
          )
        results <- data.frame(
          "df1" = character(),
          "RT" = numeric(),
          "MZ" = numeric(),
          "Intensity" = numeric(),
          "df2" = character(),
          "RT_2" = numeric(),
          "MZ_2" = numeric(),
          "Intensity_2" = numeric()
        )
        for (row in 1:nrow(df1)) {
          pb$tick()
          df2_filter <- df2 |>
            dplyr::filter(
              .data$RT > df1$RT[row] + rt_minus &
                .data$RT < df1$RT[row] + rt_plus &
                .data$MZ > df1$MZ[row] + mz_minus * df1$MZ[row] / 1e6 &
                .data$MZ < df1$MZ[row] + mz_plus * df1$MZ[row] / 1e6
            )
          if (nrow(df2_filter) > 0) {
            for (row_2 in 1:nrow(df2_filter)) {
              res_add <- data.frame(
                df1 = df1$Compound_ID[row],
                RT = df1$RT[row],
                MZ = df1$MZ[row],
                Intensity = df1$Intensity[row],
                df2 = df2_filter$Compound_ID[row_2],
                RT_2 = df2_filter$RT[row_2],
                MZ_2 = df2_filter$MZ[row_2],
                Intensity_2 = df2_filter$Intensity[row_2]
              )
              results <- results |>
                rbind(res_add)
            }
          }
        }
      } else {
        pb <-
          progress::progress_bar$new(
            format = "Matching all features from datasets [:bar] :percent :eta",
            total = nrow(df1),
            clear = F
          )
        results <- data.frame(
          "df1" = character(),
          "RT" = numeric(),
          "MZ" = numeric(),
          "df2" = character(),
          "RT_2" = numeric(),
          "MZ_2" = numeric()
        )

        for (row in 1:nrow(df1)) {
          pb$tick()
          df2_filter <- df2 |>
            dplyr::filter(
              .data$RT > (df1[row, "RT"] + rt_minus),
              .data$RT < (df1[row, "RT"] + rt_plus),
              .data$MZ > (df1[row, "MZ"] + mz_minus / 1e6),
              .data$MZ < (df1[row, "MZ"] + mz_plus / 1e6)
            )

          if (nrow(df2_filter) > 0) {
            for (row_2 in 1:nrow(df2_filter)) {
              results <- results |>
                dplyr::bind_rows(
                  data.frame(
                    "df1" = df1[row, "Compound_ID"],
                    "RT" = df1[row, "RT"],
                    "MZ" = df1[row, "MZ"],
                    "df2" = df2_filter[row_2, "Compound_ID"],
                    "RT_2" = df2_filter[row_2, "RT"],
                    "MZ_2" = df2_filter[row_2, "MZ"]
                  )
                )
            }
          }
        }
      }
    } else if (match_method == "supervised") {
      stopifnot("Metabolite" %in% colnames(df1) &
        "Metabolite" %in% colnames(df2))
      vec_1 <- df1 |>
        dplyr::rename(df1 = .data$Compound_ID) |>
        dplyr::filter(.data$Metabolite != "")
      vec_2 <- df2 |>
        dplyr::rename(
          RT_2 = .data$RT,
          MZ_2 = .data$MZ,
          Intensity_2 = .data$Intensity,
          df2 = .data$Compound_ID
        ) |>
        dplyr::filter(.data$Metabolite != "")
      results <- vec_1 |>
        dplyr::inner_join(vec_2, by = c("Metabolite"))
    }
    iso_matched(align_ms_obj) <- results
    return(align_ms_obj)
  }
