align_isolated_compounds <-
  function(df1,
           df2,
           rt_minus = -.5,
           rt_plus = .5,
           mz_minus = -15,
           mz_plus = 15,
           keep = FALSE) {
    if ("Intensity" %in% names(df1)) {
      pb <- progress::progress_bar$new(format = "Matching isolated features from datasets [:bar] :percent :eta",
                                       total = nrow(df1),
                                       clear = F)
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
          dplyr::filter(# RT > df1[row, "RT"] + rt_minus &
            #   RT < df1[row, "RT"] + rt_plus &
            MZ > df1[row, "MZ"] + mz_minus * df1[row, "MZ"] / 1e6 &
              MZ < df1[row, "MZ"] + mz_plus * df1[row, "MZ"] / 1e6)
        if (nrow(df2_filter) == 0) {
          next
        }
        for (row_2 in 1:nrow(df2_filter)) {
          res_add <- data.frame(
            "df1" = df1[row, "Compound_ID"],
            "RT" = df1[row, "RT"],
            "MZ" = df1[row, "MZ"],
            "Intensity" = df1[row, "Intensity"],
            "df2" = df2_filter[row_2, "Compound_ID"],
            "RT_2" = df2_filter[row_2, "RT"],
            "MZ_2" = df2_filter[row_2, "MZ"],
            "Intensity_2" = df2_filter[row_2, "Intensity"]
          )
          results <- results |>
            rbind(res_add)
        }
      }
      return(results)
    } else {
      results <- data.frame(
        "df1" = character(),
        "RT" = numeric(),
        "MZ" = numeric(),
        "df2" = character(),
        "RT_2" = numeric(),
        "MZ_2" = numeric(),
        "del RT" = numeric(),
        "del ppm" = numeric(),
        "multiple matches" = character(),
        "multiple queries" = character()
      )

      for (row in 1:nrow(df1)) {
        df2_filter <- df2 |>
          dplyr::filter(
            RT > (df1[row, "RT"] + rt_minus),
            RT < (df1[row, "RT"] + rt_plus),
            MZ > (df1[row, "MZ"] + mz_minus / 1e6),
            MZ < (df1[row, "MZ"] + mz_plus / 1e6)
          )

        if (nrow(df2_filter) == 1) {
          del_RT <- df2_filter$RT - x$RT
          del_ppm <- (df2_filter$MZ - x$MZ) / x$MZ * 1000000
          results <- results |>
            rbind(
              c(
                x$Compound_ID,
                x$RT,
                x$MZ,
                df2_filter$Compound_ID,
                df2_filter$RT,
                df2_filter$MZ,
                del_RT,
                del_ppm,
                NA,
                NA
              )
            )
        } else if (nrow(df2_filter > 1)) {
          for (row_2 in 1:nrow(df2_filter)) {
            del_RT <- df2_filter[row_2, "RT"] - df1[row, "RT"]
            del_ppm <-
              (df2_filter[row_2, "MZ"] - df1[row, "MZ"]) / df1[row, "MZ"] * 1e6
            results <- results |>
              dplyr::bind_rows(
                c(
                  df1[row, "Compound_ID"],
                  df1[row, "RT"],
                  df1[row, "MZ"],
                  df2_filter[row_2, "Compound_ID"],
                  df2_filter[row_2, "RT"],
                  df2_filter[row_2, "MZ"],
                  del_RT,
                  del_ppm,
                  "Multiple",
                  NA
                )
              )
          }
        } else {
          results <- results |>
            dplyr::bind_rows(c(df1[row, "Compound_ID"],
                               df1[row, "RT"],
                               df1[row, "MZ"],
                               NA,
                               NA,
                               NA,
                               NA,
                               NA,
                               NA,
                               NA))
        }
      }
    }
    return(results)
  }
