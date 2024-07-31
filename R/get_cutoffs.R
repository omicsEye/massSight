get_cutoffs <-
  function(df1,
           df2,
           has_int = TRUE) {
    data_rt <- df2$RT - df1$RT
    data_mz <- (df2$MZ - df1$MZ) / df1$MZ * 1e6

    if (has_int) {
      data_int <- log10(df2$Intensity) - log10(df1$Intensity)
      not_outliers <- !mad_based_outlier(data_rt) &
        !mad_based_outlier(data_mz) &
        !mad_based_outlier(data_int)
      data_int <- replace(data_int, data_int %in% c(Inf, -Inf), NA)
      cutoffs <- c(
        stats::sd(data_rt[not_outliers]),
        stats::sd(data_mz[not_outliers]),
        stats::sd(data_int[not_outliers])
      )
    } else {
      not_outliers <- !mad_based_outlier(data_rt) &
        !mad_based_outlier(data_mz)
      cutoffs <- c(
        stats::sd(data_rt[not_outliers]),
        stats::sd(data_mz[not_outliers])
      )
    }

    # TODO Fix below
    # rt_outliers <- df1[mad_based_outlier(data_rt)] |>
    #   rownames()
    # mz_outliers <- df1[mad_based_outlier(data_mz)] |>
    #   rownames()
    # int_outliers <- df1[mad_based_outlier(data_int)] |>
    #   rownames()
    # outliers <- c(rt_outliers, mz_outliers, int_outliers)
    outliers <- "tmp"
    return(list(
      "cutoffs" = cutoffs,
      "outliers" = outliers
    ))
  }
