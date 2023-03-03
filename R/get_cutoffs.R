get_cutoffs <-
  function(df1,
           df2,
           columns = c("RT", "MZ", "Intensity")) {
    rt_col <- columns[1]
    mz_col <- columns[2]
    int_col <- columns[3]

    data_rt <- df2[rt_col] - df1[rt_col]
    data_mz <- (df2[mz_col] - df1[mz_col]) / df1[mz_col] * 1e6
    data_int <- log10(df2[int_col]) - log10(df1[int_col])

    not_outliers <- !mad_based_outlier(data_rt) &
      !mad_based_outlier(data_mz) &
      !mad_based_outlier(data_int)

    data_int <- replace(data_int, data_int %in% c(Inf,-Inf), NA)
    cutoffs <- c(sd(data_rt[not_outliers,]),
                 sd(data_mz[not_outliers,]),
                 sd(data_int[not_outliers,]))

    # TODO Fix below
    # rt_outliers <- df1[mad_based_outlier(data_rt)] |>
    #   rownames()
    # mz_outliers <- df1[mad_based_outlier(data_mz)] |>
    #   rownames()
    # int_outliers <- df1[mad_based_outlier(data_int)] |>
    #   rownames()
    # outliers <- c(rt_outliers, mz_outliers, int_outliers)
    outliers <- "tmp"
    return(list("cutoffs" = cutoffs,
                "outliers" = outliers))
  }
