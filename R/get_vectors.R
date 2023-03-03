get_vectors <- function(df, rt_sim, mz_sim) {
  rt_metabolites <- c()
  mz_metabolites <- c()
  df_rt <- df |>
    dplyr::arrange(RT)
  df_mz <- df |>
    dplyr::arrange(MZ)
  for (row in 1:(nrow(df_rt) - 1)) {
    rt <- df_rt[row, "RT"]
    if (row == 1) {
      diff_rt <- df_rt[row + 1, "RT"] - rt
      if (diff_rt > rt_sim) {
        rt_metabolites <- c(rt_metabolites, df_rt[row, "Compound_ID"])
      }
    } else {
      diff_down_rt <- df_rt[row + 1, "RT"] - rt
      diff_up_rt <- rt - df_rt[row - 1, "RT"]
      if (diff_down_rt >= rt_sim &
          diff_up_rt >= rt_sim) {
        rt_metabolites <- c(rt_metabolites, df_rt[row, "Compound_ID"])
      }
    }
  }
  for (row in 1:(nrow(df_mz) - 1)) {
    mz <- df_mz[row, "MZ"]
    if (row == 1) {
      diff_mz <- df_mz[row + 1, "MZ"] - df_mz[row, "MZ"]
      if (diff_mz > mz_sim * df_mz[row, "MZ"] / 1e6) {
        mz_metabolites <- c(mz_metabolites, df_mz[row, "Compound_ID"])
      }
    } else {
      diff_down_mz <- df_mz[row + 1, "MZ"] - mz
      diff_up_mz <- mz - df_mz[row - 1, "MZ"]
      if (diff_down_mz > mz_sim * mz / 1e6 &
          diff_up_mz > mz_sim * mz / 1e6) {
        mz_metabolites <- c(mz_metabolites, df_mz[row, "Compound_ID"])
      }
    }
  }
  metabolites <- mz_metabolites
  if (is.null(metabolites)) {
    stop("Not enough features isolated. Consider lowering RT and MZ thresholds.")
  }
  return(metabolites)
}
