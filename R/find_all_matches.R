find_all_matches <- function(ref, query, rt_threshold, mz_threshold) {
  matches <- data.frame()
  ref <- ref |>
    dplyr::rename(
      RT_1 = RT,
      MZ_1 = MZ,
      Compound_ID_1 = Compound_ID
    )
  query <- query |>
    dplyr::rename(
      RT_2 = RT,
      MZ_2 = MZ,
      Compound_ID_2 = Compound_ID
    )
  for (i in 1:nrow(ref)) {
    rt <- ref$RT_1[i]
    mz <- ref$MZ_1[i]
    rt_lower <- rt - rt_threshold
    rt_upper <- rt + rt_threshold
    mz_lower <- mz - mz_threshold
    mz_upper <- mz + mz_threshold
    query_matches <- query |>
      dplyr::filter(RT_2 >= rt_lower &
        RT_2 <= rt_upper &
        MZ_2 >= mz_lower &
        MZ_2 <= mz_upper)
    if (nrow(query_matches) > 0) {
      matches <- matches |>
        dplyr::bind_rows(dplyr::bind_cols(ref[i, ], query_matches))
    }
  }
  return(matches)
}
