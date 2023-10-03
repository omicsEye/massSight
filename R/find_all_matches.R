find_all_matches <- function(ref, query, rt_threshold, mz_threshold) {
  matches <- data.frame()
  ref <- ref |>
    dplyr::rename(
      RT_1 = .data$RT,
      MZ_1 = .data$MZ,
      Compound_ID_1 = .data$Compound_ID
    )
  query <- query |>
    dplyr::rename(
      RT_2 = .data$RT,
      MZ_2 = .data$MZ,
      Compound_ID_2 = .data$Compound_ID
    )
  for (i in 1:nrow(ref)) {
    rt <- ref$RT_1[i]
    mz <- ref$MZ_1[i]
    rt_lower <- rt - rt_threshold
    rt_upper <- rt + rt_threshold
    mz_lower <- mz - mz_threshold
    mz_upper <- mz + mz_threshold
    query_matches <- query |>
      dplyr::filter(.data$RT_2 >= rt_lower &
        .data$RT_2 <= rt_upper &
        .data$MZ_2 >= mz_lower &
        .data$MZ_2 <= mz_upper)
    if (nrow(query_matches) > 0) {
      matches <- matches |>
        dplyr::bind_rows(dplyr::bind_cols(ref[i, ], query_matches))
    }
  }
  return(matches)
}
