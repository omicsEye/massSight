find_closest_match <-
  function(query,
           ref,
           stds) {
    ref_index <- ref$Compound_ID

    rt_hits <- (ref$RT <= (query$RT + .5)) &
      (ref$RT >= (query$RT - .5))

    mz_hits <- (ref$MZ < (query$MZ + .1)) &
      (ref$MZ > (query$MZ - .1))

    combined_hits <- rt_hits & mz_hits

    if (!(TRUE %in% (combined_hits))) {
      return(NULL)
    }

    hits <- ref |>
      dplyr::select(dplyr::any_of(c("Compound_ID", "RT", "MZ", "Intensity"))) |>
      dplyr::filter(rt_hits & mz_hits)
    hits_index <- ref_index[rt_hits & mz_hits]
    hits_results <- c()
    for (i in seq_len(nrow(hits))) {
      score <- rms(
        query,
        hits[i, ],
        stds
      )
      hits_results <- c(hits_results, score)
    }
    return(hits_index[hits_results == min(hits_results)])
  }
