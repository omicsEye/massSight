find_closest_match <-
  function(query,
           ref,
           stds) {
    ref_index <- ref$Compound_ID
    cutoffs <- stds

    rt_hits <- (ref$RT <= (query$RT + .5)) &
      (ref$RT >= (query$RT - .25))

    mz_hits <- (ref$MZ < (query$MZ + .5)) &
      (ref$MZ > (query$MZ - .5))

    if (!(TRUE %in% (rt_hits & mz_hits))) {
      return(NULL)
    }

    hits <- ref |>
      dplyr::select(dplyr::any_of(c("Compound_ID", "RT", "MZ", "Intensity"))) |>
      dplyr::filter(rt_hits & mz_hits)
    hits_index <- ref_index[rt_hits & mz_hits]
    hits_results <- c()
    for (i in 1:nrow(hits)) {
      score <- rms(
        query,
        hits[i, ],
        stds
      )
      hits_results <- c(hits_results, score)
    }
    return(hits_index[hits_results == min(hits_results)])
  }
