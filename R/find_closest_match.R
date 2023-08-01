find_closest_match <-
  function(query,
           ref,
           stds) {
    ref_index <- ref$Compound_ID
    cutoffs <- stds
    if (cutoffs[1] > 0.0) {
      rt_hits <- ref$RT <= (query$RT + cutoffs[1]) &
        ref$RT >= (query$RT - cutoffs[1])
    } else {
      rt_hits <- ref$RT <= (query$RT + .05) &
        ref$RT >= (query$RT - .05)
    }

    if (cutoffs[2] > 0) {
      mass_plus <- query$MZ * cutoffs[2] / 10000
      mass_minus <-
        query$MZ * cutoffs[2] / (10000 + cutoffs[2])
      mz_hits <- ref$MZ < (query$MZ + mass_plus) &
        ref$MZ > (query$MZ - mass_minus)
    } else {
      mz_hits <- ref$MZ < (query$MZ + 0.005) &
        ref$MZ > (query$MZ - 0.005)
    }

    if (length(cutoffs) > 2 & cutoffs[3] > 0) {
      cutoffs_2 <- 10.00000**cutoffs[3]
      int_hits <-
        ref$Intensity < (query$Intensity * cutoffs_2) &
          ref$Intensity > (query$Intensity / cutoffs_2)
    } else {
      int_hits <- rep(TRUE, nrow(ref))
    }
    if (!(TRUE %in% (rt_hits & mz_hits & int_hits))) {
      return(NULL)
    }
    hits <- ref |>
      dplyr::select(any_of(c("Compound_ID", "RT", "MZ", "Intensity"))) |>
      dplyr::filter(rt_hits & mz_hits & int_hits)
    hits_index <- ref_index[rt_hits & mz_hits & int_hits]
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
