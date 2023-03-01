find_closest_match <-
  function(rt_mz_int,
           database,
           stds,
           multipliers,
           weights) {
    database_index <- database$Compound_ID
    cutoffs <- multipliers * stds
    if (cutoffs[1] > 0.0) {
      rt_hits <- database$RT <= (rt_mz_int$RT + cutoffs[1]) &
        database$RT >= (rt_mz_int$RT - cutoffs[1])
    } else {
      rt_hits <- database$RT <= (rt_mz_int$RT + .05) &
        database$RT >= (rt_mz_int$RT - .05)
    }

    if (cutoffs[2] > 0) {
      mass_plus <- rt_mz_int$MZ * cutoffs[2] / 10000
      mass_minus <-
        rt_mz_int$MZ * cutoffs[2] / (10000 + cutoffs[2])
      mz_hits <- database$MZ < (rt_mz_int$MZ + mass_plus) &
        database$MZ > (rt_mz_int$MZ - mass_minus)
    } else {
      mz_hits <- database$MZ < (rt_mz_int$MZ + 0.005) &
        database$MZ > (rt_mz_int$MZ - 0.005)
    }

    if (cutoffs[3] > 0) {
      cutoffs_2 <- 10.00000**cutoffs[3]
      int_hits <-
        database$Intensity < (rt_mz_int$Intensity * cutoffs_2) &
          database$Intensity > (rt_mz_int$Intensity / cutoffs_2)
    } else {
      int_hits <- rep(TRUE, nrow(database))
    }
    hits <-
      database[rt_hits &
        mz_hits & int_hits, c("RT", "MZ", "Intensity")]
    hits_index <- database_index[rt_hits & mz_hits & int_hits]
    if (nrow(hits) == 0) {
      return(NULL)
    }
    hits_results <- c()
    for (i in 1:nrow(hits)) {
      score <- rms(
        rt_mz_int,
        hits[i, ],
        weights,
        stds
      )
      hits_results <- c(hits_results, score)
    }
    return(hits_index[hits_results == min(hits_results)])
  }
