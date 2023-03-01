dedup <- function(cols, item) {
  item_locations <- which(cols == item)
  for (i in 1:length(item_locations)) {
    if (i != 1) {
      cols[item_locations[i]] <- paste0(item, "_", i)
    }
  }
  return(cols)
}

mad_based_outlier <- function(points, thresh = 5) {
  points <- points[[1]]
  diff <- sqrt((points - median(points))^2)
  med_abs_deviation <- median(diff)
  if (med_abs_deviation == 0) {
    mod_z_score <- rep(0, length(diff))
  } else {
    mod_z_score <- .6745 * diff / med_abs_deviation
  }
  return(mod_z_score > thresh)
}

scale_intensity <- function(data, intensity) {
  return(data / intensity)
}

scale_intensity_parameters <-
  function(data_int1,
           data_int2,
           int_col = "Intensity",
           min_int = 0) {
    min_int <- log10(min_int)
    int_key <- (data_int1 > min_int) & (data_int2 > min_int)
    fit <- mean(data_int2[int_key] / data_int1[int_key])
    return(fit)
  }

calculate_mean <- function(values, w, method = "arit") {
  if (sum(values) == 0) {
    return(0)
  }
  if (method == "whm") {
    num <- replace(w, values < 0, 0)
    denom <- replace(w, values < 0, 0)
    denom <- replace(denom, values > 0, denom / values)
    w_h_mean <- sum(num) / sum(denom)
    nonzeroes <- sum(values != 0)
    mean_value <- w_h_mean / (nonzeroes + 1)
  } else if (method == "arit") {
    mean_value <- sum(w * values) / sum(w)
  }
  return(mean_value)
}

rms <- function(a, b, weights, std) {
  rt_difference <- a$RT - b$RT
  if (std[1] == 0) {
    rt_contribution <- abs(rt_difference)
  } else {
    rt_contribution <- (rt_difference / std[1])**2
  }

  # convert to mzs to ppm to score
  ppm_difference <- (a$MZ - b$MZ) * 1e6 / mean(c(a$MZ, b$MZ))
  if (std[2] == 0) {
    mz_contribution <- abs(ppm_difference)
  } else {
    mz_contribution <- (ppm_difference / std[2])^2
  }

  int_difference <- log10(a$Intensity) - log10(b$Intensity)
  if (std[3] == 0) {
    int_contribution <- abs(int_difference)
  } else {
    int_contribution <- (int_difference / std[3])^2
  }
  score <-
    calculate_mean(
      c(rt_contribution, mz_contribution, int_contribution),
      weights
    )
  return(score)
}

scale_smooth <- function(query_values, smooth_x, smooth_y) {

  suppressWarnings(f <- approx(
    x = smooth_x,
    y = smooth_y,
    xout = query_values,
    rule = 2
  ))

  query_out <- query_values - f$y
  # smooth_min <- smooth_y[which.min(smooth_x)]
  # smooth_max <- smooth_y[which.max(smooth_x)]
  #
  # query_values <- na.omit(query_values)
  #
  # query_values[query_values <= min(smooth_x)] <-
  #   query_values[query_values <= min(smooth_x)] - smooth_min
  # query_values[query_values >= max(smooth_x)] <-
  #   query_values[query_values >= max(smooth_x)] - smooth_max
  #
  # smooth_df <- data.frame("x" = smooth_x,
  #                         "y" = smooth_y) |>
  #   dplyr::arrange(x)
  #
  # query_values[query_values > min(smooth_x) &
  #                query_values < max(smooth_x)] <-
  #   pracma::interp1(smooth_df$x,
  #                   smooth_df$y,
  #                   query_values[query_values > min(smooth_x) &
  #                                  query_values < max(smooth_x)])

  return(query_out)
}
