#' Detect Adducts in Mass Spectrometry Data
#'
#' This function identifies potential adduct pairs in mass spectrometry data based on mass differences and retention times.
#'
#' @param data A data frame containing at least the columns \code{RT} (retention time in minutes) and \code{MZ} (mass-to-charge ratio).
#' @param ppm_tolerance Mass accuracy tolerance in parts per million (ppm) for matching mass differences. Default is 5 ppm.
#' @param rt_tolerance Retention time tolerance in minutes for considering peaks as potential adducts. Default is 0.1 minutes.
#' @param ion_mode Ionization mode, either \code{"pos"} for positive or \code{"neg"} for negative. Default is \code{"pos"}.
#' @param adducts Optional data frame of adduct definitions. If not provided, a default list of common adducts is used based on \code{ion_mode}.
#' @return A data frame of potential adduct pairs with details of the matches.
#' @examples
#' # Sample data
#' data <- data.frame(
#'   RT = c(1.0, 1.0, 1.0, 1.0),
#'   MZ = c(100.0000, 101.0073, 122.9892, 138.9632)
#' )
#' # Detect adducts in positive mode
#' results <- detect_adducts(data, ion_mode = "pos")
#' print(results)
#' @export
detect_adducts <- function(data, ppm_tolerance = 5, rt_tolerance = 0.1, ion_mode = "pos", adducts = NULL) {
  # Ensure that required columns are present
  if (!all(c("RT", "MZ") %in% colnames(data))) {
    stop("Data frame must contain 'RT' and 'MZ' columns.")
  }

  # Define default adduct lists
  default_adducts_pos <- data.frame(
    adduct = c("[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+"),
    massdiff = c(1.007276, 22.989218, 38.963158, 18.033823)
  )

  default_adducts_neg <- data.frame(
    adduct = c("[M-H]-", "[M+Cl]-", "[M+FA-H]-", "[M+CH3COO]-"),
    massdiff = c(-1.007276, 34.969402, 44.998201, 59.013851)
  )

  # Select adduct list based on ion_mode
  if (is.null(adducts)) {
    if (ion_mode == "pos") {
      adducts <- default_adducts_pos
    } else if (ion_mode == "neg") {
      adducts <- default_adducts_neg
    } else {
      stop("Invalid ion_mode. Use 'pos' for positive or 'neg' for negative ionization mode.")
    }
  }

  # Convert data to matrix for faster computations
  mz_values <- data$MZ
  rt_values <- data$RT

  # Create all combinations of peaks
  comb_indices <- which(upper.tri(matrix(1, nrow = length(mz_values), ncol = length(mz_values))), arr.ind = TRUE)
  idx1 <- comb_indices[, 1]
  idx2 <- comb_indices[, 2]

  mz1 <- mz_values[idx1]
  mz2 <- mz_values[idx2]
  rt1 <- rt_values[idx1]
  rt2 <- rt_values[idx2]

  # Calculate differences
  mz_diff <- mz2 - mz1
  rt_diff <- abs(rt2 - rt1)

  # Filter based on retention time tolerance
  rt_mask <- rt_diff <= rt_tolerance

  if (!any(rt_mask)) {
    message("No peak pairs within the specified retention time tolerance.")
    return(NULL)
  }

  # Subset the data based on rt_mask
  mz1 <- mz1[rt_mask]
  mz2 <- mz2[rt_mask]
  rt1 <- rt1[rt_mask]
  rt2 <- rt2[rt_mask]
  mz_diff <- mz_diff[rt_mask]

  # Initialize result list
  results_list <- list()
  result_index <- 1

  # Vectorized operation over adducts
  for (k in 1:nrow(adducts)) {
    theoretical_diff <- adducts$massdiff[k]
    ppm_error <- abs((mz_diff - theoretical_diff) / theoretical_diff) * 1e6
    ppm_mask <- ppm_error <= ppm_tolerance

    if (any(ppm_mask)) {
      matches <- data.frame(
        mz1 = mz1[ppm_mask],
        rt1 = rt1[ppm_mask],
        mz2 = mz2[ppm_mask],
        rt2 = rt2[ppm_mask],
        observed_diff = mz_diff[ppm_mask],
        adduct = adducts$adduct[k],
        theoretical_diff = theoretical_diff,
        ppm_error = ppm_error[ppm_mask]
      )
      results_list[[result_index]] <- matches
      result_index <- result_index + 1
    }
  }

  # Combine results into a data frame
  if (length(results_list) > 0) {
    results <- do.call(rbind, results_list)
    # Sort results by ppm_error
    results <- results[order(results$ppm_error), ]
    return(results)
  } else {
    message("No adducts detected within the specified tolerances.")
    return(NULL)
  }
}
