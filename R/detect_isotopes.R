#' Detect and label isotopes in an MSObject
#' @param ms_obj An MSObject containing mass spectrometry data
#' @param ppm_tolerance Mass tolerance in PPM for matching isotopes
#' @param rt_tolerance RT window to look for isotopes
#' @return An MSObject with isotopic labels added
#' 
#' @export 
detect_isotopes <- function(ms_obj, ppm_tolerance = 5, rt_tolerance = 0.1) {
  # Extract raw data from the MSObject
  raw_data <- raw_df(ms_obj)
  
  # Known mass differences for common isotopes
  isotope_differences <- list(
    C13 = 1.003355,   # 13C vs 12C
    N15 = 0.997035,   # 15N vs 14N
    O18 = 2.004244,   # 18O vs 16O
    S34 = 1.995796    # 34S vs 32S
  )
  
  # Function to find isotopic peaks for a given compound
  find_isotopes <- function(compound_row) {
    base_mz <- compound_row$MZ
    compound_id <- compound_row$Compound_ID
    rt <- compound_row$RT
    
    isotopes <- list()
    
    for (isotope in names(isotope_differences)) {
      mass_diff <- isotope_differences[[isotope]]
      expected_mz <- base_mz + mass_diff
      ppm_window <- (ppm_tolerance * expected_mz) / 1e6
      
      nearby_peaks <- raw_data %>%
        dplyr::filter(
          abs(RT - rt) <= rt_tolerance,
          abs(MZ - expected_mz) <= ppm_window
        )
      
      if (nrow(nearby_peaks) > 0) {
        isotopes[[isotope]] <- nearby_peaks$Compound_ID[1]
      }
    }
    
    return(isotopes)
  }
  
  # Apply isotope detection to each compound
  isotope_labels <- raw_data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Isotopes = list(find_isotopes(cur_data())))
  
  # Add isotope labels to the MSObject
  raw_df(ms_obj) <- isotope_labels
  
  return(ms_obj)
}


#' Group isotopes with their parent compounds and prepare for matching
#' @param ms_obj MSObject with isotope labels
#' @return MSObject with consolidated isotope groups
#' @export
consolidate_isotopes <- function(ms_obj) {
    df <- raw_df(ms_obj)
    
    # Create groups of related compounds (parent + isotopes)
    isotope_groups <- df %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            group_id = Compound_ID,  # Use parent compound ID as group ID
            is_parent = length(Isotopes) > 0,  # True if compound has isotopes
            isotope_count = length(unlist(Isotopes))
        ) %>%
        dplyr::ungroup()
    
    # Update the MSObject
    raw_df(ms_obj) <- isotope_groups
    consolidated(ms_obj) <- TRUE
    
    return(ms_obj)
}