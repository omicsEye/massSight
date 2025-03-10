#' Advanced ML Match
#'
#' @description Trains and applies a semi-supervised machine learning model for matching
#' metabolites across LC-MS studies, leveraging both labeled and unlabeled data.
#'
#' @param ms1 A `massSight` object representing the first LC-MS dataset (reference)
#' @param ms2 A `massSight` object representing the second LC-MS dataset (query)
#' @param mz_thresh `numeric` Mass-to-charge threshold in ppm used to limit potential matches
#' @param rt_thresh `numeric` Retention time threshold in minutes used to limit potential matches
#' @param prob_thresh `numeric` Probability threshold for considering a match (default: 0.5)
#' @param use_unlabeled `logical` Whether to use unlabeled data for semi-supervised learning (default: TRUE)
#' @param confidence_threshold `numeric` Confidence threshold for adding predictions to training set (default: 0.8)
#' @param max_iterations `numeric` Maximum iterations for semi-supervised learning (default: 3)
#' @param neg_examples_per_pos Number of negative examples to generate per positive example (default: 3)
#' @param use_semisupervised Whether to use semi-supervised learning (default: TRUE)
#' @param use_intensity Whether to use intensity-based features (default: TRUE)
#' @param seed `numeric` Random seed for reproducibility
#'
#' @return A list containing:
#'   \itemize{
#'     \item matched_data: Data frame of predicted matches
#'     \item model: The trained model
#'     \item training_data: The training data used
#'     \item features: List containing extracted features from both datasets
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' result <- ml_match(ms1, ms2, mz_thresh = 15, rt_thresh = 0.5, seed = 123)
#' head(result$matched_data)
#' }
ml_match <- function(ms1, 
                     ms2, 
                     mz_thresh = 15, 
                     rt_thresh = 1,
                     prob_thresh = 0.5,
                     use_unlabeled = TRUE,
                     confidence_threshold = 0.8,
                     max_iterations = 3,
                     neg_examples_per_pos = 3,
                     use_semisupervised = TRUE,
                     use_intensity = TRUE,
                     seed = 72) {
  
  # Start by setting the random seed
  set.seed(seed)
  
  message("Starting advanced ML-based metabolite matching")
  
  # Validate inputs
  if (!methods::is(ms1, "MSObject") || !methods::is(ms2, "MSObject")) {
    warning("Inputs must be MSObject instances")
    return(NULL)
  }
  
  if (nrow(raw_df(ms1)) == 0 || nrow(raw_df(ms2)) == 0) {
    warning("One or both datasets are empty")
    return(NULL)
  }
  
  # Required columns
  required_cols <- c("Compound_ID", "MZ", "RT", "Metabolite")
  
  if (!all(required_cols %in% colnames(raw_df(ms1))) || 
      !all(required_cols %in% colnames(raw_df(ms2)))) {
    warning("Required columns (Compound_ID, MZ, RT, Metabolite) missing in one or both datasets")
    return(NULL)
  }
  
  # Get initial labeled data from shared metabolites
  message("Getting shared metabolites between datasets")
  labeled_data <- get_shared_metabolites(ms1, ms2)
  
  # Check if we have any shared metabolites
  if (nrow(labeled_data$known) == 0) {
    warning("No shared metabolites found between datasets. Cannot train model.")
    return(NULL)
  }
  
  # Extract rich features from both datasets
  message("Extracting advanced features from ms1")
  features_ms1 <- extract_advanced_features(ms1)
  message("Extracting advanced features from ms2")
  features_ms2 <- extract_advanced_features(ms2)
  
  # Create initial training data from labeled pairs
  message("Creating enhanced training data")
  training_data <- create_enhanced_training_data(
    labeled_data, 
    features_ms1, 
    features_ms2, 
    mz_thresh, 
    rt_thresh,
    neg_examples_per_pos,
    use_intensity
  )
  
  # Only proceed if we have enough training data
  if (nrow(training_data) < 10) {
    warning("Insufficient training data. Need at least 10 labeled metabolite pairs.")
    return(NULL)
  }
  
  # Train initial model with error handling
  message("Training initial GLMNet model with 5-fold cross-validation")
  tryCatch({
    model_results <- train_logistic_model(training_data, seed)
    model <- model_results$model
    cv_metrics <- model_results$cv_metrics
    message(sprintf(
      "Cross-validation AUC: %.3f (SE: %.3f)",
      cv_metrics$cv_auc,
      cv_metrics$cv_auc_se
    ))
  }, error = function(e) {
    warning(paste("Error training model:", e$message))
    return(NULL)
  })
  
  if (is.null(model)) {
    warning("Failed to train GLMNet model")
    return(NULL)
  }
  
  # Generate all potential pairs within thresholds
  message("Generating potential metabolite pairs")
  all_pairs <- generate_potential_pairs(ms1, ms2, mz_thresh, rt_thresh)
  
  if (nrow(all_pairs) == 0) {
    warning("No potential metabolite pairs found within thresholds")
    return(NULL)
  }
  
  # Apply semi-supervised learning if enabled
  if (use_unlabeled && use_semisupervised && max_iterations > 0) {
    message("Applying semi-supervised learning")
    tryCatch({
      model_results <- apply_semisupervised_learning(
        model, 
        training_data,
        all_pairs,
        features_ms1,
        features_ms2,
        confidence_threshold,
        max_iterations,
        seed
      )
      model <- model_results$model
      training_data <- model_results$training_data
    }, error = function(e) {
      warning(paste("Error in semi-supervised learning:", e$message))
      # Continue with original model and training data
    })
  }
  
  # Make final predictions with constraints
  message("Making final predictions with one-to-one constraints")
  predictions <- predict_with_constraints(
    model, 
    all_pairs, 
    prob_thresh, 
    rt_thresh = rt_thresh,
    mz_thresh = mz_thresh,
    use_intensity = use_intensity
  )
  
  # Generate quality metrics
  message("Computing match quality metrics")
  match_metrics <- compute_match_metrics(predictions, labeled_data$known)
  message(sprintf(
    "Found %d matches with %.1f%% of known matches correctly identified",
    nrow(predictions),
    ifelse(is.na(match_metrics$recall), 0, match_metrics$recall * 100)
  ))
  
  # Return results
  return(list(
    "matched_data" = predictions,
    "model" = model,
    "training_data" = training_data,
    "features" = list(ms1 = features_ms1, ms2 = features_ms2),
    "metrics" = list(
      cv = cv_metrics,  # Add CV metrics
      final = match_metrics  # Keep final metrics on full dataset
    )
  ))
}

#' Extract Advanced Features from MS Object
#'
#' @param ms_obj A massSight object
#' @return A data frame with compound IDs and extracted features
extract_advanced_features <- function(ms_obj) {
  # Get raw data
  ms_df <- raw_df(ms_obj)
  
  # Handle empty or single-element cases
  if (nrow(ms_df) <= 1) {
    if (nrow(ms_df) == 0) {
      return(data.frame())
    } else {
      # Single row case - set ranks to 0.5 and deciles to 5
      ms_features <- ms_df
      ms_features$RT_rank <- 0.5
      ms_features$RT_decile <- 5
      ms_features$MZ_rank <- 0.5
      ms_features$MZ_decile <- 5
      
      if ("Intensity" %in% colnames(ms_df)) {
        ms_features$log_intensity <- log10(ms_features$Intensity)
        ms_features$intensity_rank <- 0.5
        ms_features$intensity_decile <- 5
        ms_features$local_intensity_mean <- ms_features$log_intensity
        ms_features$local_intensity_sd <- 0
        ms_features$rel_intensity <- 1
      }
      
      return(ms_features)
    }
  }
  
  # Add basic RT and MZ ranks
  ms_features <- ms_df %>%
    dplyr::arrange(RT) %>%
    dplyr::mutate(
      # Normalized RT rank (0-1)
      RT_rank = (dplyr::row_number() - 1) / max(dplyr::n() - 1, 1),
      # RT percentile groups from 1-10
      RT_decile = pmin(floor(RT_rank * 10) + 1, 10)
    ) %>%
    dplyr::arrange(MZ) %>%
    dplyr::mutate(
      # Normalized MZ rank (0-1)
      MZ_rank = (dplyr::row_number() - 1) / max(dplyr::n() - 1, 1),
      # MZ percentile groups from 1-10
      MZ_decile = pmin(floor(MZ_rank * 10) + 1, 10)
    )
  
  # Add intensity features if available
  if ("Intensity" %in% colnames(ms_df)) {
    ms_features <- ms_features %>%
      dplyr::arrange(Intensity) %>%
      dplyr::mutate(
        # Log intensity with safety for zeros
        log_intensity = log10(pmax(Intensity, 1e-10)),
        # Normalized intensity rank (0-1)
        intensity_rank = (dplyr::row_number() - 1) / max(dplyr::n() - 1, 1),
        # Intensity percentile groups from 1-10
        intensity_decile = pmin(floor(intensity_rank * 10) + 1, 10)
      )
    
    # Calculate neighborhood features
    ms_features <- ms_features %>%
      dplyr::arrange(RT) %>%
      dplyr::mutate(
        # Local intensity metrics within RT neighborhood
        local_intensity_mean = calculate_local_mean(RT, log_intensity, window = 0.5),
        local_intensity_sd = calculate_local_sd(RT, log_intensity, window = 0.5),
        # Intensity relative to local mean with safety for zero mean
        rel_intensity = log_intensity / pmax(local_intensity_mean, 1e-10)
      )
  }
  
  # Return features
  return(ms_features)
}

#' Calculate local mean for a variable within a sliding window
#'
#' @param x Vector of positions
#' @param y Vector of values
#' @param window Window size
#' @return Vector of local means
calculate_local_mean <- function(x, y, window) {
  n <- length(x)
  result <- numeric(n)
  
  for (i in 1:n) {
    # Find points within window
    in_window <- abs(x - x[i]) <= window
    # Calculate mean excluding NA values
    result[i] <- mean(y[in_window], na.rm = TRUE)
  }
  
  return(result)
}

#' Calculate local standard deviation for a variable within a sliding window
#'
#' @param x Vector of positions
#' @param y Vector of values
#' @param window Window size
#' @return Vector of local standard deviations
calculate_local_sd <- function(x, y, window) {
  n <- length(x)
  result <- numeric(n)
  
  for (i in 1:n) {
    # Find points within window
    in_window <- abs(x - x[i]) <= window
    # Calculate sd excluding NA values
    result[i] <- stats::sd(y[in_window], na.rm = TRUE)
  }
  
  return(result)
}


#' Create Enhanced Training Data
#'
#' @param labeled_data List of labeled metabolite data
#' @param features_ms1 Features from first dataset
#' @param features_ms2 Features from second dataset
#' @param mz_thresh MZ threshold for negative examples
#' @param rt_thresh RT threshold for negative examples
#' @param neg_examples_per_pos Number of negative examples to generate per positive example
#' @param use_intensity Whether to use intensity-based features
#' @return Enhanced training dataset
create_enhanced_training_data <- function(labeled_data, features_ms1, features_ms2, 
                                        mz_thresh, rt_thresh, neg_examples_per_pos = 3,
                                        use_intensity = TRUE) {
  # Get known matches
  known_matches <- labeled_data$known
  
  # If no known matches, return empty data frame
  if (nrow(known_matches) == 0) {
    return(data.frame(
      Compound_ID_1 = character(0),
      Compound_ID_2 = character(0),
      MZ_1 = numeric(0),
      MZ_2 = numeric(0),
      RT_1 = numeric(0),
      RT_2 = numeric(0),
      Class = factor(character(0), levels = c("unmatched", "matched"))
    ))
  }
  
  # Extract needed columns for positive examples with explicit selection
  positive_examples <- data.frame(
    Compound_ID_1 = known_matches$Compound_ID,
    Compound_ID_2 = known_matches$Compound_ID.y,
    MZ_1 = known_matches$MZ,
    MZ_2 = known_matches$MZ.y,
    RT_1 = known_matches$RT,
    RT_2 = known_matches$RT.y,
    Metabolite_1 = known_matches$Metabolite,
    Metabolite_2 = known_matches$Metabolite,
    Class = "matched",
    stringsAsFactors = FALSE
  )
  
  # Create negative examples: for each metabolite in ms1, find multiple nearest non-matching metabolites in ms2
  negative_examples <- data.frame()
  
  # For each known metabolite in ms1
  for (i in 1:nrow(labeled_data$ms1_known)) {
    # Get the current metabolite
    ms1_metabolite <- labeled_data$ms1_known[i, ]
    
    # Get all metabolites in ms2 with different names
    ms2_candidates <- features_ms2 %>%
      dplyr::filter(
        Metabolite != ms1_metabolite$Metabolite,
        # Within threshold for consideration
        abs(MZ - ms1_metabolite$MZ) <= mz_thresh,
        abs(RT - ms1_metabolite$RT) <= rt_thresh
      )
    
    # If no candidates, skip
    if (nrow(ms2_candidates) == 0) next
    
    # Calculate distance to each candidate
    ms2_candidates$distance <- sqrt(
      ((ms2_candidates$MZ - ms1_metabolite$MZ) / ms1_metabolite$MZ * 1e6)^2 + 
      ((ms2_candidates$RT - ms1_metabolite$RT) / rt_thresh)^2
    )
    
    # Get the closest N candidates
    closest_candidates <- ms2_candidates %>%
      dplyr::arrange(distance) %>%
      head(neg_examples_per_pos)
    
    # Create negative examples
    for (j in 1:nrow(closest_candidates)) {
      candidate <- closest_candidates[j, ]
      negative_example <- data.frame(
        Compound_ID_1 = ms1_metabolite$Compound_ID,
        Compound_ID_2 = candidate$Compound_ID,
        MZ_1 = ms1_metabolite$MZ,
        MZ_2 = candidate$MZ,
        RT_1 = ms1_metabolite$RT,
        RT_2 = candidate$RT,
        Metabolite_1 = ms1_metabolite$Metabolite,
        Metabolite_2 = candidate$Metabolite,
        Class = "unmatched"
      )
      negative_examples <- rbind(negative_examples, negative_example)
    }
  }
  
  # Join negative examples with enhanced features
  if (nrow(negative_examples) > 0) {
    training_data <- rbind(positive_examples, negative_examples)
  } else {
    training_data <- positive_examples
  }
  
  # Add derived features
  training_data <- training_data %>%
    dplyr::mutate(
      # Basic differences
      delta_RT = RT_1 - RT_2,
      delta_MZ = (MZ_1 - MZ_2) / MZ_1 * 1e6,  # in ppm
      
      # Normalized differences
      rel_delta_RT = delta_RT / rt_thresh,
      rel_delta_MZ = delta_MZ / mz_thresh,
      
      # Combined distance metrics
      euclidean_dist = sqrt(rel_delta_RT^2 + rel_delta_MZ^2),
      manhattan_dist = abs(rel_delta_RT) + abs(rel_delta_MZ),
      
      # Squared differences (for non-linear relationships)
      delta_RT_squared = delta_RT^2,
      delta_MZ_squared = delta_MZ^2,
      
      # Log-transformed MZ ratios
      log_MZ_ratio = log(MZ_1 / MZ_2),
      
      # RT position features
      RT_sum = RT_1 + RT_2,
      RT_product = RT_1 * RT_2,
      RT_min = pmin(RT_1, RT_2),
      RT_max = pmax(RT_1, RT_2),
      RT_range = RT_max - RT_min,
      
      # MZ position features
      MZ_sum = MZ_1 + MZ_2,
      MZ_product = MZ_1 * MZ_2,
      MZ_min = pmin(MZ_1, MZ_2),
      MZ_max = pmax(MZ_1, MZ_2),
      MZ_range = MZ_max - MZ_min
    )
  
  # Add rank-based features if they exist in the features data
  if ("RT_rank" %in% colnames(features_ms1) && "RT_rank" %in% colnames(features_ms2)) {
    training_data <- training_data %>%
      dplyr::left_join(
        features_ms1 %>% 
          dplyr::select(Compound_ID, RT_rank, RT_decile) %>%
          dplyr::rename(
            Compound_ID_1 = Compound_ID,
            RT_rank_1 = RT_rank,
            RT_decile_1 = RT_decile
          ),
        by = "Compound_ID_1"
      ) %>%
      dplyr::left_join(
        features_ms2 %>% 
          dplyr::select(Compound_ID, RT_rank, RT_decile) %>%
          dplyr::rename(
            Compound_ID_2 = Compound_ID,
            RT_rank_2 = RT_rank,
            RT_decile_2 = RT_decile
          ),
        by = "Compound_ID_2"
      ) %>%
      dplyr::mutate(
        delta_RT_rank = abs(RT_rank_1 - RT_rank_2),
        same_RT_decileTRUE = as.numeric(RT_decile_1 == RT_decile_2)
      )
  }
  
  # Add MZ rank features
  if ("MZ_rank" %in% colnames(features_ms1) && "MZ_rank" %in% colnames(features_ms2)) {
    training_data <- training_data %>%
      dplyr::left_join(
        features_ms1 %>% 
          dplyr::select(Compound_ID, MZ_rank, MZ_decile) %>%
          dplyr::rename(
            Compound_ID_1 = Compound_ID,
            MZ_rank_1 = MZ_rank,
            MZ_decile_1 = MZ_decile
          ),
        by = "Compound_ID_1"
      ) %>%
      dplyr::left_join(
        features_ms2 %>% 
          dplyr::select(Compound_ID, MZ_rank, MZ_decile) %>%
          dplyr::rename(
            Compound_ID_2 = Compound_ID,
            MZ_rank_2 = MZ_rank,
            MZ_decile_2 = MZ_decile
          ),
        by = "Compound_ID_2"
      ) %>%
      dplyr::mutate(
        delta_MZ_rank = abs(MZ_rank_1 - MZ_rank_2),
        same_MZ_decileTRUE = as.numeric(MZ_decile_1 == MZ_decile_2)
      )
  }
  
  # Add intensity-based features if available and enabled
  if (use_intensity && 
      all(c("log_intensity", "intensity_rank") %in% colnames(features_ms1)) &&
      all(c("log_intensity", "intensity_rank") %in% colnames(features_ms2))) {
    training_data <- training_data %>%
      dplyr::left_join(
        features_ms1 %>% 
          dplyr::select(Compound_ID, log_intensity, intensity_rank, intensity_decile,
                       local_intensity_mean, local_intensity_sd, rel_intensity) %>%
          dplyr::rename_with(~ paste0(.x, "_1"), -Compound_ID) %>%
          dplyr::rename(Compound_ID_1 = Compound_ID),
        by = "Compound_ID_1"
      ) %>%
      dplyr::left_join(
        features_ms2 %>% 
          dplyr::select(Compound_ID, log_intensity, intensity_rank, intensity_decile,
                       local_intensity_mean, local_intensity_sd, rel_intensity) %>%
          dplyr::rename_with(~ paste0(.x, "_2"), -Compound_ID) %>%
          dplyr::rename(Compound_ID_2 = Compound_ID),
        by = "Compound_ID_2"
      ) %>%
      dplyr::mutate(
        # Basic intensity differences
        delta_log_intensity = log_intensity_1 - log_intensity_2,
        delta_intensity_rank = abs(intensity_rank_1 - intensity_rank_2),
        same_intensity_decileTRUE = as.numeric(intensity_decile_1 == intensity_decile_2),
        
        # Advanced intensity features
        intensity_ratio = 10^(log_intensity_1 - log_intensity_2),
        log_intensity_product = log_intensity_1 * log_intensity_2,
        intensity_rank_product = intensity_rank_1 * intensity_rank_2,
        
        # Local intensity context features
        delta_local_intensity_mean = local_intensity_mean_1 - local_intensity_mean_2,
        delta_local_intensity_sd = local_intensity_sd_1 - local_intensity_sd_2,
        delta_rel_intensity = rel_intensity_1 - rel_intensity_2,
        
        # Combined intensity and position features
        intensity_RT_correlation = delta_log_intensity * delta_RT,
        intensity_MZ_correlation = delta_log_intensity * delta_MZ
      )
  }
  
  # Convert class to factor
  training_data$Class <- factor(training_data$Class, levels = c("unmatched", "matched"))
  
  return(training_data)
}

#' Get Shared Metabolites Between Datasets
#'
#' @param ms1 First MSObject
#' @param ms2 Second MSObject
#' @return List containing matched metabolites and related data
get_shared_metabolites <- function(ms1, ms2) {
  ms1_df <- raw_df(ms1)
  ms2_df <- raw_df(ms2)
  
  # Check if Metabolite column exists
  if (!"Metabolite" %in% colnames(ms1_df) || !"Metabolite" %in% colnames(ms2_df)) {
    warning("Metabolite column missing in one or both datasets")
    return(list(
      "known" = data.frame(),
      "ms1_known" = data.frame(),
      "ms2_known" = data.frame()
    ))
  }
  
  # Extract known metabolites from both datasets
  ms1_known <- ms1_df %>%
    dplyr::filter(
      !is.na(Metabolite),
      Metabolite != ""
    )
  
  if (nrow(ms1_known) == 0) {
    warning("No known metabolites in first dataset")
    return(list(
      "known" = data.frame(),
      "ms1_known" = data.frame(),
      "ms2_known" = data.frame()
    ))
  }
  
  ms2_known <- ms2_df %>%
    dplyr::filter(
      !is.na(Metabolite),
      Metabolite != ""
    )
  
  if (nrow(ms2_known) == 0) {
    warning("No known metabolites in second dataset")
    return(list(
      "known" = data.frame(),
      "ms1_known" = ms1_known,
      "ms2_known" = data.frame()
    ))
  }
  
  # Find shared metabolites
  shared_metabolites <- intersect(
    ms1_known$Metabolite,
    ms2_known$Metabolite
  )
  
  if (length(shared_metabolites) == 0) {
    warning("No shared metabolites found between datasets")
    return(list(
      "known" = data.frame(),
      "ms1_known" = ms1_known,
      "ms2_known" = ms2_known
    ))
  }
  
  # Create pairs with shared metabolites
  known_pairs <- data.frame(
    Metabolite = shared_metabolites
  )
  
  # Add data from ms1
  known_pairs <- known_pairs %>%
    dplyr::left_join(
      ms1_known %>% 
        dplyr::select(Metabolite, Compound_ID, MZ, RT),
      by = "Metabolite"
    )
  
  # Add data from ms2
  known_pairs <- known_pairs %>%
    dplyr::left_join(
      ms2_known %>% 
        dplyr::select(Metabolite, Compound_ID, MZ, RT),
      by = "Metabolite",
      suffix = c("", ".y")
    )
  
  return(list(
    "known" = known_pairs,
    "ms1_known" = ms1_known,
    "ms2_known" = ms2_known
  ))
}

#' Train Logistic Model
#'
#' @param training_data Training dataset
#' @param seed Random seed
#' @return Trained model and CV metrics
train_logistic_model <- function(training_data, seed = 42) {
  set.seed(seed)
  
  # Remove non-feature columns
  feature_data <- training_data %>%
    dplyr::select(-Compound_ID_1, -Compound_ID_2, -Metabolite_1, -Metabolite_2)
  
  # Convert response to numeric (0/1)
  y <- as.numeric(feature_data$Class) - 1
  x <- as.matrix(feature_data %>% dplyr::select(-Class))
  
  # Debug: Look at training data
  message("Training data summary:")
  print(summary(x))
  
  # Calculate and store scaling parameters
  scale_params <- list(
    center = apply(x, 2, mean),
    scale = apply(x, 2, sd)
  )
  
  # Debug: Look at scaling parameters
  message("Scaling parameters:")
  print(data.frame(
    feature = names(scale_params$center),
    mean = scale_params$center,
    sd = scale_params$scale
  ))
  
  # Standardize features
  x <- scale(x, center = scale_params$center, scale = scale_params$scale)
  
  # Debug: Look at scaled training data
  message("Scaled training data summary:")
  print(summary(x))
  
  # Fit model with cross-validation
  cv_fit <- glmnet::cv.glmnet(
    x = x,
    y = y,
    family = "binomial",
    alpha = 0.5,  # Elastic net mixing parameter
    nfolds = 5,
    type.measure = "auc"  # Use AUC for CV metric
  )
  
  # Debug: Look at coefficients
  message("Model coefficients at lambda.min:")
  print(coef(cv_fit, s = "lambda.min"))
  
  # Add scaling parameters to model object
  cv_fit$scale_params <- scale_params
  
  # Calculate CV metrics
  cv_metrics <- list(
    cv_auc = max(cv_fit$cvm),  # Best CV AUC
    cv_auc_se = cv_fit$cvsd[which.max(cv_fit$cvm)],  # SE of best CV AUC
    lambda_best = cv_fit$lambda.min  # Best lambda value
  )
  
  return(list(
    model = cv_fit,
    cv_metrics = cv_metrics
  ))
}

#' Generate Potential Pairs
#'
#' @param ms1 First MSObject
#' @param ms2 Second MSObject
#' @param mz_thresh MZ threshold in ppm
#' @param rt_thresh RT threshold in minutes
#' @return Data frame of potential pairs within thresholds
generate_potential_pairs <- function(ms1, ms2, mz_thresh, rt_thresh) {
  # Extract features from both datasets
  ms1_features <- extract_advanced_features(ms1)
  ms2_features <- extract_advanced_features(ms2)
  
  # Create search ranges for each compound in ms1
  ms1_df <- ms1_features %>%
    dplyr::mutate(
      MZ_lower = MZ - (MZ * mz_thresh / 1e6),
      MZ_upper = MZ + (MZ * mz_thresh / 1e6),
      RT_lower = RT - rt_thresh,
      RT_upper = RT + rt_thresh
    )
  
  # Generate pairs using SQL-style joins for efficiency
  pairs <- dplyr::inner_join(
    ms1_df %>% dplyr::select(Compound_ID, MZ, RT, MZ_lower, MZ_upper, RT_lower, RT_upper),
    ms2_features %>% dplyr::select(Compound_ID, MZ, RT),
    by = character(),
    suffix = c("_1", "_2")
  ) %>%
    dplyr::filter(
      # Apply thresholds
      MZ_2 >= MZ_lower,
      MZ_2 <= MZ_upper,
      RT_2 >= RT_lower,
      RT_2 <= RT_upper
    ) %>%
    dplyr::select(-MZ_lower, -MZ_upper, -RT_lower, -RT_upper)
  
  # Join all features from ms1 and ms2 using the Compound_IDs
  pairs <- pairs %>%
    dplyr::left_join(
      ms1_features %>%
        dplyr::rename_with(~ paste0(.x, "_1"), -Compound_ID) %>%
        dplyr::rename(Compound_ID_1 = Compound_ID),
      by = c("Compound_ID_1", "MZ_1", "RT_1")
    ) %>%
    dplyr::left_join(
      ms2_features %>%
        dplyr::rename_with(~ paste0(.x, "_2"), -Compound_ID) %>%
        dplyr::rename(Compound_ID_2 = Compound_ID),
      by = c("Compound_ID_2", "MZ_2", "RT_2")
    )
  
  # Calculate derived features
  pairs <- pairs %>%
    dplyr::mutate(
      # Basic differences
      delta_RT = RT_1 - RT_2,
      delta_MZ = (MZ_1 - MZ_2) / MZ_1 * 1e6,  # in ppm
      
      # Normalized differences
      rel_delta_RT = delta_RT / rt_thresh,
      rel_delta_MZ = delta_MZ / mz_thresh,
      
      # Combined distance metrics
      euclidean_dist = sqrt(rel_delta_RT^2 + rel_delta_MZ^2),
      manhattan_dist = abs(rel_delta_RT) + abs(rel_delta_MZ),
      
      # Squared differences
      delta_RT_squared = delta_RT^2,
      delta_MZ_squared = delta_MZ^2,
      
      # Log-transformed MZ ratios
      log_MZ_ratio = log(MZ_1 / MZ_2),
      
      # RT position features
      RT_sum = RT_1 + RT_2,
      RT_product = RT_1 * RT_2,
      RT_min = pmin(RT_1, RT_2),
      RT_max = pmax(RT_1, RT_2),
      RT_range = RT_max - RT_min,
      
      # MZ position features
      MZ_sum = MZ_1 + MZ_2,
      MZ_product = MZ_1 * MZ_2,
      MZ_min = pmin(MZ_1, MZ_2),
      MZ_max = pmax(MZ_1, MZ_2),
      MZ_range = MZ_max - MZ_min
    )
  
  # Add rank-based feature calculations if they exist
  if (all(c("RT_rank", "RT_decile") %in% colnames(ms1_features)) &&
      all(c("RT_rank", "RT_decile") %in% colnames(ms2_features))) {
    pairs <- pairs %>%
      dplyr::mutate(
        delta_RT_rank = abs(RT_rank_1 - RT_rank_2),
        same_RT_decileTRUE = as.numeric(RT_decile_1 == RT_decile_2)
      )
  }
  
  # Add MZ rank-based feature calculations if they exist
  if (all(c("MZ_rank", "MZ_decile") %in% colnames(ms1_features)) &&
      all(c("MZ_rank", "MZ_decile") %in% colnames(ms2_features))) {
    pairs <- pairs %>%
      dplyr::mutate(
        delta_MZ_rank = abs(MZ_rank_1 - MZ_rank_2),
        same_MZ_decileTRUE = as.numeric(MZ_decile_1 == MZ_decile_2)
      )
  }
  
  # If we have too many pairs, filter out the obvious non-matches
  if (nrow(pairs) > 100000) {
    pairs <- pairs %>%
      dplyr::filter(euclidean_dist <= 1.5)
  }
  
  return(pairs)
}

#' Apply Semi-supervised Learning
#'
#' @param initial_model Initial trained model
#' @param training_data Initial training data
#' @param all_pairs All potential pairs
#' @param features_ms1 Features from first dataset
#' @param features_ms2 Features from second dataset
#' @param confidence_threshold Confidence threshold for adding predictions to training
#' @param max_iterations Maximum iterations for semi-supervised learning
#' @param seed Random seed for reproducibility
#' @return List with updated model and training data
apply_semisupervised_learning <- function(
  initial_model,
  training_data,
  all_pairs,
  features_ms1,
  features_ms2,
  confidence_threshold,
  max_iterations,
  seed
) {
  
  current_model <- initial_model
  current_training_data <- training_data
  
  # Get the feature names from the initial model
  model_features <- colnames(current_model$glmnet.fit$beta)
  
  # For each iteration
  for (iter in 1:max_iterations) {
    message(sprintf("Semi-supervised learning iteration %d/%d", iter, max_iterations))
    
    # Add features to all pairs
    prediction_data <- all_pairs %>%
      # Remove already used compounds
      dplyr::filter(
        !Compound_ID_1 %in% current_training_data$Compound_ID_1,
        !Compound_ID_2 %in% current_training_data$Compound_ID_2
      )
    
    
    # Join with features
    prediction_data <- prediction_data %>%
      dplyr::left_join(
        features_ms1 %>% 
          dplyr::select(-MZ, -RT) %>%
          dplyr::rename_with(~ paste0(.x, "_1"), -Compound_ID) %>%
          dplyr::rename(Compound_ID_1 = Compound_ID),
        by = "Compound_ID_1"
      ) %>%
      dplyr::left_join(
        features_ms2 %>% 
          dplyr::select(-MZ, -RT) %>%
          dplyr::rename_with(~ paste0(.x, "_2"), -Compound_ID) %>%
          dplyr::rename(Compound_ID_2 = Compound_ID),
        by = "Compound_ID_2"
      )
    
    # Prepare prediction data with all necessary features
    prediction_data <- prepare_prediction_data(prediction_data, model_features)
    
    # Make predictions
    x_pred <- as.matrix(prediction_data)
    predictions <- predict(current_model, newx = x_pred, type = "response", s = "lambda.min")
    prediction_data$prob_matched <- as.vector(predictions)
    
    # Select high-confidence predictions
    high_confidence <- prediction_data %>%
      dplyr::filter(
        prob_matched >= confidence_threshold | 
          prob_matched <= (1 - confidence_threshold)
      ) %>%
      dplyr::mutate(
        Class = factor(ifelse(prob_matched >= confidence_threshold, 
                            "matched", "unmatched"),
                      levels = levels(current_training_data$Class))
      )
    
    if (nrow(high_confidence) == 0) {
      message("No high-confidence predictions found")
      break
    }
    
    # Balance classes
    pos_count <- sum(high_confidence$Class == "matched")
    neg_count <- sum(high_confidence$Class == "unmatched")
    target_count <- min(pos_count, neg_count, 100)
    
    set.seed(seed + iter)
    
    if (pos_count > 0) {
      positive_samples <- high_confidence %>%
        dplyr::filter(Class == "matched") %>%
        dplyr::slice_sample(n = min(pos_count, target_count))
    } else {
      positive_samples <- NULL
    }
    
    if (neg_count > 0) {
      negative_samples <- high_confidence %>%
        dplyr::filter(Class == "unmatched") %>%
        dplyr::slice_sample(n = min(neg_count, target_count))
    } else {
      negative_samples <- NULL
    }
    
    # Combine samples
    new_training <- dplyr::bind_rows(positive_samples, negative_samples)
    
    if (nrow(new_training) > 0) {
      # Update training data
      current_training_data <- dplyr::bind_rows(current_training_data, new_training)
      
      # Re-train model
      current_model <- train_logistic_model(current_training_data, seed + iter)
      
      message(sprintf(
        "Added %d new training examples (total: %d)",
        nrow(new_training),
        nrow(current_training_data)
      ))
    } else {
      message("No new training examples added")
      break
    }
  }
  
  return(list(
    "model" = current_model,
    "training_data" = current_training_data
  ))
}

#' Predict With Constraints
#'
#' @param model Trained model
#' @param all_pairs All potential pairs
#' @param prob_thresh Probability threshold for matches
#' @param rt_thresh RT threshold (unused, kept for backward compatibility)
#' @param mz_thresh MZ threshold (unused, kept for backward compatibility)
#' @param use_intensity Whether to use intensity features (unused, kept for backward compatibility)
#' @return Data frame of predicted matches with constraints applied
predict_with_constraints <- function(model, all_pairs, prob_thresh = 0.5, 
                                   rt_thresh = NULL, mz_thresh = NULL, use_intensity = TRUE) {
  
  # Store ID columns before feature preparation
  id_cols <- all_pairs %>%
    dplyr::select(Compound_ID_1, Compound_ID_2, MZ_1, MZ_2, RT_1, RT_2)
  
  # Get the exact feature names from the model
  model_features <- rownames(model$glmnet.fit$beta)
  
  # Debug: Show model features
  message("Model features:")
  print(model_features)
  
  # Prepare prediction data with the exact features needed by the model
  prediction_data <- prepare_prediction_data(
    all_pairs, 
    model_features = model_features
  )
  
  # Debug: Look at prediction data before scaling
  message("Prediction data summary before scaling:")
  print(summary(prediction_data))
  
  if (nrow(prediction_data) == 0) {
    warning("No pairs to predict")
    return(tibble::tibble())
  }
  
  # Scale prediction data using the same parameters as training
  x_pred <- as.matrix(prediction_data)
  x_pred <- scale(x_pred, 
                 center = model$scale_params$center,
                 scale = model$scale_params$scale)
  
  # Debug: Look at scaled prediction data
  message("Scaled prediction data summary:")
  print(summary(x_pred))
  
  # Make predictions
  linear_predictions <- predict(model, newx = x_pred, s = "lambda.min")
  
  # Debug: Look at predictions
  message("Linear prediction summary:")
  print(summary(linear_predictions))
  
  probabilities <- 1 / (1 + exp(-linear_predictions))
  
  message("Probability summary:")
  print(summary(probabilities))
  
  # Combine predictions with ID columns
  prediction_data_with_ids <- id_cols %>%
    dplyr::mutate(prob_matched = as.vector(probabilities))
  
  # Process candidates
  candidates <- prediction_data_with_ids %>%
    dplyr::filter(prob_matched >= prob_thresh) %>%
    dplyr::arrange(dplyr::desc(prob_matched))
  
  # Apply greedy matching algorithm
  matches <- tibble::tibble()
  used_ids <- tibble::tibble(
    Compound_ID_1 = character(0),
    Compound_ID_2 = character(0)
  )
  
  while (nrow(candidates) > 0) {
    # Get best match
    best_match <- candidates %>%
      dplyr::slice(1)
    
    # Add to matches
    matches <- dplyr::bind_rows(matches, best_match)
    
    # Update used IDs
    used_ids <- dplyr::bind_rows(
      used_ids,
      tibble::tibble(
        Compound_ID_1 = best_match$Compound_ID_1,
        Compound_ID_2 = best_match$Compound_ID_2
      )
    )
    
    # Filter remaining candidates
    candidates <- candidates %>%
      dplyr::filter(
        !Compound_ID_1 %in% used_ids$Compound_ID_1,
        !Compound_ID_2 %in% used_ids$Compound_ID_2
      )
  }
  
  # Return matches sorted by probability
  matches %>%
    dplyr::arrange(dplyr::desc(prob_matched))
}

#' Prepare Prediction Data
#'
#' @param all_pairs All potential pairs with pre-calculated features
#' @param model_features Feature names from model
#' @return Data frame ready for prediction
prepare_prediction_data <- function(all_pairs, model_features) {
  # Select only the features needed by the model
  prediction_data <- all_pairs %>%
    dplyr::select(dplyr::all_of(model_features))
  
  return(prediction_data)
}

#' Compute Match Metrics
#'
#' @param predictions Predicted matches
#' @param true_matches Known true matches
#' @return List of quality metrics
compute_match_metrics <- function(predictions, true_matches) {
  # If no true matches, return NA metrics
  if (nrow(true_matches) == 0) {
    return(list(
      precision = NA_real_,
      recall = NA_real_,
      f1_score = NA_real_
    ))
  }
  
  # Prepare predicted and true match sets
  pred_matches <- paste(predictions$Compound_ID_1, predictions$Compound_ID_2, sep = "_")
  
  # Prepare true match set
  true_matches_set <- paste(true_matches$Compound_ID, true_matches$Compound_ID.y, sep = "_")
  
  # Calculate metrics
  true_positives <- sum(pred_matches %in% true_matches_set)
  false_positives <- sum(!pred_matches %in% true_matches_set)
  false_negatives <- sum(!true_matches_set %in% pred_matches)
  
  precision <- true_positives / (true_positives + false_positives)
  recall <- true_positives / (true_positives + false_negatives)
  
  f1_score <- if (precision + recall > 0) {
    2 * (precision * recall) / (precision + recall)
  } else {
    0
  }
  
  return(list(
    precision = precision,
    recall = recall,
    f1_score = f1_score
  ))
}

#' Plot Feature Importance
#'
#' @param ml_match_result Result from ml_match function
#' @param top_n Number of top features to show
#' @return ggplot object showing feature importance
#' @export
plot_feature_importance <- function(ml_match_result, top_n = 20) {
  # Extract importance measures
  importance <- randomForest::importance(ml_match_result$model$finalModel)
  imp_df <- data.frame(
    Feature = rownames(importance),
    Importance = importance[, "MeanDecreaseGini"]
  )
  
  # Sort and select top features
  imp_df <- imp_df %>%
    dplyr::arrange(desc(Importance)) %>%
    head(top_n)
  
  # Create plot
  p <- ggplot2::ggplot(imp_df, ggplot2::aes(x = reorder(Feature, Importance), y = Importance)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Feature Importance",
      x = "Feature",
      y = "Importance (Mean Decrease in Gini)"
    )
  
  return(p)
}

#' Plot Decision Boundary
#'
#' @param ml_match_result Result from ml_match function
#' @return ggplot object showing decision boundary
#' @export
plot_decision_boundary <- function(ml_match_result) {
  # Extract matched data
  matched_data <- ml_match_result$matched_data
  
  # Get range for delta_RT and delta_MZ
  delta_rt_range <- range(matched_data$delta_RT)
  delta_mz_range <- range(matched_data$delta_MZ)
  
  # Create grid for prediction
  grid_size <- 100
  pred_grid <- expand.grid(
    delta_RT = seq(delta_rt_range[1], delta_rt_range[2], length.out = grid_size),
    delta_MZ = seq(delta_mz_range[1], delta_mz_range[2], length.out = grid_size)
  )
  
  # Add mean values for other features
  mean_features <- ml_match_result$training_data %>%
    dplyr::select(-delta_RT, -delta_MZ, -Class, -Compound_ID_1, -Compound_ID_2, 
                  -Metabolite_1, -Metabolite_2) %>%
    dplyr::summarise(across(everything(), mean, na.rm = TRUE))
  
  # Replicate means for grid
  for (col in names(mean_features)) {
    pred_grid[[col]] <- mean_features[[col]]
  }
  
  # Prepare for prediction
  pred_grid <- prepare_prediction_data(pred_grid, ml_match_result$model$finalModel$xNames)
  
  # Make predictions
  pred_grid$prob <- stats::predict(ml_match_result$model, pred_grid, type = "prob")$matched
  
  # Create contour plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = pred_grid, 
                      ggplot2::aes(x = delta_RT, y = delta_MZ, fill = prob)) +
    ggplot2::geom_contour(data = pred_grid,
                         ggplot2::aes(x = delta_RT, y = delta_MZ, z = prob),
                         breaks = c(0.5)) +
    ggplot2::geom_point(data = matched_data,
                       ggplot2::aes(x = delta_RT, y = delta_MZ),
                       color = "black", size = 2) +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                                 midpoint = 0.5, limits = c(0, 1)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Decision Boundary",
      x = "Delta RT (minutes)",
      y = "Delta MZ (ppm)",
      fill = "Match Probability"
    )
  
  return(p)
}

#' ML Match without Semi-supervised Learning
#'
#' @param ms1 First MSObject
#' @param ms2 Second MSObject
#' @param mz_thresh MZ threshold in ppm
#' @param rt_thresh RT threshold in minutes
#' @param prob_thresh Probability threshold for matches
#' @param neg_examples_per_pos Number of negative examples per positive
#' @param use_intensity Whether to use intensity features
#' @param n_folds Number of CV folds
#' @param seed Random seed
#' @return List containing model and predictions
ml_match_simple <- function(ms1, 
                           ms2, 
                           mz_thresh = 15, 
                           rt_thresh = 1,
                           prob_thresh = 0.5,
                           neg_examples_per_pos = 3,
                           use_intensity = TRUE,
                           n_folds = 3,
                           seed = 72) {
  
  # Get shared metabolites for training
  message("Getting shared metabolites")
  labeled_data <- get_shared_metabolites(ms1, ms2)
  
  # Extract features
  message("Extracting features")
  features_ms1 <- extract_advanced_features(ms1)
  features_ms2 <- extract_advanced_features(ms2)
  
  # Create training data
  message("Creating training data")
  training_data <- create_enhanced_training_data(
    labeled_data, 
    features_ms1, 
    features_ms2, 
    mz_thresh, 
    rt_thresh,
    neg_examples_per_pos,
    use_intensity
  )
  
  # Train model with cross-validation
  message("Training model with cross-validation")
  model <- train_logistic_model(training_data, seed)
  
  # Generate all potential pairs for prediction
  message("Generating potential pairs")
  all_pairs <- generate_potential_pairs(ms1, ms2, mz_thresh, rt_thresh)
  
  # Make predictions with constraints
  message("Making final predictions")
  predictions <- predict_with_constraints(
    model, 
    all_pairs, 
    prob_thresh, 
    rt_thresh = rt_thresh,
    mz_thresh = mz_thresh,
    use_intensity = use_intensity
  )
  
  # Calculate metrics if we have known matches
  metrics <- compute_match_metrics(predictions, labeled_data$known)
  
  # Print results
  message(sprintf(
    "Found %d matches with %.1f%% of known matches correctly identified",
    nrow(predictions),
    metrics$recall * 100
  ))
  
  return(list(
    model = model,
    predictions = predictions,
    metrics = metrics,
    training_data = training_data
  ))
}