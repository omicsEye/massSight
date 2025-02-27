test_that("mass_combine creates a MergedMSObject", {
  # Create test MSObjects
  ms1 <- create_ms_obj(
    df = hp1,
    name = "hp1",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity",
    metab_name = "Metabolite"
  )
  ms2 <- create_ms_obj(
    df = hp2,
    name = "hp2",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity",
    metab_name = "Metabolite"
  )
  
  # Test with default parameters
  result <- mass_combine(ms1, ms2, optimize = FALSE)
  
  # Basic checks
  expect_s4_class(result, "MergedMSObject")
  expect_s4_class(ms1(result), "MSObject")
  expect_s4_class(ms2(result), "MSObject")
  expect_equal(name(ms1(result)), "hp1")
  expect_equal(name(ms2(result)), "hp2")
})

test_that("mass_combine with optimization returns valid results", {
  # Create test MSObjects
  ms1 <- create_ms_obj(
    df = hp1,
    name = "hp1",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity",
    metab_name = "Metabolite"
  )
  ms2 <- create_ms_obj(
    df = hp2,
    name = "hp2",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity",
    metab_name = "Metabolite"
  )
  
  # Test with optimization (limiting iterations for test speed)
  result <- mass_combine(ms1, ms2, optimize = TRUE, n_iter = 2)
  
  # Check optimization attributes
  expect_true(!is.null(attr(result, "optimization")))
  expect_true(is.list(attr(result, "optimization")))
  expect_true("parameters" %in% names(attr(result, "optimization")))
  expect_true("final_score" %in% names(attr(result, "optimization")))
  
  # Check that optimized parameters are within expected ranges
  opt_params <- attr(result, "optimization")$parameters
  expect_true(opt_params[["rt_delta"]] >= 0.1 & opt_params[["rt_delta"]] <= 1.0)
  expect_true(opt_params[["mz_delta"]] >= 1 & opt_params[["mz_delta"]] <= 20)
})

test_that("mass_combine handles different isolation methods correctly", {
  # Create test MSObjects
  ms1 <- create_ms_obj(
    df = hp1,
    name = "hp1",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity",
    metab_name = "Metabolite"
  )
  ms2 <- create_ms_obj(
    df = hp2,
    name = "hp2",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity",
    metab_name = "Metabolite"
  )
  
  # Test with manual isolation method
  result_manual <- mass_combine(
    ms1, ms2, 
    optimize = FALSE, 
    iso_method = "manual", 
    rt_iso_threshold = 0.02,
    mz_iso_threshold = 3
  )
  
  # Test with DBSCAN isolation method
  result_dbscan <- mass_combine(
    ms1, ms2, 
    optimize = FALSE, 
    iso_method = "dbscan", 
    eps = 0.15
  )
  
  # Both should return valid MergedMSObjects
  expect_s4_class(result_manual, "MergedMSObject")
  expect_s4_class(result_dbscan, "MergedMSObject")
})

test_that("mass_combine handles different matching methods correctly", {
  # Create test MSObjects
  ms1 <- create_ms_obj(
    df = hp1,
    name = "hp1",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity",
    metab_name = "Metabolite"
  )
  ms2 <- create_ms_obj(
    df = hp2,
    name = "hp2",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity",
    metab_name = "Metabolite"
  )
  
  # Test with unsupervised matching
  result_unsup <- mass_combine(
    ms1, ms2, 
    optimize = FALSE, 
    match_method = "unsupervised"
  )
  
  # Test with supervised matching
  result_sup <- mass_combine(
    ms1, ms2, 
    optimize = FALSE, 
    match_method = "supervised"
  )
  
  # Both should return valid MergedMSObjects
  expect_s4_class(result_unsup, "MergedMSObject")
  expect_s4_class(result_sup, "MergedMSObject")
})

test_that("mass_combine handles different smoothing methods correctly", {
  skip_on_cran() # Skip on CRAN due to computation time
  
  # Create test MSObjects
  ms1 <- create_ms_obj(
    df = hp1,
    name = "hp1",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity",
    metab_name = "Metabolite"
  )
  ms2 <- create_ms_obj(
    df = hp2,
    name = "hp2",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity",
    metab_name = "Metabolite"
  )
  
  # Test with GAM smoothing
  result_gam <- mass_combine(
    ms1, ms2, 
    optimize = FALSE, 
    smooth_method = "gam"
  )
  
  # Test with linear model smoothing
  result_lm <- mass_combine(
    ms1, ms2, 
    optimize = FALSE, 
    smooth_method = "lm"
  )
  
  # Both should return valid MergedMSObjects
  expect_s4_class(result_gam, "MergedMSObject")
  expect_s4_class(result_lm, "MergedMSObject")
})

test_that("mass_combine validates input parameters correctly", {
  # Create test MSObjects
  ms1 <- create_ms_obj(
    df = hp1,
    name = "hp1",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity",
    metab_name = "Metabolite"
  )
  ms2 <- create_ms_obj(
    df = hp2,
    name = "hp2",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity",
    metab_name = "Metabolite"
  )
  
  # Test with invalid iso_method
  expect_error(
    mass_combine(ms1, ms2, optimize = FALSE, iso_method = "invalid_method"),
    "iso_method must be one of: manual, dbscan"
  )
  
  # Test with invalid match_method
  expect_error(
    mass_combine(ms1, ms2, optimize = FALSE, match_method = "invalid_method"),
    "match_method must be one of: unsupervised, supervised"
  )
  
  # Test with invalid smooth_method
  expect_error(
    mass_combine(ms1, ms2, optimize = FALSE, smooth_method = "invalid_method"),
    "smooth_method must be one of: gam, bayesian_gam, lm, gp"
  )
  
  # Test with invalid minimum_intensity
  expect_error(
    mass_combine(ms1, ms2, optimize = FALSE, minimum_intensity = -10),
    "minimum_intensity must be greater than 0"
  )
})

test_that("self-alignment ID's match", {
  ms1 <- create_ms_obj(
    df = hp2,
    name = "hp2_1",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity",
    metab_name = "Metabolite"
  )
  ms2 <-
    create_ms_obj(
      df = hp2,
      name = "hp2_2",
      id_name = "Compound_ID",
      rt_name = "RT",
      mz_name = "MZ",
      int_name = "Intensity",
      metab_name = "Metabolite"
    )

  aligned <- mass_combine(
    ms1,
    ms2,
    log = NULL
  )
  
  # Use dynamic column names based on object names
  col1 <- paste0("Compound_ID_", name(ms1))
  col2 <- paste0("Compound_ID_", name(ms2))
  
  # Get matched compounds and check they're identical
  matches <- get_unique_matches(aligned)
  expect_equal(matches[[col1]], matches[[col2]])
})

test_that("all compound IDs from ms1 and ms2 are in the final merged object", {
  # Create two massSight objects
  ms1 <- create_ms_obj(
    df = hp1,
    name = "hp1",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity",
    metab_name = "Metabolite"
  )
  ms2 <- create_ms_obj(
    df = hp2,
    name = "hp2",
    id_name = "Compound_ID",
    rt_name = "RT",
    mz_name = "MZ",
    int_name = "Intensity",
    metab_name = "Metabolite"
  )

  # Perform mass_combine
  aligned <- mass_combine(ms1, ms2, log = NULL)

  # Extract compound IDs from ms1, ms2, and the aligned object
  ms1_ids <- ms1@df$Compound_ID
  ms2_ids <- ms2@df$Compound_ID
  
  # Get merged results
  merged_results <- all_matched(aligned)
  
  # Check column names in merged results
  col1 <- paste0("Compound_ID_", name(ms1))
  col2 <- paste0("Compound_ID_", name(ms2))
  
  # Check that all IDs from ms1 and ms2 are in the merged results
  expect_true(all(ms1_ids %in% merged_results[[col1]] | is.na(merged_results[[col1]])))
  expect_true(all(ms2_ids %in% merged_results[[col2]] | is.na(merged_results[[col2]])))
  
  # Count total unique compounds
  total_compounds <- length(unique(c(ms1_ids, ms2_ids)))
  merged_compounds <- nrow(merged_results)
  
  # The merged results should contain at least the total unique compounds
  expect_gte(merged_compounds, total_compounds)
})