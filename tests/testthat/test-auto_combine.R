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
  out <- all_matched(aligned)
  expect_equal(out$df1, out$df2)
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
  aligned_ids <- aligned@all_matched$ref_Compound_ID

  # Check if all IDs from ms1 and ms2 are in the aligned object
  expect_true(all(ms1_ids %in% aligned_ids))
  expect_true(all(ms2_ids %in% aligned_ids))

  # Check if the number of unique IDs in aligned object matches the total unique IDs from ms1 and ms2
  expected_unique_ids <- length(unique(c(ms1_ids, ms2_ids)))
  actual_unique_ids <- length(unique(aligned_ids))
  expect_equal(actual_unique_ids, expected_unique_ids)
})
