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

  aligned <- auto_combine(
    ms1,
    ms2,
    log = NULL
  )
  out <- all_matched(aligned)
  expect_equal(out$df1, out$df2)
})
