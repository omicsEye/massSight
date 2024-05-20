#' @export
#' @title Create MS Object
#' @description Create an MSObject from a data frame.
#'
#' @param df A data frame containing the raw data.
#' @param name A character indicating the name of the experiment.
#' @param id_name A character indicating the name of the column containing the
#' compound IDs.
#' @param rt_name A character indicating the name of the column containing the
#' retention times.
#' @param mz_name A character indicating the name of the column containing the
#' m/z values.
#' @param int_name A character indicating the name of the column containing the
#' intensities
#' @param metab_name An optional character indicating the name of the column
#' containing the metabolite annotations
#' @importFrom data.table :=
#' @return An MSObject.
create_ms_obj <- function(df,
                          name,
                          id_name = "Compound_ID",
                          rt_name = "RT",
                          mz_name = "MZ",
                          int_name = "Intensity",
                          metab_name = "Metabolite") {
  compound_ids <- df[[id_name]]
  stopifnot("Compound IDs must be unique." = !any(duplicated(compound_ids))) |>
    try()
  ms <- methods::new("MSObject")
  name(ms) <- name
  consolidated(ms) <- FALSE
  if (typeof(df) != "data.table") {
    df = data.table::as.data.table(df)
  }
  raw_data <- data.table::setDT(df)[, c(
    Compound_ID = data.table::.SD[[id_name]],
    Metabolite = data.table::.SD[[metab_name]],
    RT = data.table::.SD[[rt_name]],
    MZ = data.table::.SD[[mz_name]],
    Intensity = data.table::.SD[[int_name]]
  ), .SDcols = c(id_name, metab_name, rt_name, mz_name, int_name)]
  raw_df(ms) <- raw_data

  meta_data <- data.table::copy(raw_data)
  meta_data[, c(metab_name, rt_name, mz_name, int_name) := NULL]

  if (ncol(meta_data) != 1) {
    metadata(ms) <- meta_data
  }
  return(ms)
}

create_aligned_ms_obj <- function(ms1, ms2) {
  ms <- methods::new(
    "MergedMSObject",
    ms1 = ms1,
    ms2 = ms2,
    all_matched = NULL,
    iso_matched = NULL,
    scaled_values = NULL,
    adjusted_df = NULL,
    cutoffs = NULL,
    aligned = NULL,
    metadata = NULL,
    smooth_method = NULL
  )
  return(ms)
}
