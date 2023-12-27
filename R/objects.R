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

  raw_data <- df |>
    dplyr::select(dplyr::all_of(c(id_name, metab_name, rt_name, mz_name, int_name))) |>
    dplyr::rename(
      Compound_ID = id_name,
      Metabolite = metab_name,
      RT = rt_name,
      MZ = mz_name,
      Intensity = int_name
    )
  raw_df(ms) <- raw_data

  meta_data <- df |>
    dplyr::select(-dplyr::any_of(c(rt_name, mz_name, int_name, metab_name))) |>
    dplyr::rename(Compound_ID = id_name)
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
