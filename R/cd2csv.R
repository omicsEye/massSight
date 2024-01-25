#' @title CD-To-Csv
#' @description
#' Converts output from Compound Discoverer into a csv format suitable
#' for `massSight`
#'
#'
#' @param path The path to the Compound Discoverer Excel (.xlsx) file
#' @param gen_id If `TRUE`, this function will generate a unique id
#' for each compound based on retention time and mass to charge ratio.
#' @param output_file The name of the file to save the formatted data. If
#' `NULL` (default), a csv file will not be written.
#'
#' @return A dataframe that can then be used by `massSight`'s `auto_combine()`
#' or `auto_scale()`
#' @export
cd2csv <- function(path, gen_id = TRUE, output_file = NULL) {
  df <- suppressMessages(readxl::read_excel(path))

  blue_ind <- which(df$Checked)
  orange_ind <- which(df$Name == "TRUE")

  df_blue <- df |>
    dplyr::slice(blue_ind) |>
    dplyr::mutate(
      `Calc. MW` = as.numeric(`Calc. MW`),
      `m/z` = as.numeric(`m/z`),
      `RT [min]` = as.numeric(`RT [min]`),
      `Area (Max.)` = as.numeric(`Area (Max.)`)
    )
  df_orange <- df |> dplyr::slice(orange_ind)
  new_col <- unique(df_orange |> dplyr::pull(7))

  blue_ind <- append(blue_ind, Inf)

  list_of_dfs <- make_df_list(blue_ind, orange_ind, orange_df, df)

  orange_final_df <- dplyr::bind_rows(list_of_dfs)
  final_df <- dplyr::bind_cols(df_blue, orange_final_df)
  if (gen_id) {
    final_df <- final_df |>
      dplyr::mutate(Compound_ID = paste(round(`m/z`, 2),
        round(`RT [min]`, 2),
        sep = "_"
      )) |>
      dplyr::select(-c(7, 8, 9)) |>
      dplyr::select(Compound_ID, dplyr::everything())
  }
  if (!is.null(output_file)) {
    readr::write_csv(final_df, output_file)
  }
  return(final_df)
}

make_df_list <- function(blue_ind, orange_ind, orange_df, df) {
  list_of_dfs <- purrr::map(
    1:(length(blue_ind) - 1),
    ~ {
      orange_df <- df[orange_ind[orange_ind < blue_ind[.x + 1] &
        orange_ind > blue_ind[.x]], c(5, 7)]
      orange_df[, 1] <- as.numeric(orange_df[[1]])
      orange_df |>
        tidyr::pivot_wider(
          names_from = 2,
          values_from = 1,
          values_fn = mean
        )
    }
  )
  return(list_of_dfs)
}
