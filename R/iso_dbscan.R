iso_dbscan <- function(df, eps) {
  dbscan_out <- df |>
    dplyr::select(.data$RT, .data$MZ) |>
    dbscan::dbscan(eps = eps, minPts = 1)

  singleton_clusters <- which(dbscan_out$cluster |> table() == 1) |>
    as.vector()

  df <- df[dbscan_out$cluster %in% singleton_clusters, ]
  return(df)
}
