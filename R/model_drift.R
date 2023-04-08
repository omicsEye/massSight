model_drift <- function(aligned_ms_obj, smooth_method) {
  rt <- iso_matched(aligned_ms_obj)$RT
  delta_rt <- iso_matched(aligned_ms_obj)$RT_2 - iso_matched(aligned_ms_obj)$RT
  mz <- iso_matched(aligned_ms_obj)$MZ
  delta_mz <- iso_matched(aligned_ms_obj)$MZ_2 - iso_matched(aligned_ms_obj)$MZ

  if (smooth_method == "loess") {
    rt_smooth <- loess(delta_rt ~ rt)
    mz_smooth <- loess(delta_mz ~ mz)
  } else if (smooth_method == "gam") {
    rt_smooth <- mgcv::gam(delta_rt ~ rt)
    mz_smooth <- mgcv::gam(delta_mz ~ mz)
  }

  delta_rt_smooth <- predict(rt_smooth, rt)
  delta_mz_smooth <- predict(rt_smooth, mz)
  iso_matched(aligned_ms_obj) <- iso_matched(aligned_ms_obj) |>
    dplyr::mutate(
      delta_rt = delta_rt,
      delta_mz = delta_mz,
      delta_rt_smooth = delta_rt_smooth,
      delta_mz_smooth = delta_mz_smooth
    )

  smooth(aligned_ms_obj) <- list(
    "rt_smooth" = rt_smooth,
    "mz_smooth" = mz_smooth
  )

  return(aligned_ms_obj)
}
