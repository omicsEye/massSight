# create_base_dataset <- function(n = 15000, seed = 123) {
#   set.seed(123)
#   base_dataset <- dplyr::tibble(
#     Compound_ID = paste0("Compound_", 1:n),
#     RT = runif(n, 1, 10),
#     MZ = runif(n, 100, 700),
#     Metabolite = paste0("Metabolite_", 1:n),
#     Intensity = rep(2000, n)
#   )
#   return(base_dataset)
# }
#
# add_constant_drift <- function(base,
#                                mz_beta = 5,
#                                rt_beta = .1,
#                                mz_eps = 1,
#                                rt_eps = .05) {
#   n <- nrow(base)
#   mz_add <- rnorm(n, mz_beta, mz_eps)
#   rt_add <- rnorm(n, rt_beta, rt_eps)
#   drift_dataset <- base |>
#     dplyr::mutate(MZ = MZ + mz_add * MZ / 1e6,
#                   RT = RT + rt_add)
#   return(drift_dataset)
# }
#
# add_linear_drift <- function(base,
#                              mz_beta = 5,
#                              rt_beta = .2,
#                              mz_eps = 2,
#                              rt_eps = .2) {
#   n <- nrow(base)
#   mz_add <- rep(mz_beta, n) * base$MZ + rnorm(n, 0, mz_eps)
#   rt_add <- rep(rt_beta) * base$RT + rnorm(n, 0, rt_eps)
#   drift_dataset <- base |>
#     dplyr::mutate(MZ = MZ + mz_add * MZ / 1e6,
#                   RT = RT + rt_add)
#   return(drift_dataset)
# }
#
# base_dataset <- create_base_dataset()
# drift_dataset <- add_constant_drift(base_dataset)
# lin_drift_dataset <- add_linear_drift(base_dataset)
#
# ms1 <-
#   create_ms_obj(
#     df = base_dataset,
#     name = "base",
#     id_name = "Compound_ID",
#     rt_name = "RT",
#     mz_name = "MZ",
#     int_name = "Intensity"
#   )
#
# ms2 <-
#   create_ms_obj(
#     df = lin_drift_dataset,
#     name = "drift",
#     id_name = "Compound_ID",
#     rt_name = "RT",
#     mz_name = "MZ",
#     int_name = "Intensity"
#   )
#
# out <-
#   auto_combine(
#     ms1,
#     ms2,
#     rt_lower = -.5,
#     rt_upper = .5,
#     mz_lower = -10,
#     mz_upper = 10
#   )
#
# final_plots(out)
#
# out@all_matched$df1 == out@all_matched$df2
#
# samp <- sample(1:nrow(base_dataset), size = 14900)
# base_dataset$Metabolite[samp] <- NA
# lin_drift_dataset$Metabolite[samp] <- NA
#
# ms1 <-
#   create_ms_obj(
#     df = base_dataset,
#     name = "base",
#     id_name = "Compound_ID",
#     rt_name = "RT",
#     mz_name = "MZ",
#     int_name = "Intensity"
#   )
#
# ms2 <-
#   create_ms_obj(
#     df = lin_drift_dataset,
#     name = "drift",
#     id_name = "Compound_ID",
#     rt_name = "RT",
#     mz_name = "MZ",
#     int_name = "Intensity"
#   )
#
# out <-
#   auto_combine(
#     ms1,
#     ms2,
#     rt_lower = -.5,
#     rt_upper = .5,
#     mz_lower = -10,
#     mz_upper = 10,
#     match_method = "supervised"
#   )
#
# final_plots(out)
#
#
# base_dataset <- base_dataset |>
#   arrange(RT)
#
# lin_drift_dataset <- lin_drift_dataset |>
#   arrange(RT)
#
# drift <- lin_drift_dataset$RT - base_dataset$RT
#
# plot(base_dataset$RT, drift)
#
#
# for (i in 1:1000) {
#   base <- sample(base_dataset$RT, )
# }
#
