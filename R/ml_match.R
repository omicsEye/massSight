get_shared_metabolites <- function(ms1, ms2) {
  ms1_df <- raw_df(ms1)
  ms2_df <- raw_df(ms2)
  ms1_known <- raw_df(ms1) |>
    filter(Metabolite != "" &
             Compound_ID != Metabolite)
  ms2_known <- raw_df(ms2) |>
    filter(Metabolite != "" &
             Compound_ID != Metabolite)
  known <- inner_join(ms1_known, ms2_known, by = "Metabolite") |>
    select(MZ.x, MZ.y, RT.x, RT.y, Metabolite) |>
    rename(
      MZ_1 = MZ.x,
      MZ_2 = MZ.y,
      RT_1 = RT.x,
      RT_2 = RT.y,
      Metabolite_1 = Metabolite
    ) |>
    mutate(Class = "matched",
           Metabolite_2 = Metabolite)
  return(known)
}


fit_model <- function(ms1, ms2, seed = 72) {
  known <- get_shared_metabolites(ms1, ms2)

  for (i in 1:nrow(known)) {
    ms1_sample <- ms1_known |> slice_sample(n = 1)
    ms2_sample <- ms2_known |>
      filter(Metabolite != ms1_sample[1, "Metabolite"]) |>
      slice_sample(n = 1)
    row <- data.frame(
      MZ_1 = ms1_sample$MZ,
      MZ_2 = ms2_sample$MZ,
      RT_1 = ms1_sample$RT,
      RT_2 = ms2_sample$RT,
      Metabolite_1 = ms1_sample$Metabolite,
      Metabolite_2 = ms2_sample$Metabolite,
      Class = "unmatched"
    )
    known <- bind_rows(known, row)
  }
  known <- known |>
    mutate(
      Class = as.factor(Class),
      delta_RT = RT.x - RT.y,
      delta_MZ = MZ.x - MZ.y
    )

  set.seed(seed)
}
#
#
# set.seed(1234)
# rf <- train(
#   class ~ delta_RT + delta_MZ + RT.x + RT.y + MZ.x + MZ.y,
#   data = known,
#   method = "rf",
#   ntree = 501,
#   tuneGrid = data.frame(mtry = 1:6),
#   trControl = trainControl(method = "cv", number = 10)
# )
#
# rf$finalModel$err.rate |>
#   as_tibble() |>
#   ggplot(aes(x = seq_along(OOB), y = OOB)) +
#   geom_point() +
#   geom_smooth(se = F) +
#   labs(x = "# Trees")
#
#
# b_test <- b |>
#   filter(Metabolite == "")
#
# a_test <- a |>
#   filter(Metabolite == "" |
#            Metabolite == Compound_ID)
#
# matched <- data.frame(
#   "matched_a" = character(),
#   "matched_b" = character(),
#   "prob" = numeric()
# )
# pb <- progress::progress_bar$new(total = nrow(a_test),
#                                  format = "matching [:bar] :percent eta: :eta")
# for (row in 1:nrow(a_test)) {
#   query <- a_test[row, ] |>
#     dplyr::slice(rep(1:n(), each = nrow(b_test)))
#   full <- data.frame(
#     MZ.x = query$MZ,
#     RT.x = query$RT,
#     MZ.y = b_test$MZ,
#     RT.y = b_test$RT
#   ) |>
#     mutate(delta_RT = RT.x - RT.y,
#            delta_MZ = MZ.x - MZ.y)
#   out <- predict(rf$finalModel, full, type = "prob")
#   matched <- matched |>
#     bind_rows(
#       data.frame(
#         "matched_a" = a_test[row, "Compound_ID"],
#         "MZ_a" = a_test[row, "MZ"],
#         "RT_a" = a_test[row, "RT"],
#         "matched_b" = b_test[which.max(out[, 1]), "Compound_ID"],
#         "MZ_b" = b_test[which.max(out[, 1]), "MZ"],
#         "RT_b" = b_test[which.max(out[, 1]), "RT"],
#         "prob" = max(out[, 1])
#       )
#     )
#   pb$tick()
# }
