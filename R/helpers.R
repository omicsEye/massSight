validate_parameters <-
  function(iso_method,
           match_method,
           smooth_method,
           minimum_intensity) {
    checkmate::assert_choice(smooth_method, c("loess", "gam"), .var.name = "smooth_method")
    checkmate::assert_choice(iso_method, c("manual", "dbscan"), .var.name = "iso_method")
    checkmate::assert_choice(match_method, c("supervised", "unsupervised"), .var.name = "match_method")
    checkmate::assert_numeric(minimum_intensity)
  }
