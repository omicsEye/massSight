initialize_log <- function(call) {
  logr::log_open(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M"),
                        ".log"),
                 traceback = FALSE,
                 show_notes = FALSE)

  header_text <- "massSight Run and Parameters "
  rep_count <- 80 - nchar(header_text)
  header <- paste0(header_text, strrep("-", rep_count))
  logr::log_print(header,
                  console = F,
                  hide_notes = T)
  logr::log_print(call, console = F, hide_notes = T)
}
