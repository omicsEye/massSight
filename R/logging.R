check_jsonlite <- function() {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    message("The 'jsonlite' package is required for logging.")
    message("1: Install 'jsonlite")
    message("2: Do not install and exit")
    choice <- readline(prompt = "Enter your choice (1 or 2): ")

    if (choice == "1") {
      install.packages("jsonlite")
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("Package 'jsonlite' failed to install. Please try installing manually.",
             call. = FALSE)
      }
    } else if (choice == "2") {
      stop("Package 'jsonlite' is not installed.", call. = FALSE)
    } else {
      stop("Invalid choice. Please restart the function and enter 1 or 2.",
           call. = FALSE)
    }
  }
}

initialize_log <- function(call) {
  logr::log_open(paste0(format(Sys.time(), "%Y-%m-%d_%H-%M"), ".log"),
                 traceback = FALSE,
                 show_notes = FALSE)

  header_text <- "massSight Run and Parameters "
  rep_count <- 80 - nchar(header_text)
  header <- paste0(header_text, strrep("-", rep_count))
  logr::log_print(header, console = F, hide_notes = T)
  logr::log_print(call, console = F, hide_notes = T)
}


modify_call <- function(call) {
  defaults <- as.list(formals(mass_combine))
  combined <- modifyList(defaults, as.list(call)[-1])
  call <- as.call(c(as.symbol(deparse(call[[1]])), combined))
  return(call)
}
