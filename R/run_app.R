#' Launch massSight Shiny App
#'
#' @description
#' Launches a Shiny web application that provides a graphical user interface for
#' the mass_combine() function. This allows users to interactively align and
#' analyze LC-MS data without writing code.
#'
#' @return Launches a Shiny application in the user's default web browser
#' @importFrom shiny fluidPage titlePanel sidebarLayout sidebarPanel mainPanel
#' @export
run_massSight_app <- function() {
  # Check if required packages are installed
  rlang::check_installed(c("shiny", "DT", "shinyjs", "ggplot2", "ggExtra", "cowplot"))

  # Define UI
  ui <- shiny::fluidPage(
    shinyjs::useShinyjs(),
    # Enable shinyjs
    shiny::titlePanel("massSight - LC-MS Data Alignment"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        # Dataset 1 section
        shiny::h4("Dataset 1"),
        shiny::fileInput("file1", "Upload CSV file"),
        shiny::textInput("name1", "Dataset Name:", value = "dataset1"),
        shiny::h5("Column Mapping"),
        shiny::selectInput(
          "id_name1",
          "Compound ID Column:",
          choices = c("Select Column" = ""),
          selected = ""
        ),
        shiny::selectInput(
          "rt_name1",
          "Retention Time Column:",
          choices = c("Select Column" = ""),
          selected = ""
        ),
        shiny::selectInput(
          "mz_name1",
          "M/Z Column:",
          choices = c("Select Column" = ""),
          selected = ""
        ),
        shiny::selectInput(
          "int_name1",
          "Intensity Column:",
          choices = c("Select Column" = ""),
          selected = ""
        ),
        shiny::selectInput(
          "metab_name1",
          "Metabolite Name Column:",
          choices = c("Select Column" = ""),
          selected = ""
        ),

        shiny::hr(),

        # Dataset 2 section
        shiny::h4("Dataset 2"),
        shiny::fileInput("file2", "Upload CSV file"),
        shiny::textInput("name2", "Dataset Name:", value = "dataset2"),
        shiny::h5("Column Mapping"),
        shiny::selectInput(
          "id_name2",
          "Compound ID Column:",
          choices = c("Select Column" = ""),
          selected = ""
        ),
        shiny::selectInput(
          "rt_name2",
          "Retention Time Column:",
          choices = c("Select Column" = ""),
          selected = ""
        ),
        shiny::selectInput(
          "mz_name2",
          "M/Z Column:",
          choices = c("Select Column" = ""),
          selected = ""
        ),
        shiny::selectInput(
          "int_name2",
          "Intensity Column:",
          choices = c("Select Column" = ""),
          selected = ""
        ),
        shiny::selectInput(
          "metab_name2",
          "Metabolite Name Column:",
          choices = c("Select Column" = ""),
          selected = ""
        ),

        shiny::hr(),

        # Create MS Objects button
        shiny::actionButton(
          "create_ms",
          "Create MS Objects",
          class = "btn-primary",
          style = "width: 100%; margin-top: 10px;"
        ),

        # Analysis parameters section
        shiny::conditionalPanel(
          condition = "input.create_ms > 0",
          shiny::hr(),
          shiny::h4("Analysis Parameters"),
          shiny::checkboxInput("optimize", "Use Parameter Optimization", value = TRUE),

          # Basic parameters when not optimizing
          shiny::conditionalPanel(
            condition = "!input.optimize",
            shiny::numericInput(
              "rt_delta",
              "RT Delta (min):",
              value = 0.5,
              min = 0.1,
              max = 1.0
            ),
            shiny::numericInput(
              "mz_delta",
              "M/Z Delta (ppm):",
              value = 15,
              min = 1,
              max = 20
            ),
            shiny::numericInput(
              "minimum_intensity",
              "Minimum Intensity:",
              value = 1000,
              min = 0
            ),
            shiny::selectInput(
              "iso_method",
              "Isolation Method:",
              choices = c("manual", "dbscan"),
              selected = "manual"
            ),
            shiny::numericInput(
              "rt_iso_threshold",
              "RT Isolation Threshold:",
              value = 0.01,
              min = 0.001,
              max = 0.1
            ),
            shiny::numericInput(
              "mz_iso_threshold",
              "MZ Isolation Threshold:",
              value = 2,
              min = 0.1,
              max = 5
            )
          ),

          # Optimization parameters
          shiny::conditionalPanel(
            condition = "input.optimize",
            shiny::numericInput(
              "n_iter",
              "Number of Optimization Iterations:",
              value = 50,
              min = 10,
              max = 200
            )
          ),

          # Common parameters
          shiny::selectInput(
            "match_method",
            "Matching Method:",
            choices = c("supervised", "unsupervised"),
            selected = "unsupervised"
          ),
          shiny::selectInput(
            "smooth_method",
            "Smoothing Method:",
            choices = c("gam", "bayesian_gam", "gp", "lm"),
            selected = "gam"
          ),

          # Run analysis button
          shiny::actionButton("run", "Run Analysis", class = "btn-primary", style = "width: 100%; margin-top: 10px;")
        )
      ),

      shiny::mainPanel(
        width = 9,
        shiny::tabsetPanel(
          shiny::tabPanel(
            "Upload & Parameters",
            shiny::h4("Upload your datasets and configure parameters"),
            shiny::p(
              "Once you've uploaded both datasets and selected the appropriate columns, click 'Create MS Objects' to proceed."
            )
          ),
          shiny::tabPanel("Data Preview", DT::DTOutput("dataPreview")),
          shiny::tabPanel(
            "Distribution Plots",
            shiny::plotOutput("distPlot1", height = "400px"),
            shiny::plotOutput("distPlot2", height = "400px")
          ),
          shiny::tabPanel(
            "Alignment Results",
            shiny::plotOutput("alignPlot", height = "800px")
          ),
          shiny::tabPanel("Results Summary", shiny::verbatimTextOutput("summary")),
          shiny::tabPanel(
            "Optimization Progress",
            shiny::conditionalPanel(
              condition = "input.optimize",
              shiny::plotOutput("optProgress"),
              shiny::verbatimTextOutput("optSummary")
            )
          )
        )
      )
    )
  )

  # Server logic
  server <- function(input, output, session) {
    # Reactive values to store data and results
    rv <- shiny::reactiveValues(
      df1 = NULL,
      df2 = NULL,
      ms1 = NULL,
      ms2 = NULL,
      results = NULL,
      opt_history = NULL,
      ms_objects_ready = FALSE
    )

    # Update column choices when files are uploaded - separate observers for each dataset
    observe({
      req(input$file1)  # Only require file1
      tryCatch({
        rv$df1 <- utils::read.csv(input$file1$datapath)
        cols1 <- names(rv$df1)

        # Update column selection inputs for dataset 1
        updateSelectInput(session,
                          "id_name1",
                          choices = c("Select Column" = "", cols1))
        updateSelectInput(session,
                          "rt_name1",
                          choices = c("Select Column" = "", cols1))
        updateSelectInput(session,
                          "mz_name1",
                          choices = c("Select Column" = "", cols1))
        updateSelectInput(session,
                          "int_name1",
                          choices = c("Select Column" = "", cols1))
        updateSelectInput(session,
                          "metab_name1",
                          choices = c("Select Column" = "", cols1))

        showNotification("Dataset 1 loaded successfully!", type = "message")
      }, error = function(e) {
        showNotification(paste("Error loading Dataset 1:", e$message), type = "error")
      })
    })

    observe({
      req(input$file2)  # Only require file2
      tryCatch({
        rv$df2 <- utils::read.csv(input$file2$datapath)
        cols2 <- names(rv$df2)

        # Update column selection inputs for dataset 2
        updateSelectInput(session,
                          "id_name2",
                          choices = c("Select Column" = "", cols2))
        updateSelectInput(session,
                          "rt_name2",
                          choices = c("Select Column" = "", cols2))
        updateSelectInput(session,
                          "mz_name2",
                          choices = c("Select Column" = "", cols2))
        updateSelectInput(session,
                          "int_name2",
                          choices = c("Select Column" = "", cols2))
        updateSelectInput(session,
                          "metab_name2",
                          choices = c("Select Column" = "", cols2))

        showNotification("Dataset 2 loaded successfully!", type = "message")
      }, error = function(e) {
        showNotification(paste("Error loading Dataset 2:", e$message), type = "error")
      })
    })

    # Create MS objects when button is clicked
    observeEvent(input$create_ms, {
      # Validate column selections
      if (input$id_name1 == "" || input$rt_name1 == "" ||
          input$mz_name1 == "" || input$int_name1 == "" ||
          input$metab_name1 == "") {
        showNotification("Please select all required columns for Dataset 1", type = "warning")
        return()
      }

      if (input$id_name2 == "" || input$rt_name2 == "" ||
          input$mz_name2 == "" || input$int_name2 == "" ||
          input$metab_name2 == "") {
        showNotification("Please select all required columns for Dataset 2", type = "warning")
        return()
      }

      withProgress(message = 'Creating MS objects...', {
        tryCatch({
          # Create first MS object
          rv$ms1 <- create_ms_obj(
            df = rv$df1,
            name = input$name1,
            id_name = input$id_name1,
            rt_name = input$rt_name1,
            mz_name = input$mz_name1,
            int_name = input$int_name1,
            metab_name = input$metab_name1
          )

          # Create second MS object
          rv$ms2 <- create_ms_obj(
            df = rv$df2,
            name = input$name2,
            id_name = input$id_name2,
            rt_name = input$rt_name2,
            mz_name = input$mz_name2,
            int_name = input$int_name2,
            metab_name = input$metab_name2
          )

          rv$ms_objects_ready <- TRUE
          showNotification("MS objects created successfully!", type = "message")

        }, error = function(e) {
          showNotification(
            paste("Error creating MS objects:", e$message),
            type = "error",
            duration = NULL
          )
          rv$ms_objects_ready <- FALSE
        })
      })
    })

    # Output to control analysis parameters visibility
    output$ms_objects_ready <- reactive({
      return(rv$ms_objects_ready)
    })
    outputOptions(output, "ms_objects_ready", suspendWhenHidden = FALSE)

    # Distribution plots
    output$distPlot1 <- shiny::renderPlot({
      # Add validation message when files are uploaded but MS objects aren't created yet
      validate(
        need(
          input$file1 && input$file2,
          "Please upload both datasets first"
        ),
        need(
          rv$ms_objects_ready,
          "Please click 'Create MS Objects' button to generate plots"
        )
      )

      withProgress(message = 'Generating distribution plot for Dataset 1...', {
        massSight::distribution_plot(rv$ms1)
      })
    })

    output$distPlot2 <- shiny::renderPlot({
      validate(
        need(
          input$file1 && input$file2,
          "Please upload both datasets first"
        ),
        need(
          rv$ms_objects_ready,
          "Please click 'Create MS Objects' button to generate plots"
        )
      )

      withProgress(message = 'Generating distribution plot for Dataset 2...', {
        massSight::distribution_plot(rv$ms2)
      })
    })

    # Data Preview
    output$dataPreview <- DT::renderDT({
      req(rv$df1, rv$df2)

      df1_preview <- head(rv$df1, 5) |>
        dplyr::mutate(Source = input$name1, .before = 1)

      df2_preview <- head(rv$df2, 5) |>
        dplyr::mutate(Source = input$name2, .before = 1)

      combined_preview <- dplyr::bind_rows(df1_preview, df2_preview)

      DT::datatable(
        combined_preview,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'ftip'
        ),
        filter = 'top',
        style = 'bootstrap'
      )
    })

    # Enable/disable run button based on MS objects status
    observe({
      if (rv$ms_objects_ready) {
        shinyjs::enable("run")
      } else {
        shinyjs::disable("run")
      }
    })

    # Run analysis when button is clicked
    observeEvent(input$run, {
      req(rv$ms1, rv$ms2, rv$ms_objects_ready)

      withProgress(message = 'Running analysis...', {
        if (input$optimize) {
          rv$results <- mass_combine(
            ms1 = rv$ms1,
            ms2 = rv$ms2,
            optimize = TRUE,
            n_iter = input$n_iter,
            match_method = input$match_method,
            smooth_method = input$smooth_method,
            minimum_intensity = input$minimum_intensity
          )

          rv$opt_history <- attr(rv$results, "optimization")

        } else {
          rv$results <- mass_combine(
            ms1 = rv$ms1,
            ms2 = rv$ms2,
            optimize = FALSE,
            rt_delta = input$rt_delta,
            mz_delta = input$mz_delta,
            minimum_intensity = input$minimum_intensity,
            iso_method = input$iso_method,
            rt_iso_threshold = input$rt_iso_threshold,
            mz_iso_threshold = input$mz_iso_threshold,
            match_method = input$match_method,
            smooth_method = input$smooth_method
          )
        }
      })
    })

    # Alignment plot output
    output$alignPlot <- shiny::renderPlot({
      req(rv$results)
      final_plots(rv$results)
    })

    # Results summary output
    output$summary <- shiny::renderPrint({
      req(rv$results)
      cat("Alignment Summary:\n")
      cat("Number of matched pairs:", nrow(all_matched(rv$results)), "\n")
      cat("Number of isolated pairs:", nrow(iso_matched(rv$results)), "\n")

      if (!is.null(rv$opt_history)) {
        cat("\nOptimization Results:\n")
        cat("Final score:",
            round(rv$opt_history$final_score, 4),
            "\n")
        cat("\nOptimal parameters:\n")
        for (param in names(rv$opt_history$parameters)) {
          cat(sprintf("  %s: %.4f\n", param, rv$opt_history$parameters[[param]]))
        }
      }
    })

    # Optimization progress plot
    output$optProgress <- shiny::renderPlot({
      req(rv$opt_history)
      history_df <- rv$opt_history$history

      ggplot2::ggplot(history_df, ggplot2::aes(x = seq_len(nrow(history_df)), y = y)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::labs(x = "Iteration", y = "Score", title = "Optimization Progress") +
        ggplot2::theme_minimal()
    })
  }

  # Run the app
  shiny::shinyApp(ui = ui, server = server)
}
