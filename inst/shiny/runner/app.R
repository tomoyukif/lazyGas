library(shiny)
library(lazyGas)

.default_gds_path <- function() {
  opt <- getOption("lazygas.runner.gds_path", "")
  if (nzchar(opt) && file.exists(opt)) {
    return(opt)
  }
  demo <- file.path(getwd(), "demo_output", "explorer", "sample.gds")
  if (file.exists(demo)) {
    return(normalizePath(demo, winslash = "/", mustWork = FALSE))
  }
  ""
}

.default_pheno_path <- function() {
  multitrait <- system.file("extdata", "demo_pheno_multitrait.csv", package = "lazyGas")
  if (nzchar(multitrait) && file.exists(multitrait)) {
    return(multitrait)
  }
  pheno <- system.file("extdata", "pheno.csv", package = "lazyGas")
  if (nzchar(pheno) && file.exists(pheno)) {
    return(pheno)
  }
  ""
}

.default_gff_path <- function() {
  demo <- system.file("extdata", "demo_annotation.gff", package = "lazyGas")
  if (nzchar(demo) && file.exists(demo)) {
    return(demo)
  }
  ""
}

.default_ann_path <- function() {
  demo <- system.file("extdata", "demo_ann.csv", package = "lazyGas")
  if (nzchar(demo) && file.exists(demo)) {
    return(demo)
  }
  ""
}

.default_snpeff_path <- function() {
  demo <- file.path(getwd(), "demo_output", "explorer", "demo_snpeff.gds")
  if (file.exists(demo)) {
    return(normalizePath(demo, winslash = "/", mustWork = FALSE))
  }
  ""
}

.parse_rename <- function(text, n_traits) {
  text <- trimws(text %||% "")
  if (!nzchar(text)) {
    return(NULL)
  }
  parts <- trimws(strsplit(text, ",", fixed = TRUE)[[1L]])
  parts <- parts[nzchar(parts)]
  if (length(parts) != n_traits) {
    stop(
      "Trait rename count (", length(parts), ") must match phenotype columns (",
      n_traits, "). Leave blank to keep CSV column names.",
      call. = FALSE
    )
  }
  parts
}

.has_trait_data <- function(lg, dataset, pheno) {
  res <- tryCatch(
    lazyData(object = lg, dataset = dataset, pheno = pheno),
    error = function(e) NULL
  )
  !is.null(res) && (!is.data.frame(res) || nrow(res) > 0L)
}

.step_status <- function(lg, phenos, steps) {
  rows <- lapply(steps, function(step) {
    ok <- if (step == "candidate") {
      length(phenos) > 0L && all(vapply(
        phenos,
        function(pn) .has_trait_data(lg, "candidate", pn),
        logical(1L)
      ))
    } else if (length(phenos)) {
      .has_trait_data(lg, step, phenos[1L])
    } else {
      FALSE
    }
    data.frame(step = step, status = if (ok) "ok" else "missing", stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

ui <- fluidPage(
  titlePanel("lazyGas pipeline runner"),
  helpText(
    "Run scan → peakcall → recalc → candidate from filesystem paths. ",
    "When finished, open ",
    tags$code("runLazyGasExplorer()"),
    " to explore and rank genes."
  ),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("Input files"),
      textInput(
        "gds_path",
        "GDS file path",
        value = .default_gds_path(),
        placeholder = "path/to/sample.gds"
      ),
      textInput(
        "pheno_path",
        "Phenotype CSV",
        value = .default_pheno_path(),
        placeholder = "must contain id or ID column"
      ),
      textInput(
        "trait_rename",
        "Trait names (optional, comma-separated)",
        value = "",
        placeholder = "e.g. fruit_weight,fruit_length"
      ),
      textInput(
        "gff_path",
        "GFF annotation",
        value = .default_gff_path(),
        placeholder = "required for candidate step"
      ),
      textInput(
        "ann_path",
        "Gene annotation CSV (optional)",
        value = .default_ann_path(),
        placeholder = "demo_ann.csv"
      ),
      textInput(
        "snpeff_path",
        "SnpEff GDS (optional)",
        value = .default_snpeff_path(),
        placeholder = "uses stored snpeff when empty"
      ),
      helpText(
        "Paths must point to files on disk (not uploads). ",
        "The companion folder ",
        tags$code("{stem}.lazygas"),
        " is created beside the GDS."
      ),
      hr(),
      h4("Pipeline"),
      checkboxGroupInput(
        "steps",
        NULL,
        choices = c(
          "Association scan" = "scan",
          "Peak calling" = "peakcall",
          "Recalculate peaks" = "recalc",
          "List candidates" = "candidate"
        ),
        selected = c("scan", "peakcall", "recalc", "candidate")
      ),
      checkboxInput("overwrite", "Overwrite existing store results", value = FALSE),
      tags$details(
        tags$summary(style = "cursor: pointer; margin-bottom: 8px;", "Advanced settings"),
        selectInput(
          "formula",
          "Association formula",
          choices = c("add + dom" = "add + dom", "add" = "add"),
          selected = "add + dom"
        ),
        numericInput("n_threads", "Threads (recalc)", value = 1L, min = 1L, max = 64L),
        numericInput("limit_peakcall", "Max peaks per trait", value = 3L, min = 1L, max = 50L)
      ),
      actionButton("run_btn", "Run pipeline", class = "btn-primary btn-block"),
      hr(),
      uiOutput("next_step_panel")
    ),
    mainPanel(
      width = 8,
      h4("Log"),
      div(
        style = "max-height: 420px; overflow-y: auto; background: #f8f8f8; padding: 8px;",
        verbatimTextOutput("log", placeholder = TRUE)
      ),
      hr(),
      h4("Step status"),
      tableOutput("status_table"),
      hr(),
      h4("Assigned phenotypes"),
      verbatimTextOutput("pheno_summary")
    )
  )
)

server <- function(input, output, session) {
  session$userData$lg <- NULL
  log_lines <- reactiveVal(character())
  lg_obj <- reactiveVal(NULL)
  last_gds <- reactiveVal("")
  running <- reactiveVal(FALSE)

  append_log <- function(...) {
    line <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", paste(..., collapse = ""))
    log_lines(c(log_lines(), line))
  }

  output$log <- renderText({
    lines <- log_lines()
    if (!length(lines)) {
      return("Ready. Set paths and click Run pipeline.")
    }
    paste(lines, collapse = "\n")
  })

  output$pheno_summary <- renderPrint({
    lg <- lg_obj()
    if (is.null(lg)) {
      cat("No project loaded yet.\n")
      return(invisible(NULL))
    }
    phenos <- getPheno(lg)$pheno_names
    if (!length(phenos)) {
      cat("No phenotypes assigned.\n")
    } else {
      cat(paste(phenos, collapse = ", "), "\n")
    }
  })

  output$status_table <- renderTable({
    lg <- lg_obj()
    req(lg)
    phenos <- getPheno(lg)$pheno_names
    if (!length(phenos)) {
      return(data.frame(note = "No phenotypes"))
    }
    steps <- c("scan", "peakcall", "recalc", "candidate")
    .step_status(lg, phenos, steps)
  }, rownames = FALSE)

  output$next_step_panel <- renderUI({
    gds <- last_gds()
    if (!nzchar(gds)) {
      return(NULL)
    }
    tags$div(
      class = "well",
      tags$strong("Next: explore results"),
      tags$p(
        "In your R console (new session or after stopping this app):"
      ),
      tags$pre(
        style = "white-space: pre-wrap;",
        paste0("runLazyGasExplorer(gds_path = \"", gds, "\")")
      )
    )
  })

  observeEvent(input$run_btn, {
    if (isTRUE(running())) {
      return()
    }

    gds_raw <- trimws(input$gds_path)
    pheno_raw <- trimws(input$pheno_path)
    gff_raw <- trimws(input$gff_path)
    steps <- input$steps

    if (!nzchar(gds_raw)) {
      showNotification("Enter a GDS file path.", type = "error")
      return()
    }
    if (!file.exists(gds_raw)) {
      showNotification("GDS file not found.", type = "error")
      return()
    }
    if (!nzchar(pheno_raw) || !file.exists(pheno_raw)) {
      showNotification("Phenotype CSV not found.", type = "error")
      return()
    }
    if ("candidate" %in% steps && (!nzchar(gff_raw) || !file.exists(gff_raw))) {
      showNotification("GFF path is required for the candidate step.", type = "error")
      return()
    }
    if (!length(steps)) {
      showNotification("Select at least one pipeline step.", type = "error")
      return()
    }

    running(TRUE)
    log_lines(character())
    last_gds("")

    old_lg <- session$userData$lg
    if (!is.null(old_lg)) {
      try(closeGDS(old_lg, verbose = FALSE), silent = TRUE)
    }
    session$userData$lg <- NULL
    lg_obj(NULL)

    shiny::withProgress(
      message = "Running lazyGas pipeline",
      value = 0,
      {
        err <- NULL
        lg <- NULL
        snpeff <- NULL
        gds_norm <- normalizePath(gds_raw, winslash = "/", mustWork = TRUE)

        tryCatch({
          append_log("buildLazyGas: ", gds_norm)
          incProgress(0.05, detail = "Loading GDS")
          lg <- buildLazyGas(
            gds_fn = gds_norm,
            load_filter = TRUE,
            overwrite = isTRUE(input$overwrite),
            lazygas_store = "parquet"
          )

          append_log("Reading phenotype CSV: ", pheno_raw)
          pheno_df <- read.csv(pheno_raw, stringsAsFactors = FALSE, check.names = FALSE)
          id_col <- grep("^ID$|^id$", names(pheno_df))
          n_traits <- ncol(pheno_df) - length(id_col)
          rename_vec <- .parse_rename(input$trait_rename, n_traits)

          withCallingHandlers(
            {
              append_log("assignPheno")
              incProgress(0.1, detail = "Assigning phenotypes")
              lg <- assignPheno(object = lg, pheno = pheno_df, rename = rename_vec)
            },
            message = function(m) {
              append_log(conditionMessage(m))
              invokeRestart("muffleMessage")
            }
          )

          phenos <- getPheno(lg)$pheno_names
          append_log("Traits: ", paste(phenos, collapse = ", "))

          gff <- NULL
          ann <- NULL
          snpeff <- NULL
          if ("candidate" %in% steps) {
            append_log("Loading GFF: ", gff_raw)
            gff <- rtracklayer::import.gff(gff_raw)
            ann_raw <- trimws(input$ann_path %||% "")
            if (nzchar(ann_raw) && file.exists(ann_raw)) {
              append_log("Loading annotation: ", ann_raw)
              ann <- read.csv(ann_raw, stringsAsFactors = FALSE)
            }
            snpeff_raw <- trimws(input$snpeff_path %||% "")
            if (nzchar(snpeff_raw) && file.exists(snpeff_raw)) {
              append_log("Opening SnpEff GDS: ", snpeff_raw)
              snpeff <- open_snpeff(snpeff_raw)
            }
          }

          conv_fun <- makeConvFun(geno_format = "dosage", n_levels = 3L)
          append_log("runLazyGas: ", paste(steps, collapse = " → "))

          withCallingHandlers(
            {
              incProgress(0.2, detail = "Running pipeline")
              lg <- runLazyGas(
                object = lg,
                steps = steps,
                resume = !isTRUE(input$overwrite),
                gff = gff,
                snpeff = snpeff,
                ann = ann,
                recalc = TRUE,
                formula = input$formula,
                conv_fun = conv_fun,
                geno_format = "dosage",
                limit_peakcall = as.integer(input$limit_peakcall),
                n_threads = as.integer(input$n_threads)
              )
            },
            message = function(m) {
              append_log(conditionMessage(m))
              invokeRestart("muffleMessage")
            }
          )

          hist <- tryCatch(
            lazyData(object = lg, dataset = "pipeline"),
            error = function(e) NULL
          )
          if (!is.null(hist) && nrow(hist)) {
            last <- hist[nrow(hist), , drop = FALSE]
            append_log(
              "History: ran ", last$steps_ran,
              if (nzchar(last$steps_skipped)) paste0(" (skipped: ", last$steps_skipped, ")")
            )
          }

          append_log("Done.")
          incProgress(1, detail = "Complete")
        }, error = function(e) {
          err <<- conditionMessage(e)
          append_log("ERROR: ", err)
          if (!is.null(lg)) {
            try(closeGDS(lg, verbose = FALSE), silent = TRUE)
            lg <<- NULL
            session$userData$lg <- NULL
          }
        }, finally = {
          if (!is.null(snpeff)) {
            try(gdsfmt::closefn.gds(snpeff), silent = TRUE)
          }
        })

        if (!is.null(err)) {
          showNotification(err, type = "error", duration = NULL)
        } else {
          session$userData$lg <- lg
          lg_obj(lg)
          last_gds(gds_norm)
          showNotification("Pipeline finished. Open the explorer to continue.", type = "message")
        }
        running(FALSE)
      }
    )
  })

  session$onSessionEnded(function() {
    lg <- session$userData$lg
    if (!is.null(lg)) {
      try(closeGDS(lg, verbose = FALSE), silent = TRUE)
      session$userData$lg <- NULL
    }
  })
}

`%||%` <- function(x, y) if (is.null(x)) y else x

shinyApp(ui, server)
