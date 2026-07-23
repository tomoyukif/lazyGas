library(shiny)
library(lazyGas)

.default_gds_path <- function() {
  opt <- getOption("lazygas.explorer.gds_path", "")
  if (nzchar(opt) && file.exists(opt)) {
    return(opt)
  }
  demo <- file.path(getwd(), "demo_output", "explorer", "sample.gds")
  if (file.exists(demo)) {
    return(normalizePath(demo, winslash = "/", mustWork = FALSE))
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

.default_snpeff_gds_path <- function() {
  demo <- file.path(getwd(), "demo_output", "explorer", "demo_snpeff.gds")
  if (file.exists(demo)) {
    return(normalizePath(demo, winslash = "/", mustWork = FALSE))
  }
  ""
}

.ggplotly_safe <- function(p, ...) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for explorer plots.", call. = FALSE)
  }
  plotly::ggplotly(p, ...)
}

.manhattan_png_tag <- function(p, alt = "Manhattan plot", width = 800, height = 500) {
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  ggplot2::ggsave(
    filename = tmp,
    plot = p,
    width = width / 96,
    height = height / 96,
    dpi = 96,
    bg = "white"
  )
  png_raw <- readBin(tmp, "raw", n = file.info(tmp)$size)
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    return(tags$p(em(
      "Install 'base64enc' for static Manhattan plots, or switch to Interactive."
    )))
  }
  b64 <- base64enc::base64encode(png_raw)
  tags$img(
    src = paste0("data:image/png;base64,", b64),
    alt = alt,
    width = width,
    height = height,
    style = "max-width:100%;height:auto;"
  )
}

.missing_dataset_msg <- function(label) {
  tags$p(em(paste0("No ", label, " data for this phenotype. Run the pipeline step first.")))
}

.sidebar_llm_controls <- function() {
  tagList(
    checkboxInput("use_llm", "Use local LLM for query parsing", value = FALSE),
    textInput("llm_model", "LLM model", value = Sys.getenv("LAZYGAS_LLM_MODEL", "gemma4-31b-64k:latest")),
    textInput("llm_url", "LLM URL", value = Sys.getenv("LAZYGAS_LLM_URL", "http://127.0.0.1:11435")),
    hr(),
    h5("Server actions"),
    actionButton(
      "start_ollama_btn",
      "Start Ollama",
      class = "btn-warning btn-sm btn-block",
      style = "margin-bottom: 4px;"
    ),
    helpText("Start the Ollama API server at the URL above."),
    actionButton(
      "pull_ollama_btn",
      "Pull model",
      class = "btn-default btn-sm btn-block",
      style = "margin-bottom: 4px;"
    ),
    helpText("Download the selected model if it is not installed yet."),
    actionButton(
      "stop_ollama_btn",
      "Stop Ollama (lazyGas)",
      class = "btn-default btn-sm btn-block",
      style = "margin-bottom: 4px;"
    ),
    helpText("Stop only Ollama processes started from this R session."),
    actionButton(
      "setup_apptainer_btn",
      "Setup Ollama (Apptainer)",
      class = "btn-primary btn-sm btn-block",
      style = "margin-bottom: 4px;"
    ),
    helpText("Build or configure the Apptainer image for containerized Ollama."),
    actionButton(
      "refresh_ollama_btn",
      "Refresh status",
      class = "btn-default btn-sm btn-block",
      style = "margin-bottom: 4px;"
    ),
    helpText("Reload server and model status in the main panel.")
  )
}

ui <- fluidPage(
  tags$head(
    tags$script(HTML(
      "window.lazyGasLockButtons = function(ids, busy) {
         if (!Array.isArray(ids)) return;
         ids.forEach(function(id) {
           var el = document.getElementById(id);
           if (!el) return;
           el.disabled = !!busy;
           if (busy) el.classList.add('disabled'); else el.classList.remove('disabled');
         });
       };
       Shiny.addCustomMessageHandler('setBusyButtons', function(msg) {
         if (!msg || !Array.isArray(msg.ids)) return;
         window.lazyGasLockButtons(msg.ids, !!msg.busy);
       });
       Shiny.addCustomMessageHandler('armBusyClickLock', function(msg) {
         if (!msg || !Array.isArray(msg.trigger_ids) || !Array.isArray(msg.lock_ids)) return;
         msg.trigger_ids.forEach(function(id) {
           var el = document.getElementById(id);
           if (!el || el.dataset.lazygasBusyBound === '1') return;
           el.dataset.lazygasBusyBound = '1';
           el.addEventListener('click', function() {
             if (!el.disabled) window.lazyGasLockButtons(msg.lock_ids, true);
           }, true);
         });
       });"
    ))
  ),
  titlePanel("lazyGas phenotype explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      conditionalPanel(
        condition = "input.main_tabs == 'GWAS overview'",
        h4("Project"),
        textInput(
          "gds_path",
          "GDS file path",
          value = .default_gds_path(),
          placeholder = "e.g. demo_output/explorer/sample.gds"
        ),
        textInput(
          "companion_path",
          "Companion store (optional)",
          value = "",
          placeholder = "auto-detect {stem}.lazygas next to GDS"
        ),
        actionButton("load_btn", "Load project", class = "btn-primary"),
        helpText(
          "Use the filesystem path to the GDS file, not a file-upload copy. ",
          "The companion folder (sample.lazygas) must sit beside the GDS. ",
          "To run scan → candidate first, use ",
          tags$code("runLazyGasRunner()"),
          "."
        ),
        hr(),
        h4("Phenotype"),
        selectInput("pheno_name", "Trait", choices = character()),
        hr(),
        h4("GWAS plots"),
        radioButtons(
          "manhattan_mode",
          "Manhattan plot",
          choices = c("Static PNG" = "png", "Interactive" = "plotly"),
          selected = "png",
          inline = TRUE
        )
      ),
      conditionalPanel(
        condition = "input.main_tabs == 'Local LLM'",
        h4("Local LLM"),
        helpText(
          "Optional. When enabled, the LLM parses phenotype queries during ranking ",
          "and powers chat explanations."
        ),
        .sidebar_llm_controls()
      ),
      conditionalPanel(
        condition = "input.main_tabs == 'Ranked genes'",
        h4("Candidate ranking"),
        helpText(
          "Uses the phenotype selected on the GWAS overview tab. ",
          "Configure the local LLM on the Local LLM tab if needed."
        ),
        textAreaInput(
          "trait_text",
          "Phenotype description",
          value = "fruit weight at maturity",
          rows = 3
        ),
        textInput("tissues", "Tissues / organs (comma-separated)", value = "fruit"),
        textInput("stage", "Developmental stage", value = "maturation"),
        textInput("conditions", "Conditions (optional)", value = ""),
        numericInput("top_n", "Top N genes", value = 15, min = 1, max = 200),
        actionButton("rank_btn", "Rank candidates", class = "btn-success")
      ),
      conditionalPanel(
        condition = "input.main_tabs == 'Evidence'",
        h4("Evidence"),
        helpText(
          "Rank candidates on the Ranked genes tab, select a row in the table, ",
          "then open this tab to inspect per-source evidence."
        )
      ),
      conditionalPanel(
        condition = "input.main_tabs == 'Locus / variants'",
        h4("Annotation"),
        textInput(
          "gff_path",
          "GFF path (variant viewer)",
          value = .default_gff_path(),
          placeholder = "path to annotation.gff"
        ),
        textInput(
          "snpeff_gds_path",
          "SnpEff GDS (optional)",
          value = .default_snpeff_gds_path(),
          placeholder = "uses stored snpeff when empty"
        ),
        helpText(
          "Select a ranked gene in the Ranked genes tab to populate haplotype ",
          "and variant views."
        )
      ),
      conditionalPanel(
        condition = "input.main_tabs == 'Chat / explanation'",
        h4("Explanation"),
        helpText(
          "Uses LLM settings from the Local LLM tab when enabled. ",
          "Top N genes follows the value set on the Ranked genes tab."
        ),
        selectInput("explain_lang", "Explanation language", choices = c("ja", "en")),
        actionButton("explain_btn", "Generate explanation", class = "btn-info")
      )
    ),
    mainPanel(
      width = 9,
      uiOutput("global_task_status"),
      tabsetPanel(
        id = "main_tabs",
        tabPanel(
          "GWAS overview",
          br(),
          uiOutput("overview_status"),
          h4("Phenotype"),
          plotly::plotlyOutput("pheno_plot", height = "320px"),
          h4("Manhattan"),
          uiOutput("manhattan_panel"),
          h4("Peakcall"),
          uiOutput("peakcall_panel"),
          h4("Recalculated peaks"),
          uiOutput("recalc_panel"),
          h4("Peak groups"),
          uiOutput("groups_panel")
        ),
        tabPanel(
          "Local LLM",
          br(),
          uiOutput("llm_status_panel")
        ),
        tabPanel(
          "Ranked genes",
          br(),
          uiOutput("rank_status"),
          reactable::reactableOutput("rank_table"),
          br(),
          verbatimTextOutput("query_summary")
        ),
        tabPanel(
          "Evidence",
          br(),
          uiOutput("evidence_panel")
        ),
        tabPanel(
          "Locus / variants",
          br(),
          uiOutput("locus_status"),
          h4("Haplotype (selected peak)"),
          uiOutput("haplo_panel"),
          h4("Variant viewer"),
          uiOutput("variant_panel")
        ),
        tabPanel(
          "Chat / explanation",
          br(),
          uiOutput("chat_status"),
          div(
            style = "max-height: 420px; overflow-y: auto; border: 1px solid #ddd; padding: 8px;",
            uiOutput("chat_history")
          ),
          br(),
          textAreaInput("chat_input", "Follow-up question", rows = 2),
          actionButton("chat_btn", "Ask", class = "btn-default")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  session$userData$lg <- NULL
  lg_obj <- reactiveVal(NULL)
  ranked <- reactiveVal(NULL)
  query_obj <- reactiveVal(NULL)
  chat_msgs <- reactiveVal(list())
  gff_cache <- reactiveVal(NULL)
  gff_cache_path <- reactiveVal("")
  ranking_busy <- reactiveVal(FALSE)
  chat_busy <- reactiveVal(NULL)
  task_active <- reactiveVal(FALSE)
  task_label <- reactiveVal("")
  task_detail <- reactiveVal("")
  task_value <- reactiveVal(0)
  ollama_status_text <- reactiveVal(
    describeOllamaStatus(
      base_url = Sys.getenv("LAZYGAS_LLM_URL", "http://127.0.0.1:11435"),
      model = Sys.getenv("LAZYGAS_LLM_MODEL", "gemma4-31b-64k:latest")
    )
  )

  refresh_ollama_status <- function() {
    ollama_status_text(
      describeOllamaStatus(base_url = input$llm_url, model = input$llm_model)
    )
  }

  button_ids_to_lock <- c(
    "load_btn", "start_ollama_btn", "pull_ollama_btn", "stop_ollama_btn",
    "setup_apptainer_btn", "refresh_ollama_btn", "rank_btn",
    "explain_btn", "chat_btn"
  )

  start_task <- function(label) {
    task_label(label)
    task_detail("")
    task_value(0)
    task_active(TRUE)
  }

  update_task <- function(value = NULL, detail = NULL) {
    if (!is.null(value)) {
      task_value(max(0, min(1, value)))
    }
    if (!is.null(detail)) {
      task_detail(detail)
    }
  }

  finish_task <- function() {
    task_active(FALSE)
    task_label("")
    task_detail("")
    task_value(0)
  }

  observe({
    busy <- isTRUE(task_active())
    session$sendCustomMessage(
      "setBusyButtons",
      list(ids = button_ids_to_lock, busy = busy)
    )
  })

  observe({
    session$sendCustomMessage(
      "armBusyClickLock",
      list(
        trigger_ids = c("rank_btn", "explain_btn", "chat_btn"),
        lock_ids = button_ids_to_lock
      )
    )
  })

  output$global_task_status <- renderUI({
    if (!isTRUE(task_active())) {
      return(NULL)
    }
    pct <- round(task_value() * 100)
    tags$div(
      class = "alert alert-info",
      style = "margin-bottom: 12px;",
      role = "status",
      tags$div(
        style = "display:flex;justify-content:space-between;align-items:center;gap:8px;",
        tags$strong(task_label())
      ),
      tags$div(
        class = "progress",
        style = "margin-top:8px;margin-bottom:8px;",
        tags$div(
          class = "progress-bar progress-bar-info progress-bar-striped active",
          role = "progressbar",
          style = paste0("width:", pct, "%;min-width:2em;"),
          paste0(pct, "%")
        )
      ),
      if (nzchar(task_detail())) tags$p(style = "margin:0;", task_detail())
    )
  })

  output$llm_status_panel <- renderUI({
    input$refresh_ollama_btn
    input$start_ollama_btn
    input$pull_ollama_btn
    input$stop_ollama_btn
    input$setup_apptainer_btn
    input$main_tabs
    input$use_llm
    input$llm_model
    input$llm_url

    running <- tryCatch(
      llmHealthCheck(base_url = input$llm_url, timeout = 2),
      error = function(e) FALSE
    )
    models <- if (isTRUE(running)) {
      tryCatch(
        listOllamaModels(base_url = input$llm_url),
        error = function(e) character()
      )
    } else {
      character()
    }
    model_avail <- isTRUE(running) &&
      nzchar(input$llm_model) &&
      input$llm_model %in% models

    tags$div(
      tags$h4("Runtime"),
      tags$table(
        class = "table table-condensed",
        style = "max-width: 640px;",
        tags$tr(tags$td(tags$strong("Base URL")), tags$td(tags$code(input$llm_url))),
        tags$tr(
          tags$td(tags$strong("Server")),
          tags$td(if (isTRUE(running)) "Running" else "Not reachable")
        ),
        tags$tr(
          tags$td(tags$strong("Use LLM")),
          tags$td(
            if (isTRUE(input$use_llm)) {
              "Enabled for ranking and chat"
            } else {
              "Disabled (rule-based fallback)"
            }
          )
        ),
        tags$tr(tags$td(tags$strong("Selected model")), tags$td(tags$code(input$llm_model))),
        tags$tr(
          tags$td(tags$strong("Model status")),
          tags$td(
            if (!isTRUE(running)) {
              "Start Ollama to check models"
            } else if (model_avail) {
              "Installed"
            } else {
              "Not installed — use Pull model"
            }
          )
        )
      ),
      tags$h4("Diagnostics"),
      tags$pre(
        style = "white-space: pre-wrap; background: #f7f7f7; padding: 10px; border-radius: 4px;",
        ollama_status_text()
      ),
      if (length(models)) {
        tagList(
          tags$h4("Installed models"),
          tags$ul(lapply(models, tags$li))
        )
      } else if (isTRUE(running)) {
        tags$p(em("No models reported by the Ollama API."))
      }
    )
  })

  observeEvent(input$main_tabs, {
    if (identical(input$main_tabs, "Local LLM")) {
      refresh_ollama_status()
    }
  }, ignoreInit = TRUE)

  observeEvent(input$refresh_ollama_btn, {
    refresh_ollama_status()
  })

  observeEvent(input$setup_apptainer_btn, {
    tryCatch({
      configureLazyGasOllama(
        base_url = input$llm_url,
        model = input$llm_model,
        build_sif_if_missing = TRUE
      )
      refresh_ollama_status()
      showNotification("Apptainer Ollama configured.", type = "message")
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = NULL)
      refresh_ollama_status()
    })
  })

  observeEvent(input$start_ollama_btn, {
    tryCatch({
      startLocalLLM(
        model = input$llm_model,
        base_url = input$llm_url,
        pull_model = FALSE,
        start_server = TRUE
      )
      refresh_ollama_status()
      showNotification("Ollama server is ready.", type = "message")
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = NULL)
      refresh_ollama_status()
    })
  })

  observeEvent(input$pull_ollama_btn, {
    tryCatch({
      if (!llmHealthCheck(base_url = input$llm_url, timeout = 2)) {
        startOllamaServer(base_url = input$llm_url, wait_seconds = 45)
      }
      pullOllamaModel(model = input$llm_model, base_url = input$llm_url)
      refresh_ollama_status()
      showNotification(paste("Model", input$llm_model, "is ready."), type = "message")
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = NULL)
      refresh_ollama_status()
    })
  })

  observeEvent(input$stop_ollama_btn, {
    tryCatch({
      stopOllamaServer(base_url = input$llm_url, only_started_by_lazygas = TRUE)
      refresh_ollama_status()
      showNotification("Stop request sent.", type = "message")
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error")
      refresh_ollama_status()
    })
  })

  observeEvent(input$load_btn, {
    gds_raw <- trimws(input$gds_path)
    if (!nzchar(gds_raw)) {
      showNotification("Enter the path to a GDS file.", type = "error")
      return()
    }

    paths <- tryCatch(
      resolveLazyGasPaths(
        gds_fn = gds_raw,
        companion_path = if (nzchar(trimws(input$companion_path))) {
          trimws(input$companion_path)
        } else {
          NULL
        }
      ),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = NULL)
        NULL
      }
    )
    req(paths)

    old_lg <- session$userData$lg
    if (!is.null(old_lg)) {
      try(closeGDS(old_lg, verbose = FALSE), silent = TRUE)
      session$userData$lg <- NULL
    }

    lg <- buildLazyGas(
      gds_fn = paths$gds,
      load_filter = TRUE,
      overwrite = FALSE,
      companion_path = paths$companion
    )
    lg <- restorePhenoFromStore(lg, warn_stub = TRUE)
    session$userData$lg <- lg
    lg_obj(lg)
    ranked(NULL)
    query_obj(NULL)
    chat_msgs(list())
    gff_cache(NULL)
    gff_cache_path("")

    phenos <- getPheno(lg)$pheno_names
    if (!length(phenos)) {
      showNotification(
        paste(
          "No phenotypes found in companion store.",
          if (is.null(paths$companion)) {
            "Expected a folder like sample.lazygas next to the GDS."
          } else {
            paste0("Checked: ", paths$companion)
          }
        ),
        type = "error",
        duration = NULL
      )
      return()
    }
    updateSelectInput(session, "pheno_name", choices = phenos, selected = phenos[1L])
    updateTabsetPanel(session, "main_tabs", selected = "GWAS overview")
    showNotification(
      paste0(
        "Loaded ", length(phenos), " phenotype(s).",
        if (!is.null(paths$companion)) paste0(" Store: ", paths$companion)
      ),
      type = "message"
    )
  })

  observeEvent(input$pheno_name, {
    ranked(NULL)
    query_obj(NULL)
    chat_msgs(list())
  }, ignoreInit = TRUE)

  active_pheno <- reactive({
    req(lg_obj(), nzchar(input$pheno_name))
    input$pheno_name
  })

  output$overview_status <- renderUI({
    if (is.null(lg_obj())) {
      return(tags$p(em("Load a project to view GWAS plots.")))
    }
    tags$p(
      "Phenotype: ", tags$code(active_pheno()),
      " — plots use the same helpers as makeInteractiveDashboard()."
    )
  })

  output$pheno_plot <- plotly::renderPlotly({
    lg <- lg_obj()
    pheno <- active_pheno()
    req(lg, pheno)
    p <- tryCatch(
      plotPheno(object = lg, pheno = pheno, xlab = pheno, boxplot = FALSE),
      error = function(e) NULL
    )
    req(p)
    .ggplotly_safe(p)
  })

  output$manhattan_panel <- renderUI({
    lg <- lg_obj()
    pheno <- active_pheno()
    req(lg, pheno)
    scan_df <- tryCatch(
      lazyData(object = lg, dataset = "scan", pheno = pheno),
      error = function(e) NULL
    )
    if (is.null(scan_df)) {
      return(.missing_dataset_msg("scan"))
    }
    p <- tryCatch(
      plotManhattan(object = lg, pheno = pheno),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error")
        NULL
      }
    )
    req(p)
    if (identical(input$manhattan_mode, "plotly")) {
      return(plotly::plotlyOutput("manhattan_plotly", height = "400px"))
    }
    .manhattan_png_tag(p, alt = paste("Manhattan plot for", pheno))
  })

  output$manhattan_plotly <- plotly::renderPlotly({
    lg <- lg_obj()
    pheno <- active_pheno()
    req(lg, pheno, identical(input$manhattan_mode, "plotly"))
    p <- plotManhattan(object = lg, pheno = pheno)
    .ggplotly_safe(p)
  })

  output$peakcall_panel <- renderUI({
    lg <- lg_obj()
    pheno <- active_pheno()
    req(lg, pheno)
    peak1 <- tryCatch(
      lazyData(object = lg, dataset = "peakcall", pheno = pheno),
      error = function(e) NULL
    )
    if (is.null(peak1)) {
      return(.missing_dataset_msg("peakcall"))
    }
    plotly::plotlyOutput("peakcall_plot", height = "360px")
  })

  output$peakcall_plot <- plotly::renderPlotly({
    lg <- lg_obj()
    pheno <- active_pheno()
    req(lg, pheno)
    p <- plotPeaks(object = lg, pheno = pheno, recalc = FALSE)
    .ggplotly_safe(p, tooltip = "text")
  })

  output$recalc_panel <- renderUI({
    lg <- lg_obj()
    pheno <- active_pheno()
    req(lg, pheno)
    peak2 <- tryCatch(
      lazyData(object = lg, dataset = "recalc", pheno = pheno),
      error = function(e) NULL
    )
    if (is.null(peak2)) {
      return(.missing_dataset_msg("recalc"))
    }
    plotly::plotlyOutput("recalc_plot", height = "360px")
  })

  output$recalc_plot <- plotly::renderPlotly({
    lg <- lg_obj()
    pheno <- active_pheno()
    req(lg, pheno)
    p <- plotPeaks(object = lg, pheno = pheno, recalc = TRUE)
    .ggplotly_safe(p, tooltip = "text")
  })

  output$groups_panel <- renderUI({
    lg <- lg_obj()
    pheno <- active_pheno()
    req(lg, pheno)
    groups <- tryCatch(
      lazyData(object = lg, dataset = "groups", pheno = pheno),
      error = function(e) NULL
    )
    if (is.null(groups) || !nrow(groups)) {
      return(.missing_dataset_msg("groups"))
    }
    reactable::reactableOutput("groups_table")
  })

  output$groups_table <- reactable::renderReactable({
    lg <- lg_obj()
    pheno <- active_pheno()
    req(lg, pheno)
    groups <- lazyData(object = lg, dataset = "groups", pheno = pheno)
    req(!is.null(groups), nrow(groups) > 0L)
    reactable::reactable(
      groups,
      sortable = TRUE,
      resizable = TRUE,
      filterable = TRUE,
      searchable = TRUE,
      compact = TRUE,
      defaultPageSize = 10
    )
  })

  observeEvent(input$rank_btn, {
    if (isTRUE(ranking_busy())) {
      return()
    }
    req(lg_obj(), input$pheno_name, nzchar(input$trait_text))
    lg <- lg_obj()
    ranking_busy(TRUE)
    start_task("Ranking candidates")
    on.exit(ranking_busy(FALSE), add = TRUE)
    on.exit(finish_task(), add = TRUE)

    rank_err <- NULL
    shiny::withProgress(message = "Ranking candidates", value = 0, style = "old", {
        tryCatch({
          update_task(0.05, "Loading reference data")
          shiny::incProgress(0.05, detail = "Loading reference data")
          expr_fn <- system.file("extdata", "demo_expression.csv", package = "lazyGas")
          meta_fn <- system.file("extdata", "demo_expression_meta.csv", package = "lazyGas")
          expr_mat <- NULL
          expr_meta <- NULL
          if (nzchar(expr_fn) && file.exists(expr_fn)) {
            expr_df <- read.csv(expr_fn, stringsAsFactors = FALSE, check.names = FALSE)
            rownames(expr_df) <- expr_df$Gene_ID
            expr_mat <- as.matrix(expr_df[, setdiff(names(expr_df), "Gene_ID"), drop = FALSE])
            if (nzchar(meta_fn) && file.exists(meta_fn)) {
              expr_meta <- read.csv(meta_fn, stringsAsFactors = FALSE)
              rownames(expr_meta) <- expr_meta$sample
            }
          }

          ortho_fn <- system.file("extdata", "demo_orthologs.csv", package = "lazyGas")
          ortho <- if (nzchar(ortho_fn) && file.exists(ortho_fn)) {
            read.csv(ortho_fn, stringsAsFactors = FALSE)
          } else {
            NULL
          }

          update_task(0.25, "Parsing phenotype query")
          shiny::incProgress(0.2, detail = "Parsing phenotype query")
          q <- phenotypeQuery(
            text = input$trait_text,
            tissues = input$tissues,
            stage = input$stage,
            conditions = if (nzchar(input$conditions)) input$conditions else NULL,
            use_llm = isTRUE(input$use_llm),
            llm_model = input$llm_model,
            llm_base_url = input$llm_url,
            object = lg,
            save = TRUE
          )
          query_obj(q)

          update_task(0.45, "Scoring candidates")
          shiny::incProgress(0.2, detail = "Scoring candidates")
          rank_df <- rankPhenotypeCandidates(
            object = lg,
            pheno = input$pheno_name,
            query = q,
            sources = c("annotation", "gwas", "expression", "literature", "ortholog"),
            expression_matrix = expr_mat,
            expression_meta = expr_meta,
            ortholog_table = ortho,
            top_n = input$top_n,
            use_cache = TRUE,
            save = TRUE
          )
          ranked(rank_df)
          chat_msgs(list())
          update_task(1, "Complete")
          shiny::incProgress(0.55, detail = "Complete")
          updateTabsetPanel(session, "main_tabs", selected = "Ranked genes")
          showNotification(
            paste0("Ranked ", nrow(rank_df), " candidate gene(s)."),
            type = "message"
          )
        }, error = function(e) {
          rank_err <<- e
        })
    })
    if (!is.null(rank_err)) {
      showNotification(conditionMessage(rank_err), type = "error", duration = NULL)
    }
  })

  output$rank_status <- renderUI({
    if (!isTRUE(ranking_busy())) {
      return(NULL)
    }
    tags$div(
      class = "alert alert-info",
      style = "margin-bottom: 12px;",
      role = "status",
      tags$strong("Ranking in progress…"),
      tags$p(
        style = "margin: 6px 0 0 0;",
        "Parsing the phenotype query and scoring candidates. ",
        "This may take a minute when the LLM or literature search is enabled."
      )
    )
  })

  output$rank_table <- reactable::renderReactable({
    df <- ranked()
    req(df)
    show_cols <- intersect(
      names(df),
      c(
        "Gene_ID", "Name", "peak_ID", "composite_score",
        "score_annotation", "score_gwas", "score_expression",
        "score_literature", "score_ortholog", "negLog10P"
      )
    )
    reactable::reactable(
      df[, show_cols, drop = FALSE],
      selection = "single",
      onClick = "select",
      highlight = TRUE,
      compact = TRUE,
      defaultPageSize = 10
    )
  })

  output$query_summary <- renderPrint({
    q <- query_obj()
    req(q)
    print(q)
  })

  selected_gene <- reactive({
    df <- ranked()
    req(df)
    sel <- reactable::getReactableState("rank_table", "selected")
    if (length(sel)) {
      return(as.character(df$Gene_ID[sel]))
    }
    if (nrow(df)) {
      return(as.character(df$Gene_ID[1L]))
    }
    NULL
  })

  selected_peak_id <- reactive({
    df <- ranked()
    gid <- selected_gene()
    req(df, gid)
    row <- df[df$Gene_ID == gid, , drop = FALSE]
    if (!nrow(row) || !"peak_ID" %in% names(row)) {
      return(NULL)
    }
    row$peak_ID[1L]
  })

  output$evidence_panel <- renderUI({
    df <- ranked()
    gid <- selected_gene()
    req(df, gid)
    row <- df[df$Gene_ID == gid, , drop = FALSE][1L, , drop = FALSE]
    ev <- tryCatch(
      jsonlite::fromJSON(row$evidence_json[1L], simplifyVector = FALSE),
      error = function(e) list()
    )
    tags$div(
      tags$h4("Gene: ", tags$code(gid)),
      lapply(names(ev), function(src) {
        block <- ev[[src]]
        tags$div(
          tags$h5(src, " (score ", round(as.numeric(block$score %||% 0), 3), ")"),
          if (length(block$snippets)) {
            tags$ul(lapply(block$snippets, tags$li))
          } else {
            tags$p(em("No snippets."))
          }
        )
      })
    )
  })

  get_gff <- reactive({
    path <- trimws(input$gff_path %||% "")
    if (!nzchar(path) || !file.exists(path)) {
      return(NULL)
    }
    if (identical(gff_cache_path(), path) && !is.null(gff_cache())) {
      return(gff_cache())
    }
    gff <- tryCatch(
      rtracklayer::import.gff(path),
      error = function(e) {
        showNotification(
          paste("Failed to load GFF:", conditionMessage(e)),
          type = "error",
          duration = NULL
        )
        NULL
      }
    )
    if (!is.null(gff)) {
      gff_cache(gff)
      gff_cache_path(path)
    }
    gff
  })

  output$locus_status <- renderUI({
    if (is.null(lg_obj())) {
      return(tags$p(em("Load a project first.")))
    }
    if (is.null(ranked())) {
      return(tags$p(em(
        "Rank candidates and select a gene in the Ranked genes tab ",
        "to populate haplotype and variant views."
      )))
    }
    gid <- selected_gene()
    tags$p(
      "Selected gene: ", tags$code(gid %||% "none"),
      if (!is.null(selected_peak_id())) {
        tagList(" (peak ", tags$code(as.character(selected_peak_id())), ")")
      }
    )
  })

  haplo_plot_obj <- reactive({
    lg <- lg_obj()
    pheno <- active_pheno()
    peak_id <- selected_peak_id()
    req(lg, pheno)

    use_recalc <- !is.null(tryCatch(
      lazyData(object = lg, dataset = "recalc", pheno = pheno),
      error = function(e) NULL
    ))
    plots <- tryCatch(
      haploPlot(object = lg, pheno = pheno, recalc = use_recalc),
      error = function(e) {
        message(conditionMessage(e))
        NULL
      }
    )
    if (is.null(plots) || !length(plots)) {
      return(NULL)
    }

    peak_df <- tryCatch(
      lazyData(
        object = lg,
        dataset = if (use_recalc) "recalc" else "peakcall",
        pheno = pheno
      ),
      error = function(e) NULL
    )
    if (is.null(peak_df) || is.null(peak_id) || !"peak_ID" %in% names(peak_df)) {
      return(plots[[1L]])
    }
    # Match haploPlot order (tapply over peak_ID, same as .get_peak_info)
    peak_order <- names(tapply(
      peak_df$peak_negLog10P,
      peak_df$peak_ID,
      function(x) x[[1L]]
    ))
    idx <- match(as.character(peak_id), as.character(peak_order))
    if (is.na(idx) || idx < 1L || idx > length(plots)) {
      return(plots[[1L]])
    }
    plots[[idx]]
  })

  output$haplo_panel <- renderUI({
    req(ranked(), selected_gene())
    p <- haplo_plot_obj()
    if (is.null(p)) {
      return(.missing_dataset_msg("haplotype / peak"))
    }
    plotly::plotlyOutput("haplo_plotly", height = "360px")
  })

  output$haplo_plotly <- plotly::renderPlotly({
    p <- haplo_plot_obj()
    req(p)
    .ggplotly_safe(p)
  })

  variant_data <- reactive({
    lg <- lg_obj()
    pheno <- active_pheno()
    gid <- selected_gene()
    gff <- get_gff()
    req(lg, pheno, gid, gff)

    snpeff_obj <- NULL
    snpeff_path <- trimws(input$snpeff_gds_path %||% "")
    stored_snpeff <- tryCatch(
      lazyData(object = lg, dataset = "snpeff", pheno = pheno),
      error = function(e) NULL
    )
    need_snpeff_gds <- is.null(stored_snpeff) || !nrow(stored_snpeff)
    if (need_snpeff_gds && nzchar(snpeff_path) && file.exists(snpeff_path)) {
      snpeff_obj <- tryCatch(
        open_snpeff(snpeff_path),
        error = function(e) {
          showNotification(
            paste("Failed to open SnpEff GDS:", conditionMessage(e)),
            type = "warning"
          )
          NULL
        }
      )
      if (!is.null(snpeff_obj)) {
        on.exit(try(gdsfmt::closefn.gds(snpeff_obj), silent = TRUE), add = TRUE)
      }
    }

    tryCatch(
      getVariantViewerData(
        object = lg,
        gene_id = gid,
        gff = gff,
        pheno = pheno,
        snpeff = snpeff_obj,
        peak_id = selected_peak_id(),
        recalc = !is.null(tryCatch(
          lazyData(object = lg, dataset = "recalc", pheno = pheno),
          error = function(e) NULL
        ))
      ),
      error = function(e) {
        structure(list(error = conditionMessage(e)), class = "lazygas_variant_error")
      }
    )
  })

  output$variant_panel <- renderUI({
    if (is.null(ranked()) || is.null(selected_gene())) {
      return(NULL)
    }
    if (is.null(get_gff())) {
      return(tags$p(em(
        "Set a valid GFF path in the sidebar to enable the variant viewer."
      )))
    }
    vd <- variant_data()
    if (inherits(vd, "lazygas_variant_error")) {
      return(tags$p(em(vd$error)))
    }
    tagList(
      plotly::plotlyOutput("variant_plot", height = "420px"),
      br(),
      reactable::reactableOutput("variant_table")
    )
  })

  output$variant_plot <- plotly::renderPlotly({
    req(ranked(), selected_gene(), get_gff())
    vd <- variant_data()
    req(!is.null(vd), !inherits(vd, "lazygas_variant_error"), !is.null(vd$ggplot))
    .ggplotly_safe(vd$ggplot, tooltip = "text")
  })

  output$variant_table <- reactable::renderReactable({
    req(ranked(), selected_gene(), get_gff())
    vd <- variant_data()
    req(!inherits(vd, "lazygas_variant_error"), !is.null(vd$geno_df), nrow(vd$geno_df) > 0L)
    reactable::reactable(
      vd$geno_df,
      sortable = TRUE,
      filterable = TRUE,
      searchable = TRUE,
      compact = TRUE,
      defaultPageSize = 10
    )
  })

  observeEvent(input$explain_btn, {
    if (!is.null(chat_busy()) || isTRUE(ranking_busy())) {
      return()
    }
    df <- ranked()
    q <- query_obj()
    req(df, q)

    chat_busy("Generating explanation…")
    start_task("Generating explanation")
    on.exit(chat_busy(NULL), add = TRUE)
    on.exit(finish_task(), add = TRUE)

    explain_err <- NULL
    shiny::withProgress(message = "Generating explanation", value = 0, style = "old", {
        tryCatch({
          update_task(0.2, "Collecting ranked candidates")
          shiny::incProgress(0.2, detail = "Collecting ranked candidates")
          txt <- explainPhenotypeCandidates(
            rank_result = df,
            query = q,
            object = lg_obj(),
            top_n = min(input$top_n, nrow(df)),
            use_llm = isTRUE(input$use_llm),
            language = input$explain_lang,
            model = input$llm_model,
            base_url = input$llm_url
          )
          update_task(0.85, "Formatting response")
          shiny::incProgress(0.65, detail = "Formatting response")
          msgs <- chat_msgs()
          msgs[[length(msgs) + 1L]] <- list(role = "assistant", content = txt)
          chat_msgs(msgs)
          update_task(1, "Complete")
          shiny::incProgress(0.15, detail = "Complete")
          updateTabsetPanel(session, "main_tabs", selected = "Chat / explanation")
        }, error = function(e) {
          explain_err <<- e
        })
    })
    if (!is.null(explain_err)) {
      showNotification(conditionMessage(explain_err), type = "error", duration = NULL)
    }
  })

  observeEvent(input$chat_btn, {
    if (!is.null(chat_busy()) || isTRUE(ranking_busy())) {
      return()
    }
    req(nzchar(input$chat_input))
    df <- ranked()
    q <- query_obj()
    req(df, q)
    user_q <- input$chat_input
    msgs <- chat_msgs()
    msgs[[length(msgs) + 1L]] <- list(role = "user", content = user_q)
    chat_msgs(msgs)
    updateTextAreaInput(session, "chat_input", value = "")

    chat_busy("Answering question…")
    start_task("Answering question")
    on.exit(chat_busy(NULL), add = TRUE)
    on.exit(finish_task(), add = TRUE)

    chat_err <- NULL
    shiny::withProgress(message = "Answering question", value = 0, style = "old", {
        tryCatch({
          update_task(0.2, "Preparing context")
          shiny::incProgress(0.2, detail = "Preparing context")
          use_llm <- isTRUE(input$use_llm) &&
            llmHealthCheck(base_url = input$llm_url, timeout = 3)
          update_task(0.45, if (use_llm) "Calling local LLM" else "Building rule-based answer")
          shiny::incProgress(
            0.25,
            detail = if (use_llm) "Calling local LLM" else "Building rule-based answer"
          )
          reply <- if (!use_llm) {
            answerPhenotypeQuestion(
              rank_result = df,
              question = user_q,
              query = q,
              object = lg_obj(),
              top_n = min(input$top_n, nrow(df)),
              use_llm = FALSE,
              language = input$explain_lang,
              chat_history = msgs
            )
          } else {
            answerPhenotypeQuestion(
              rank_result = df,
              question = user_q,
              query = q,
              object = lg_obj(),
              top_n = min(input$top_n, nrow(df)),
              use_llm = TRUE,
              language = input$explain_lang,
              model = input$llm_model,
              base_url = input$llm_url,
              chat_history = msgs
            )
          }
          update_task(0.9, "Updating chat")
          shiny::incProgress(0.45, detail = "Updating chat")
          msgs <- chat_msgs()
          msgs[[length(msgs) + 1L]] <- list(role = "assistant", content = reply)
          chat_msgs(msgs)
          update_task(1, "Complete")
          shiny::incProgress(0.1, detail = "Complete")
        }, error = function(e) {
          chat_err <<- e
          msgs <- chat_msgs()
          msgs[[length(msgs) + 1L]] <- list(
            role = "assistant",
            content = paste("Error:", conditionMessage(e))
          )
          chat_msgs(msgs)
        })
    })
    if (!is.null(chat_err)) {
      showNotification(conditionMessage(chat_err), type = "error", duration = NULL)
    }
  })

  output$chat_status <- renderUI({
    status <- chat_busy()
    if (is.null(status)) {
      return(NULL)
    }
    tags$div(
      class = "alert alert-info",
      style = "margin-bottom: 12px;",
      role = "status",
      tags$strong(status),
      tags$p(
        style = "margin: 6px 0 0 0;",
        if (grepl("Answering", status, fixed = TRUE)) {
          "Composing a reply from ranked candidates and chat history."
        } else {
          "Summarizing top candidates. This may take longer when the LLM is enabled."
        }
      )
    )
  })

  output$chat_history <- renderUI({
    msgs <- chat_msgs()
    if (!length(msgs)) {
      return(tags$p(em("Run ranking, then generate an explanation or ask a question.")))
    }
    tagList(lapply(msgs, function(m) {
      tags$div(
        style = paste0(
          "margin-bottom: 8px; padding: 6px; border-radius: 4px;",
          if (m$role == "user") " background: #eef;" else " background: #f7f7f7;"
        ),
        tags$strong(if (m$role == "user") "You: " else "Assistant: "),
        tags$span(style = "white-space: pre-wrap;", m$content)
      )
    }))
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
