library(shiny)
library(lazyGas)
library(ggplot2)

ui <- fluidPage(
  titlePanel("lazyGas explorer"),
  sidebarLayout(
    sidebarPanel(
      fileInput("gds_file", "GDS file", accept = c(".gds")),
      textInput("pheno_name", "Phenotype name", value = "Demo_trait"),
      actionButton("load_btn", "Load / refresh", class = "btn-primary"),
      hr(),
      selectInput("dataset", "Dataset", choices = c("scan", "peakcall", "recalc", "candidate", "qc")),
      actionButton("plot_btn", "Plot / summarize", class = "btn-success")
    ),
    mainPanel(
      plotOutput("main_plot"),
      tableOutput("summary_table")
    )
  )
)

server <- function(input, output, session) {
  lg_obj <- reactiveVal(NULL)

  observeEvent(input$load_btn, {
    req(input$gds_file)
    lg <- buildLazyGas(gds_fn = input$gds_file$datapath, load_filter = TRUE, overwrite = FALSE)
    lg_obj(lg)
  })

  observeEvent(input$plot_btn, {
    req(lg_obj())
    lg <- lg_obj()
    pn <- input$pheno_name
    if (input$dataset == "qc") {
      output$main_plot <- renderPlot({
        plotQQ(lg, pheno = pn)
      })
      output$summary_table <- renderTable({
        summarizeGWASQC(lg, pheno = pn)
      })
    } else if (input$dataset %in% c("scan", "peakcall", "recalc", "candidate")) {
      if (input$dataset == "scan") {
        output$main_plot <- renderPlot({
          plotManhattan(lg, pheno = pn)
        })
      } else {
        output$main_plot <- renderPlot({
          plotPeaks(lg, pheno = pn, recalc = input$dataset == "recalc")
        })
      }
      output$summary_table <- renderTable({
        head(lazyData(lg, dataset = input$dataset, pheno = pn), 20)
      })
    }
  })

  session$onSessionEnded(function() {
    lg <- lg_obj()
    if (!is.null(lg)) {
      try(closeGDS(lg, verbose = FALSE), silent = TRUE)
    }
  })
}

shinyApp(ui, server)
