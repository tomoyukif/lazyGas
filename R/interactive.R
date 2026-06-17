################################################################################
#' Generate interactive summary
#'
#' @param object QTLscan object
#' @param pehno Phenotype names to be drawn
#' @param out_fn Prefix of output file
#'
#' @import reactable
#' @import htmltools
#' @import plotly
#' @description
#' The `what` argument accept a character vector of one or more of the following
#' strings:
#' \itemize{
#'  \item{"scan"}{Draw an interactive Manhattan plot}
#'  \item{"scan_png"}{Draw an static Manhattan plot}
#'  \item{"peakcall"}{Draw an interactive peackcall plot}
#'  \item{"recalc"}{Draw an interactive recalculated peackll plot}
#'  \item{"groups"}{Draw an interactive peack groping list}
#'  \item{"preakcall_haplo"}{Draw an interactive haplotype-wise phenotype distribution plot for peaks before the recalculation}
#'  \item{"recalc_haplo"}{Draw an interactive haplotype-wise phenotype distribution plot for peaks after the recalculation}
#'  \item{"candidate"}{Draw an interactive candidate list}
#' }
#'
#' @seealso [makeInteractiveDashboard()]
#'
#' @export
#'
makeInteractiveSummary <- function(object, pheno,
                                   what = c("scan_png", "peakcall", "recalc", "groups", "recalc_haplo", "candidate"),
                                   out_fn){
  pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
  tag_list <- .interactive_summary_tags(
    object = object,
    pheno = pheno_name,
    what = what,
    candidate_gene_links = FALSE
  )
  save_html(tag_list, file = out_fn)
  invisible(out_fn)
}

.interactive_summary_tags <- function(object,
                                      pheno,
                                      what,
                                      candidate_gene_links = FALSE,
                                      genes_for_links = NULL) {
  what <- match.arg(
    arg = what,
    choices = c("scan", "scan_png", "peakcall", "recalc", "groups",
                "peakcall_haplo", "recalc_haplo", "candidate"),
    several.ok = TRUE
  )

  tag_list <- tagList(div(h1(pheno), style = "text-align:center"))

  plot_pheno <- plotPheno(object = object, pheno = pheno, xlab = pheno, boxplot = FALSE)
  plot_pheno <- ggplotly(plot_pheno)
  tag_list <- tagList(
    tag_list,
    div(h1("Phenotype"), style = "text-align:center"),
    div(plot_pheno, style = "margin:auto;width:80vw;")
  )

  if ("scan" %in% what) {
    plot_man <- plotManhattan(object = object, pheno = pheno)
    plot_man <- ggplotly(plot_man)
    tag_list <- tagList(
      tag_list,
      div(h1("Manhattan plot"), style = "text-align:center"),
      div(plot_man, style = "margin:auto;width:80vw;")
    )
  }

  if ("scan_png" %in% what) {
    plot_man <- plotManhattan(object = object, pheno = pheno)
    plot_man <- .interactive_ggplot_img_tag(
      plot_man,
      alt = paste0("Manhattan plot for ", pheno),
      width = 800, height = 600
    )
    tag_list <- tagList(
      tag_list,
      div(h1("Manhattan plot"), style = "text-align:center"),
      div(plot_man, style = "text-align:center;margin:auto;width:80vw;")
    )
  }

  if ("peakcall" %in% what) {
    peak1 <- lazyData(object = object, dataset = "peakcall", pheno = pheno)
    if (!is.null(peak1)) {
      plot_peak1 <- plotPeaks(object = object, pheno = pheno, recalc = FALSE)
      tag_list <- tagList(
        tag_list,
        div(h1("Peakcall plot"), style = "text-align:center"),
        div(ggplotly(plot_peak1), style = "margin:auto;width:80vw;")
      )
    }
  }

  if ("recalc" %in% what) {
    peak2 <- lazyData(object = object, dataset = "recalc", pheno = pheno)
    if (!is.null(peak2)) {
      plot_peak2 <- plotPeaks(object = object, pheno = pheno, recalc = TRUE)
      tag_list <- tagList(
        tag_list,
        div(h1("Recalculated peakcall plot"), style = "text-align:center"),
        div(ggplotly(plot_peak2), style = "margin:auto;width:80vw;")
      )
    }
  }

  if ("groups" %in% what) {
    groups <- lazyData(object = object, dataset = "groups", pheno = pheno)
    if (!is.null(groups)) {
      table <- reactable(
        data = groups,
        sortable = TRUE,
        resizable = TRUE,
        filterable = TRUE,
        searchable = TRUE,
        showPageSizeOptions = TRUE,
        wrap = TRUE,
        striped = TRUE
      )
      tag_list <- tagList(
        tag_list,
        div(h1("Peak grouping list"), style = "text-align:center"),
        div(table, style = "margin:auto;width:90vw;")
      )
    }
  }

  if ("peakcall_haplo" %in% what) {
    plot_hap <- haploPlot(object = object, pheno = pheno, recalc = FALSE)
    if (!is.null(plot_hap)) {
      tag_list <- tagList(
        tag_list,
        div(h1("Haplotype plot (all peaks)"), style = "text-align:center")
      )
      for (i in seq_along(plot_hap)) {
        tag_list <- tagList(
          tag_list,
          div(ggplotly(plot_hap[[i]]), style = "margin:auto;width:60vw;")
        )
      }
    }
  }

  if ("recalc_haplo" %in% what) {
    plot_hap <- haploPlot(object = object, pheno = pheno, recalc = TRUE)
    if (!is.null(plot_hap)) {
      tag_list <- tagList(
        tag_list,
        div(h1("Haplotype plot (recalculated peaks)"), style = "text-align:center")
      )
      for (i in seq_along(plot_hap)) {
        tag_list <- tagList(
          tag_list,
          div(ggplotly(plot_hap[[i]]), style = "margin:auto;width:60vw;")
        )
      }
    }
  }

  if ("candidate" %in% what) {
    candidate <- lazyData(object = object, dataset = "candidate", pheno = pheno)
    if (length(candidate) == 0L) {
      candidate <- NULL
    }
    if (!is.null(candidate)) {
      col_names <- colnames(candidate)
      col_def <- vector("list", length(col_names))
      names(col_def) <- col_names
      for (i in seq_along(col_names)) {
        min_width <- max(nchar(col_names[i]) * 16, 100)
        col_def[[i]] <- colDef(minWidth = min_width)
      }
      if (candidate_gene_links && "Gene_ID" %in% col_names) {
        col_def[["Gene_ID"]] <- colDef(
          name = "Gene_ID",
          minWidth = 120,
          html = TRUE,
          cell = function(value) {
            if (is.na(value) || !nzchar(value)) {
              return("")
            }
            safe <- .variant_viewer_safe_id(value)
            sprintf(
              '<a href="#" onclick="lazyGasShowGene(\'%s\'); return false;">%s</a>',
              safe,
              htmltools::htmlEscape(value)
            )
          }
        )
      }
      table <- reactable(
        data = candidate,
        columns = col_def,
        sortable = TRUE,
        resizable = TRUE,
        filterable = TRUE,
        searchable = TRUE,
        showPageSizeOptions = TRUE,
        wrap = TRUE,
        striped = TRUE
      )
      hint <- if (candidate_gene_links) {
        p("Click Gene_ID to open the variant viewer below.",
          style = "text-align:center;color:#555;margin-bottom:0.5em;")
      } else {
        NULL
      }
      tag_list <- tagList(
        tag_list,
        div(h1("Candidate list"), style = "text-align:center"),
        hint,
        div(table, style = "margin:auto;width:90vw;")
      )
    }
  }

  tag_list
}

.interactive_ggplot_img_tag <- function(p, alt = "", width = 800, height = 600) {
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
  b64 <- if (requireNamespace("base64enc", quietly = TRUE)) {
    base64enc::base64encode(png_raw)
  } else {
    stop(
      "Package 'base64enc' is required for scan_png plots. ",
      "Install it or use what = \"scan\" for an interactive Manhattan plot.",
      call. = FALSE
    )
  }
  htmltools::tags$img(
    src = paste0("data:image/png;base64,", b64),
    alt = alt,
    width = width,
    height = height,
    style = "max-width:100%;height:auto;"
  )
}
