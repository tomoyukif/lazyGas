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
#'  \item{"qq"}{Draw a QQ plot}
#'  \item{"qc"}{Draw a GWAS QC summary table}
#'  \item{"cross_trait"}{Draw cross-trait peak clustering (requires multiple phenotypes)}
#'  \item{"fine_mapping"}{Draw fine-mapping interpretation summary (conditional signals and credible-set resolution)}
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
                "peakcall_haplo", "recalc_haplo", "candidate", "qq", "qc", "cross_trait",
                "fine_mapping"),
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
        div(.ggplotly_peak_plot(plot_peak1), style = "margin:auto;width:80vw;")
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
        div(.ggplotly_peak_plot(plot_peak2), style = "margin:auto;width:80vw;")
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

  if ("qq" %in% what) {
    # Genome-wide QQ via ggplotly can be multi-GB; embed a static PNG instead.
    plot_qq <- plotQQ(object = object, pheno = pheno)
    tag_list <- tagList(
      tag_list,
      div(h1("QQ plot"), style = "text-align:center"),
      div(
        .interactive_ggplot_img_tag(
          plot_qq,
          alt = paste0("QQ plot for ", pheno),
          width = 800,
          height = 600
        ),
        style = "text-align:center;margin:auto;width:80vw;"
      )
    )
  }

  if ("qc" %in% what) {
    qc_df <- summarizeGWASQC(object = object, pheno = pheno, store = TRUE)
    if (!is.null(qc_df) && nrow(qc_df) > 0L) {
      tag_list <- tagList(
        tag_list,
        div(h1("GWAS QC summary"), style = "text-align:center"),
        div(
          reactable(
            data = qc_df,
            sortable = TRUE,
            resizable = TRUE,
            striped = TRUE
          ),
          style = "margin:auto;width:80vw;"
        )
      )
    }
  }

  if ("cross_trait" %in% what) {
    pheno_names <- getPheno(object)$pheno_names
    if (length(pheno_names) >= 2L) {
      mt <- tryCatch(
        plotMultiTraitOverview(object = object, recalc = TRUE),
        error = function(e) NULL
      )
      if (!is.null(mt)) {
        tag_list <- tagList(
          tag_list,
          div(h1("Cross-trait peak overview"), style = "text-align:center"),
          div(ggplotly(mt$heatmap), style = "margin:auto;width:80vw;")
        )
        if (nrow(mt$summary) > 0L) {
          tag_list <- tagList(
            tag_list,
            div(h2("Shared peak clusters"), style = "text-align:center"),
            div(
              reactable(data = mt$summary, sortable = TRUE, striped = TRUE),
              style = "margin:auto;width:90vw;"
            )
          )
        }
      }
    }
  }

  if ("fine_mapping" %in% what) {
    fm <- tryCatch(
      summarizeFineMapping(
        object = object,
        pheno = pheno,
        recalc = TRUE,
        run_if_missing = TRUE,
        store = TRUE
      ),
      error = function(e) NULL
    )
    if (!is.null(fm) && nrow(fm) > 0L) {
      fm_col_names <- colnames(fm)
      fm_col_def <- vector("list", length(fm_col_names))
      names(fm_col_def) <- fm_col_names
      for (i in seq_along(fm_col_names)) {
        nm <- fm_col_names[i]
        if (nm %in% c("conditional_message", "credible_set_message")) {
          fm_col_def[[i]] <- colDef(minWidth = 360, wrap = TRUE)
        } else {
          fm_col_def[[i]] <- colDef(minWidth = max(nchar(nm) * 12, 90))
        }
      }
      tag_list <- tagList(
        tag_list,
        div(h1("Fine-mapping interpretation"), style = "text-align:center"),
        p(
          "Credible-set size reflects mapping resolution, not proof of causality. ",
          "Large sets indicate the locus could not be narrowed beyond the LD block.",
          style = "text-align:center;color:#555;max-width:80vw;margin:0.5em auto 1em;"
        ),
        div(
          reactable(
            data = fm,
            columns = fm_col_def,
            sortable = TRUE,
            resizable = TRUE,
            striped = TRUE,
            wrap = TRUE,
            defaultPageSize = 10L
          ),
          style = "margin:auto;width:98vw;overflow-x:auto;"
        )
      )
    }
  }

  if ("candidate" %in% what) {
    candidate <- lazyData(object = object, dataset = "candidate", pheno = pheno)
    if (length(candidate) == 0L) {
      candidate <- NULL
    }
    if (!is.null(candidate)) {
      if (requireNamespace("DT", quietly = TRUE)) {
        table_block <- .interactive_candidate_table_dt(
          candidate = candidate,
          candidate_gene_links = candidate_gene_links
        )
      } else {
        candidate <- .interactive_candidate_display_df(candidate)
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
        table_block <- div(
          reactable(
            data = candidate,
            columns = col_def,
            sortable = TRUE,
            resizable = TRUE,
            filterable = TRUE,
            searchable = TRUE,
            showPageSizeOptions = TRUE,
            defaultPageSize = 20L,
            wrap = FALSE,
            striped = TRUE
          ),
          style = "margin:auto;width:90vw;"
        )
      }
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
        table_block
      )
    }
  }

  tag_list
}

#' Trim/prioritize candidate columns for HTML tables
#' @noRd
.interactive_candidate_display_df <- function(candidate, max_text = 120L) {
  prefer <- c(
    "peak_ID", "Gene_ID", "Gene_chr", "Gene_start", "dist2peak", "negLog10P",
    "HIGH", "MODERATE", "LOW", "MODIFIER",
    "RAPDB_geneSymbol", "Oryzabase_geneSymbol", "CGSNL_geneSymbol",
    "RAP_Note", "MSU_Note", "EGGNog_Description", "MapMan_DESCRIPTION"
  )
  keep <- unique(c(intersect(prefer, names(candidate)),
                   setdiff(names(candidate), prefer)))
  # Drop very wide ontology dumps from the default HTML table
  drop <- intersect(
    keep,
    c("EGGNog_GOs", "EGGNog_KEGG_Pathway", "InterPro_Description",
      "Oryzabase_Trait_Ontology", "Oryzabase_Plant_Ontology", "ID", "gene_id")
  )
  keep <- setdiff(keep, drop)
  out <- as.data.frame(candidate[, keep, drop = FALSE], stringsAsFactors = FALSE)
  for (nm in names(out)) {
    if (is.character(out[[nm]]) || is.factor(out[[nm]])) {
      x <- as.character(out[[nm]])
      long <- !is.na(x) & nchar(x) > max_text
      x[long] <- paste0(substr(x[long], 1L, max_text - 3L), "...")
      out[[nm]] <- x
    }
  }
  out
}

#' DT candidate table with peak filter (falls back to reactable if DT missing)
#' @noRd
.interactive_candidate_table_dt <- function(candidate,
                                            candidate_gene_links = FALSE) {
  display <- .interactive_candidate_display_df(candidate)
  table_id <- paste0(
    "lazygas_cand_",
    substr(
      gsub("[^a-zA-Z0-9]", "", paste(c(names(display), nrow(display)), collapse = "")),
      1L,
      16L
    )
  )
  peak_col0 <- if ("peak_ID" %in% names(display)) {
    which(names(display) == "peak_ID")[1L] - 1L
  } else {
    NA_integer_
  }
  filter_ui <- if (!is.na(peak_col0)) {
    peaks <- sort(unique(display$peak_ID))
    tags$div(
      style = "margin:0.5em auto 1em;max-width:90vw;",
      tags$label(`for` = paste0(table_id, "_peak"), "Filter by peak_ID: "),
      tags$select(
        id = paste0(table_id, "_peak"),
        tags$option(value = "All", "All peaks"),
        lapply(as.character(peaks), function(pid) {
          tags$option(value = pid, pid)
        }),
        onchange = sprintf(
          paste0(
            "var v=document.getElementById('%s_peak').value;",
            "var t=$('#%s table').DataTable();",
            "if(v==='All'){t.column(%d).search('').draw();}",
            "else{t.column(%d).search('^'+v+'$', true, false).draw();}"
          ),
          table_id,
          table_id,
          peak_col0,
          peak_col0
        )
      )
    )
  } else {
    NULL
  }
  if (isTRUE(candidate_gene_links) && "Gene_ID" %in% names(display)) {
    display$Gene_ID <- vapply(as.character(display$Gene_ID), function(value) {
      if (is.na(value) || !nzchar(value)) {
        return("")
      }
      safe <- .variant_viewer_safe_id(value)
      sprintf(
        '<a href="#" onclick="lazyGasShowGene(\'%s\'); return false;">%s</a>',
        safe,
        htmltools::htmlEscape(value)
      )
    }, character(1))
  }
  dt <- DT::datatable(
    display,
    extensions = c("Buttons", "ColReorder"),
    rownames = FALSE,
    escape = if (isTRUE(candidate_gene_links)) {
      setdiff(seq_along(display), which(names(display) == "Gene_ID"))
    } else {
      TRUE
    },
    elementId = table_id,
    options = list(
      dom = "Blfrtip",
      buttons = list("colvis"),
      colReorder = TRUE,
      pageLength = 20L,
      lengthMenu = list(c(10L, 20L, 50L, 100L, -1L), c("10", "20", "50", "100", "All")),
      scrollX = TRUE
    )
  )
  tagList(
    tags$p(
      style = "text-align:center;color:#555;margin:0.25em auto;",
      "Use Colvis to show/hide columns; drag headers to reorder; change page length or peak filter above."
    ),
    filter_ui,
    div(dt, style = "margin:auto;width:90vw;")
  )
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
