#!/usr/bin/env Rscript
# Run from the lazyGas package root:
#   Rscript inst/demo/run_dashboard_demo.R [output_dir]
# Or from R after devtools::load_all():
#   runDashboardDemo()

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else file.path(getwd(), "demo_output")

find_pkg_root <- function() {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(script) == 1L && nzchar(script)) {
    candidate <- normalizePath(file.path(dirname(script), "..", ".."), winslash = "/", mustWork = FALSE)
    if (dir.exists(file.path(candidate, "R"))) {
      return(candidate)
    }
  }
  if (dir.exists(file.path(getwd(), "R"))) {
    return(normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  }
  NULL
}

if (!exists("runDashboardDemo", mode = "function")) {
  pkg_root <- find_pkg_root()
  if (!is.null(pkg_root) && requireNamespace("devtools", quietly = TRUE)) {
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(pkg_root)
    suppressPackageStartupMessages(devtools::load_all(pkg_root, quiet = TRUE, export_all = FALSE))
  }
}

if (!exists("runDashboardDemo", mode = "function")) {
  stop(
    "runDashboardDemo() not found. Run from the package root with devtools::load_all(), ",
    "or install lazyGas with demo extdata files.",
    call. = FALSE
  )
}

runDashboardDemo(out_dir = out_dir)
