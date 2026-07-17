#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1L]] else file.path(getwd(), "demo_output", "explorer")
pkg_root <- if (file.exists("DESCRIPTION")) "." else system.file(package = "lazyGas")
if (file.exists(file.path(pkg_root, "DESCRIPTION")) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(pkg_root, quiet = TRUE)
} else if (!requireNamespace("lazyGas", quietly = TRUE)) {
  stop("Install lazyGas or run from the package source tree with devtools.")
}
lazyGas::runExplorerDemo(out_dir = out_dir, overwrite = TRUE, use_llm = FALSE)
