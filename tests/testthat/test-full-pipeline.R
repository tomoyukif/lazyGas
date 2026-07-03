test_that("full pipeline: build, phenotype, and genotype helpers", {
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")

  gds_fn <- .copy_sample_gds(.skip_without_sample_gds())
  on.exit(unlink(gds_fn), add = TRUE)

  lg <- lazyGas::buildLazyGas(
    gds_fn = gds_fn,
    load_filter = TRUE,
    overwrite = TRUE,
    lazygas_store = "parquet"
  )
  on.exit(.close_lg(lg), add = TRUE)

  expect_s4_class(lg, "LazyGas")

  demo <- .skip_without_demo_extdata()
  pheno <- read.csv(demo$pheno_fn)
  lg <- lazyGas::assignPheno(object = lg, pheno = pheno, rename = "test_trait")

  pheno_obj <- lazyGas::getPheno(object = lg)
  expect_equal(pheno_obj$pheno_names, "test_trait")
  expect_equal(nrow(pheno_obj$pheno), GBScleanR::nsam(lg))

  g_dos <- lazyGas::getGenoPerMarker(object = lg, geno_format = "dosage")
  expect_true(is.numeric(g_dos))
  expect_equal(length(g_dos), GBScleanR::nsam(lg))

  g_hap <- lazyGas::getGenoPerMarker(object = lg, geno_format = "haplotype")
  expect_true(is.matrix(g_hap))

  conv_fun <- lazyGas::makeConvFun(geno_format = "dosage", n_levels = 3L)
  expect_true(is.function(conv_fun))

  p_pheno <- suppressMessages(lazyGas::plotPheno(object = lg, pheno = "test_trait"))
  expect_s3_class(p_pheno, "ggplot")
})

test_that("full pipeline: scan, assignPvalues, and manhattan plot", {
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")

  gds_fn <- .copy_sample_gds(.skip_without_sample_gds())
  on.exit(unlink(gds_fn), add = TRUE)

  lg <- lazyGas::buildLazyGas(gds_fn = gds_fn, load_filter = TRUE, overwrite = TRUE)
  on.exit(.close_lg(lg), add = TRUE)

  demo <- .skip_without_demo_extdata()
  pheno <- read.csv(demo$pheno_fn)
  lg <- lazyGas::assignPheno(object = lg, pheno = pheno, rename = "test_trait")

  conv_fun <- lazyGas::makeConvFun(geno_format = "dosage", n_levels = 3L)
  suppressMessages(
    lazyGas::scanAssoc(
      object = lg,
      formula = "add + dom",
      conv_fun = conv_fun,
      geno_format = "dosage"
    )
  )

  scan_dat <- lazyGas::lazyData(object = lg, dataset = "scan", pheno = "test_trait")
  expect_equal(nrow(scan_dat), GBScleanR::nmar(lg))
  expect_true(all(c("FDR", "negLog10P", "P.model", "P.add") %in% names(scan_dat)))

  lg <- lazyGas::assignPvalues(
    object = lg,
    pheno_name = "test_trait",
    p_values = scan_dat$P.model,
    geno_format = "dosage",
    conv_fun = conv_fun,
    formula = "add + dom"
  )
  scan2 <- lazyGas::lazyData(object = lg, dataset = "scan", pheno = "test_trait")
  expect_equal(scan2$P.model, scan_dat$P.model)

  p_man <- lazyGas::plotManhattan(object = lg, pheno = "test_trait")
  expect_s3_class(p_man, "ggplot")
})

test_that("full pipeline: peak calling, recalc, plots, and haploPlot", {
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")
  skip_if_not_installed("rtracklayer")

  ctx <- .run_sample_pipeline(pheno_name = "test_trait")
  lg <- ctx$lg
  pheno_name <- ctx$pheno_name
  on.exit(.close_lg(lg), add = TRUE)

  peakcall <- lazyGas::lazyData(object = lg, dataset = "peakcall", pheno = pheno_name)
  expect_true(is.data.frame(peakcall))
  expect_gt(nrow(peakcall), 0L)
  expect_true("peak_ID" %in% names(peakcall))

  recalc <- lazyGas::lazyData(object = lg, dataset = "recalc", pheno = pheno_name)
  expect_true(is.data.frame(recalc))
  expect_gt(nrow(recalc), 0L)

  groups <- lazyGas::lazyData(object = lg, dataset = "groups", pheno = pheno_name)
  expect_true(is.data.frame(groups))

  p_peak <- suppressWarnings(lazyGas::plotPeaks(object = lg, pheno = pheno_name, recalc = FALSE))
  expect_s3_class(p_peak, "ggplot")

  p_recalc <- suppressWarnings(lazyGas::plotPeaks(object = lg, pheno = pheno_name, recalc = TRUE))
  expect_s3_class(p_recalc, "ggplot")

  hap_peak <- lazyGas::haploPlot(object = lg, pheno = pheno_name, recalc = FALSE)
  expect_true(is.list(hap_peak))
  expect_gt(length(hap_peak), 0L)
  expect_s3_class(hap_peak[[1L]], "ggplot")
  scan_dat <- lazyGas::lazyData(object = lg, dataset = "scan", pheno = pheno_name)
  expect_true("Coef.add" %in% names(scan_dat))
  hap_title <- hap_peak[[1L]]$labels$title
  expect_true(grepl("Coef.add", hap_title, fixed = TRUE))
  expect_match(hap_title, "alt allele (increases|decreases) phenotype")

  hap_recalc <- lazyGas::haploPlot(object = lg, pheno = pheno_name, recalc = TRUE)
  expect_true(is.list(hap_recalc))
  expect_gt(length(hap_recalc), 0L)
  expect_match(hap_recalc[[1L]]$labels$title, "Coef.add")
})

test_that("full pipeline: listCandidate and searchCandidateGenes", {
  skip_if_not_installed("rtracklayer")

  ctx <- .run_sample_pipeline(pheno_name = "test_trait")
  lg <- ctx$lg
  pheno_name <- ctx$pheno_name
  on.exit(.close_lg(lg), add = TRUE)

  candidate <- lazyGas::lazyData(object = lg, dataset = "candidate", pheno = pheno_name)
  expect_true(is.data.frame(candidate))
  expect_gt(nrow(candidate), 0L)
  expect_true("Gene_ID" %in% names(candidate))
  expect_true(all(nzchar(candidate$Gene_ID)))

  snpeff_tab <- lazyGas::lazyData(object = lg, dataset = "snpeff", pheno = pheno_name)
  expect_true(is.data.frame(snpeff_tab))
  expect_gt(nrow(snpeff_tab), 0L)
  expect_true("Gene_ID" %in% names(snpeff_tab))

  ranked <- lazyGas::searchCandidateGenes(
    candidate = candidate,
    query = "demo",
    mode = "keyword",
    keyword_match = "any"
  )
  expect_true(is.data.frame(ranked))
  expect_gt(nrow(ranked), 0L)
})

test_that("full pipeline: variant viewer", {
  skip_if_not_installed("rtracklayer")

  ctx <- .run_sample_pipeline(pheno_name = "test_trait")
  lg <- ctx$lg
  pheno_name <- ctx$pheno_name
  on.exit(.close_lg(lg), add = TRUE)
  gene_id <- lazyGas::lazyData(lg, "candidate", pheno_name)$Gene_ID[[1L]]

  vdata <- lazyGas::getVariantViewerData(
    object = lg,
    gene_id = gene_id,
    gff = ctx$gff,
    pheno = pheno_name,
    snpeff = ctx$snpeff,
    ann = ctx$ann,
    recalc = TRUE
  )
  expect_type(vdata, "list")
  expect_true(inherits(vdata$ggplot, "ggplot"))
  expect_gt(nrow(vdata$plot_df), 0L)
  expect_gt(nrow(vdata$geno_df), 0L)

  p <- lazyGas::plotVariantViewer(vdata)
  expect_s3_class(p, "ggplot")

  p2 <- lazyGas::plotVariantViewer(
    object = lg,
    gene_id = gene_id,
    gff = ctx$gff,
    pheno = pheno_name,
    snpeff = ctx$snpeff,
    ann = ctx$ann,
    recalc = TRUE
  )
  expect_s3_class(p2, "ggplot")
})

test_that("full pipeline: interactive HTML exports", {
  skip_if_not_installed("rtracklayer")

  ctx <- .run_sample_pipeline(pheno_name = "test_trait")
  lg <- ctx$lg
  pheno_name <- ctx$pheno_name
  on.exit(.close_lg(lg), add = TRUE)
  out_dir <- tempfile("lazygas_html_")
  dir.create(out_dir)

  summary_fn <- file.path(out_dir, "summary.html")
  lazyGas::makeInteractiveSummary(
    object = lg,
    pheno = pheno_name,
    out_fn = summary_fn,
    what = c("scan", "peakcall", "recalc", "candidate")
  )
  expect_true(file.exists(summary_fn))
  expect_gt(file.info(summary_fn)$size, 1000L)

  dash_fn <- file.path(out_dir, "dashboard.html")
  what_dash <- c("scan_png", "peakcall", "recalc", "candidate")
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    what_dash <- setdiff(what_dash, "scan_png")
  }
  lazyGas::makeInteractiveDashboard(
    object = lg,
    pheno = pheno_name,
    gff = ctx$gff,
    out_fn = dash_fn,
    snpeff = ctx$snpeff,
    ann = ctx$ann,
    recalc = TRUE,
    what = what_dash
  )
  expect_true(file.exists(dash_fn))
  expect_gt(file.info(dash_fn)$size, 1000L)
  expect_true(grepl("lazyGasShowGene", paste(readLines(dash_fn, n = 200L), collapse = "\n")))

  cand_fn <- file.path(out_dir, "candidate.html")
  lazyGas::makeCanditeList(object = lg, pheno = pheno_name, out_fn = cand_fn)
  expect_true(file.exists(cand_fn))
})

test_that("full pipeline: runDashboardDemo", {
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("rtracklayer")

  out_dir <- tempfile("lazygas_demo_")
  out_fn <- lazyGas::runDashboardDemo(out_dir = out_dir, pheno_name = "demo_trait_test")
  expect_true(file.exists(out_fn))
  expect_gt(file.info(out_fn)$size, 1000L)
})

test_that("full pipeline: snpeff2gds and open_snpeff", {
  demo <- .skip_without_demo_extdata()
  out_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(out_gds), add = TRUE)

  path <- lazyGas::snpeff2gds(vcf_fn = demo$vcf_fn, out_fn = out_gds, verbose = FALSE)
  expect_true(file.exists(path))

  se <- lazyGas::open_snpeff(out_gds)
  expect_true(inherits(se, "snpeff_gds"))
  chr <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(se$root, "chromosome"))
  expect_gt(length(chr), 0L)
})

test_that("full pipeline: importLazyGasResults from legacy GDS store", {
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")

  gds_fn <- .copy_sample_gds(.skip_without_sample_gds())
  on.exit(unlink(gds_fn), add = TRUE)
  on.exit(unlink(sub("\\.gds$", ".lazygas", gds_fn)), add = TRUE)

  lg_gds <- lazyGas::buildLazyGas(
    gds_fn = gds_fn,
    load_filter = TRUE,
    overwrite = TRUE,
    lazygas_store = "gds"
  )
  on.exit(.close_lg(lg_gds), add = TRUE)

  demo <- .skip_without_demo_extdata()
  pheno <- read.csv(demo$pheno_fn)
  lg_gds <- lazyGas::assignPheno(object = lg_gds, pheno = pheno, rename = "legacy_trait")
  suppressMessages(
    lazyGas::scanAssoc(object = lg_gds, formula = "add", geno_format = "dosage")
  )
  GBScleanR::closeGDS(lg_gds, verbose = FALSE)
  lg_gds <- NULL

  companion <- sub("\\.gds$", ".lazygas", gds_fn)
  if (dir.exists(companion)) {
    unlink(companion, recursive = TRUE)
  }

  lg_parquet <- lazyGas::buildLazyGas(
    gds_fn = gds_fn,
    load_filter = TRUE,
    overwrite = TRUE,
    lazygas_store = "parquet"
  )
  on.exit(.close_lg(lg_parquet), add = TRUE)

  lg_parquet <- lazyGas::assignPheno(object = lg_parquet, pheno = pheno, rename = "legacy_trait")
  lg_parquet <- lazyGas::importLazyGasResults(lg_parquet)
  scan_dat <- lazyGas::lazyData(lg_parquet, dataset = "scan", pheno = "legacy_trait")
  expect_equal(nrow(scan_dat), GBScleanR::nmar(lg_parquet))
})
