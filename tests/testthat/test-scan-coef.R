test_that("scan Coef.* columns use original phenotype units", {
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

  demo <- .skip_without_demo_extdata()
  pheno <- read.csv(demo$pheno_fn)
  pheno_name <- "test_trait"
  lg <- lazyGas::assignPheno(object = lg, pheno = pheno, rename = pheno_name)

  conv_fun <- lazyGas::makeConvFun(geno_format = "dosage", n_levels = 3L)
  lazyGas::scanAssoc(
    object = lg,
    formula = "add + dom",
    conv_fun = conv_fun,
    geno_format = "dosage"
  )

  scan_dat <- lazyGas::lazyData(object = lg, dataset = "scan", pheno = pheno_name)
  expect_true("Coef.add" %in% names(scan_dat))

  valid <- which(is.finite(scan_dat$Coef.add))
  skip_if(length(valid) == 0L, "no finite Coef.add in scan results")

  marker_idx <- valid[which.max(scan_dat$negLog10P[valid])]
  variant_id <- scan_dat$variant_ID[marker_idx]

  pheno_vec <- lazyGas::getPheno(lg)$pheno[[pheno_name]]
  g <- lazyGas::getGenoPerMarker(
    object = lg,
    geno_format = "dosage",
    marker_index = marker_idx
  )
  g <- as.numeric(g)
  g[g == 63] <- NA

  df <- data.frame(
    phe = pheno_vec,
    add = g,
    dom = as.numeric(g == 1)
  )
  df <- df[stats::complete.cases(df), , drop = FALSE]
  skip_if(nrow(df) < 10L, "insufficient complete cases for manual glm")

  fit <- stats::glm(phe ~ add + dom, data = df, family = "gaussian")
  expected <- stats::coef(fit)[c("add", "dom")]

  expect_equal(
    as.numeric(scan_dat[marker_idx, c("Coef.add", "Coef.dom")]),
    as.numeric(expected),
    tolerance = 1e-5
  )
  expect_equal(scan_dat$variant_ID[marker_idx], variant_id)
})
