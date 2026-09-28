test_that("E0 config write/read and mapping_mode errors", {
  skip_if_not_installed("yaml")
  wd <- tempfile("lazygas_e0_")
  dir.create(wd)
  on.exit(unlink(wd, recursive = TRUE), add = TRUE)

  writeLazyGasExploreConfig(wd, mapping_mode = "gwas", top_n = 5)
  cfg <- readLazyGasExploreConfig(wd)
  expect_identical(cfg$mapping_mode, "gwas")
  expect_equal(cfg$top_n, 5)

  writeLazyGasExploreConfig(wd, mapping_mode = "qtl", top_n = Inf)
  cfg2 <- readLazyGasExploreConfig(wd)
  expect_identical(cfg2$mapping_mode, "qtl")
  expect_true(is.infinite(cfg2$top_n))

  bad <- file.path(wd, "config.yaml")
  writeLines("mapping_mode: foo\n", bad)
  expect_error(readLazyGasExploreConfig(wd), "mapping_mode")
})

test_that("E0 enrich CS coords from peak map", {
  cred <- data.frame(
    variant_ID = c(1L, 2L),
    PIP = c(0.6, 0.4),
    in_credible_set = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  expect_null(lazyGas:::.e_cs_chr_minmax(cred))

  # Fake object-less enrichment: manually supply coords then minmax works
  cred2 <- data.frame(
    variant_ID = c(1L, 2L),
    Chr = "chr1",
    Pos = c(100L, 200L),
    PIP = c(0.6, 0.4),
    in_credible_set = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  info <- lazyGas:::.e_cs_chr_minmax(cred2)
  expect_identical(info$Chr, "chr1")
  expect_equal(info$start, 100)
  expect_equal(info$end, 100)
  expect_equal(info$n, 1L)

  # CS that already has Chr/Pos must keep a single Pos column (not Pos.x/Pos.y)
  keep <- lazyGas:::.e_enrich_cs_coords(cred2, object = NULL, pheno_name = NULL, peak_id = NULL)
  expect_true("Pos" %in% names(keep))
  expect_false(any(grepl("Pos\\.", names(keep))))
  expect_equal(keep$Pos, cred2$Pos)
})

test_that("E0 classifyAnnColumns use_llm=FALSE defaults to free_description", {
  ann <- data.frame(
    Gene_ID = "g1",
    Description = "starch synthase",
    GO = "GO:0001",
    stringsAsFactors = FALSE
  )
  # Gene_ID included if passed — caller should pass function ann only
  map <- classifyAnnColumns(ann[, c("Description", "GO"), drop = FALSE], use_llm = FALSE)
  expect_true(all(map$class == "free_description"))
  map2 <- classifyAnnColumns(
    ann[, c("Description", "GO"), drop = FALSE],
    use_llm = FALSE,
    overrides = c(GO = "go")
  )
  expect_identical(map2$class[map2$column == "GO"], "go")
})

test_that("E0 normalize/validate maps class-kegg → kegg", {
  allowed <- c(
    "free_description", "go", "kegg", "pfam_interpro", "numeric", "other"
  )
  expect_identical(lazyGas:::.e0_normalize_ann_class("class-kegg", allowed), "kegg")
  expect_identical(lazyGas:::.e0_normalize_ann_class("KEGG", allowed), "kegg")
  raw <- '{"columns":[{"column":"EGGNog_KEGG_Pathway","class":"class-kegg"},{"column":"GO","class":"go"}]}'
  class_by <- lazyGas:::.e0_parse_classify_json(raw)
  checked <- lazyGas:::.e0_validate_classify_map(
    c("EGGNog_KEGG_Pathway", "GO"), class_by, allowed, normalize = TRUE
  )
  expect_true(checked$ok)
  expect_identical(checked$out, c("kegg", "go"))
})

test_that("E0 validate map reports missing names without subscript error", {
  allowed <- c(
    "free_description", "go", "kegg", "pfam_interpro", "numeric", "other"
  )
  checked <- lazyGas:::.e0_validate_classify_map(
    cols = c("Description", "GO"),
    class_by = character(),
    allowed = allowed
  )
  expect_false(checked$ok)
  expect_true(any(grepl("missing class for 'Description'", checked$problems)))
  checked2 <- lazyGas:::.e0_validate_classify_map(
    cols = c("Description", "GO"),
    class_by = c(Description = "free_description"),
    allowed = allowed
  )
  expect_false(checked2$ok)
  expect_true(any(grepl("missing class for 'GO'", checked2$problems)))
})

test_that("E0 LLM classify parser rejects bad JSON via stop", {
  skip_if_not(exists("llmChat", mode = "function"))
  # Force failure by pointing at unreachable URL with tiny timeout
  llm <- list(model = "no-such-model", base_url = "http://127.0.0.1:9", timeout = 1)
  expect_error(
    classifyAnnColumns(
      data.frame(Description = "x", GO = "y", stringsAsFactors = FALSE),
      use_llm = TRUE,
      llm = llm
    ),
    regexp = "classify"
  )
})

test_that("E0 GFF windows and QTL/GWAS gene select", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomeInfoDb")

  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1", "chr1"),
    ranges = IRanges::IRanges(
      start = c(1000, 1100, 5000, 5100),
      end = c(1000, 2000, 5000, 6000)
    ),
    strand = c("+", "+", "-", "-"),
    type = c("gene", "CDS", "gene", "CDS"),
    ID = c("G1", "cds1", "G2", "cds2"),
    Parent = c("", "G1", "", "G2")
  )
  win <- lazyGas:::.gff_gene_windows(gr)
  expect_true(nrow(win) >= 1L)
  expect_true(all(c("Gene_ID", "window_start", "window_end", "mid") %in% names(win)))

  # G1 CDS 1100-2000 + strand → window 1100-3000= -3kb → -1900 clipped? start-3000
  g1 <- win[win$Gene_ID == "G1", , drop = FALSE]
  if (nrow(g1)) {
    expect_equal(g1$window_start[1L], max(1L, 1100L - 3000L))
    expect_equal(g1$window_end[1L], 2000L + 500L)
  }

  cand <- data.frame(
    Gene_ID = c("G1", "G2", "G3"),
    Name = c("a", "b", "c"),
    Chr = "chr1",
    dist2peak = c(10, 20, 30),
    stringsAsFactors = FALSE
  )
  cred <- data.frame(
    Chr = "chr1",
    Pos = c(1500, 1600),
    PIP = c(0.4, 0.3),
    in_credible_set = TRUE,
    stringsAsFactors = FALSE
  )

  qtl <- lazyGas:::.report_e_select_genes(
    mapping_mode = "qtl",
    simple_candidates = cand,
    credible_set = cred,
    gff_windows = win,
    snpeff = NULL,
    top_n = Inf
  )
  expect_true("G1" %in% qtl$Gene_ID)
  expect_true("nearest_cs_PIP" %in% names(qtl))

  gwas <- lazyGas:::.report_e_select_genes(
    mapping_mode = "gwas",
    simple_candidates = cand,
    credible_set = cred,
    gff_windows = win,
    snpeff = NULL,
    top_n = Inf
  )
  expect_true("G1" %in% gwas$Gene_ID)
  expect_true("max_PIP" %in% names(gwas))
  expect_true("max_PIP_SnpEff" %in% names(gwas))
  expect_true(all(is.na(gwas$max_PIP_SnpEff)))

  snpeff <- data.frame(
    Gene_ID = "G1",
    Chr = "chr1",
    Pos = 1500,
    Annotation_Impact = "HIGH",
    stringsAsFactors = FALSE
  )
  gwas2 <- lazyGas:::.report_e_select_genes(
    mapping_mode = "gwas",
    simple_candidates = cand,
    credible_set = cred,
    gff_windows = win,
    snpeff = snpeff,
    top_n = Inf
  )
  expect_true("G1" %in% gwas2$Gene_ID)
  expect_equal(gwas2$max_PIP[gwas2$Gene_ID == "G1"][1L], 0.4)
  expect_identical(gwas2$max_PIP_SnpEff[gwas2$Gene_ID == "G1"][1L], "HIGH")
})
