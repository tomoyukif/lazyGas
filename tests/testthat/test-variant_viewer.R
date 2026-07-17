make_mock_gff <- function(gene_id = "g1") {
  GenomicRanges::GRanges(
    seqnames = "1",
    ranges = IRanges::IRanges(start = 1000, end = 5000),
    strand = "+",
    type = "gene",
    ID = gene_id,
    gene_id = gene_id
  )
}

make_mock_snpeff_df <- function(gene_id = "g1") {
  data.frame(
    Allele = "A/G",
    Annotation = "missense_variant",
    Annotation_Impact = "MODERATE",
    Gene_Name = "GENE1",
    Gene_ID = gene_id,
    Feature_Type = "transcript",
    Feature_ID = "tx1",
    Transcript_BioType = "mRNA",
    Rank = "1/1",
    HGVS.c = "c.1A>G",
    HGVS.p = "p.M1V",
    `Pos.in.tx` = 1L,
    `Pos.in.CDS` = 1L,
    `Pos.in.AA` = 1L,
    Distance = 0L,
    INFO = "",
    Chr = "1",
    Pos = 1500L,
    negLog10P = 3,
    stringsAsFactors = FALSE
  )
}

test_that("variant viewer builds ggplot from mock data", {
  skip_if_not_installed("GenomicRanges")
  gff <- make_mock_gff("g1")
  cds <- GenomicRanges::GRanges(
    seqnames = "1",
    ranges = IRanges::IRanges(start = 1200, end = 1800),
    type = "CDS",
    Parent = "tx1",
    gene_id = "g1"
  )
  gff <- c(gff, cds)
  GenomeInfoDb::seqlevels(gff) <- "1"

  ann_df <- make_mock_snpeff_df("g1")
  scan_df <- data.frame(
    variant_ID = 1L,
    Chr = "1",
    Pos = 1500L,
    P.model = 0.001,
    negLog10P = 3,
    stringsAsFactors = FALSE
  )
  plot_df <- lazyGas:::.variant_viewer_build_plot_df(
    ann_df = ann_df,
    scan_df = scan_df,
    gene_id = "g1"
  )
  cds_df <- lazyGas:::.variant_viewer_cds_df(gff = gff, gene_id = "g1", plot_df = plot_df)
  p <- lazyGas:::.variant_viewer_ggplot(
    plot_df = plot_df,
    cds_df = cds_df,
    gene_id = "g1",
    pheno = "fruit_weight"
  )
  expect_true(inherits(p, "ggplot"))
  expect_match(p$labels$title, "fruit_weight")
  expect_match(p$labels$title, "g1")
})

test_that("variant_viewer_plot_title handles missing fields", {
  expect_equal(
    lazyGas:::.variant_viewer_plot_title(gene_id = "g1", pheno = "trait"),
    "g1 — trait"
  )
  expect_equal(lazyGas:::.variant_viewer_plot_title(pheno = "trait"), "trait")
  expect_equal(lazyGas:::.variant_viewer_plot_title(gene_id = "g1"), "g1")
  expect_null(lazyGas:::.variant_viewer_plot_title())
})

test_that("variant viewer tooltips distinguish coding and non-coding variants", {
  coding <- data.frame(
    Allele = "A/G",
    Position_in_CDS = 12L,
    Position_in_AA = 4L,
    Change_in_AA = "p.M1V",
    Start = 1500L,
    negLog10p = 3.1,
    P_value = 0.001,
    stringsAsFactors = FALSE
  )
  noncoding <- data.frame(
    Allele = "C/T",
    Position_in_CDS = NA_integer_,
    Position_in_AA = NA_integer_,
    Change_in_AA = NA_character_,
    Start = 1600L,
    negLog10p = 2.5,
    P_value = 0.01,
    stringsAsFactors = FALSE
  )

  coding_tip <- lazyGas:::.variant_viewer_variant_tooltip(coding)
  noncoding_tip <- lazyGas:::.variant_viewer_variant_tooltip(noncoding)
  gwas_tip <- lazyGas:::.variant_viewer_gwas_tooltip(coding)

  expect_match(coding_tip, "Allele: A/G")
  expect_match(coding_tip, "CDS pos: 12")
  expect_match(coding_tip, "AA pos: 4")
  expect_match(coding_tip, "AA change: p.M1V")
  expect_equal(noncoding_tip, "Allele: C/T")
  expect_match(gwas_tip, "Pos: 1500")
  expect_match(gwas_tip, "negLog10P:")
})

test_that("variant viewer ggplotly uses custom hover text", {
  skip_if_not_installed("plotly")
  df <- data.frame(
    Allele = "A/G",
    Position_in_CDS = 1L,
    Position_in_AA = 1L,
    Change_in_AA = "p.M1V",
    Change_in_DNA = "c.1A>G",
    Start = 1500L,
    End = 1500L,
    ID = 1L,
    P_value = 0.001,
    negLog10p = 3,
    Chr = "1",
    Gene_ID = "g1",
    Transcript_ID = "tx1",
    Effect = factor("MODERATE", levels = c("HIGH", "MODERATE", "LOW", "MODIFIER")),
    stringsAsFactors = FALSE
  )
  df$y_pos <- 1
  df$x_pos <- 1500
  df$ymin <- 1
  df$ymax <- 1.8
  cds_df <- data.frame(
    Transcript_ID = factor("tx1", levels = "tx1"),
    cds_start = 1200L,
    cds_end = 1800L,
    ymin = 1,
    ymax = 1.8,
    stringsAsFactors = FALSE
  )
  p <- lazyGas:::.variant_viewer_ggplot(
    plot_df = df,
    cds_df = cds_df,
    gene_id = "g1",
    pheno = "fruit_weight"
  )
  pg <- plotly::ggplotly(p, tooltip = "text")
  texts <- unlist(lapply(pg$x$data, `[[`, "text"))
  texts <- texts[nzchar(texts)]
  expect_true(any(grepl("Allele:", texts)))
  expect_true(any(grepl("negLog10P:", texts)))
})

test_that("variant_viewer_gene_info handles NA gene_id metadata", {
  skip_if_not_installed("GenomicRanges")
  gff <- GenomicRanges::GRanges(
    seqnames = c("2", "2"),
    ranges = IRanges::IRanges(start = c(1000L, 1000L), end = c(5000L, 5000L)),
    strand = c("+", "+"),
    type = c("gene", "mRNA"),
    ID = c("g1", "g1.t1"),
    gene_id = c(NA_character_, NA_character_)
  )
  info <- lazyGas:::.variant_viewer_gene_info(gff = gff, gene_id = "g1")
  expect_equal(info$Gene_ID, "g1")
  expect_equal(info$gene_chr, "2")
})

test_that("variant_viewer_safe_id sanitizes gene ids", {
  expect_equal(lazyGas:::.variant_viewer_safe_id("gene/A"), "gene_A")
})

test_that("getVariantViewerData and plotVariantViewer are exported", {
  skip_if_not_installed("lazyGas")
  expect_true("getVariantViewerData" %in% getNamespaceExports("lazyGas"))
  expect_true("plotVariantViewer" %in% getNamespaceExports("lazyGas"))
  expect_true("makeInteractiveDashboard" %in% getNamespaceExports("lazyGas"))
})
