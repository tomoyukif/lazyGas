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
  p <- lazyGas:::.variant_viewer_ggplot(plot_df = plot_df, cds_df = cds_df)
  expect_true(inherits(p, "ggplot"))
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
