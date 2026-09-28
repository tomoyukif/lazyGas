test_that("GO OBO parse and resolve drops unresolved IDs", {
  obo <- tempfile(fileext = ".obo")
  writeLines(c(
    "[Term]",
    "id: GO:0003844",
    "name: 1,4-alpha-glucan branching enzyme activity",
    "",
    "[Term]",
    "id: GO:000001",
    "name: padded test term"
  ), obo)
  on.exit(unlink(obo), add = TRUE)
  cache <- tempfile("go_cache_")
  dir.create(cache)
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)

  map <- lazyGas:::.go_term_name_map(cache_dir = cache, obo_path = obo, download = FALSE)
  expect_true("GO:0003844" %in% names(map))
  expect_true("GO:0000001" %in% names(map)) # zero-padded

  terms <- lazyGas:::.e5_resolve_term_names(
    c("GO:0003844", "GO:9999999", "starch synthase"),
    cache_dir = cache,
    obo_path = obo,
    download = FALSE
  )
  expect_true("1,4-alpha-glucan branching enzyme activity" %in% terms)
  expect_true("starch synthase" %in% terms)
  expect_false(any(grepl("^GO:", terms)))
})

test_that("KEGG list parse and resolve", {
  cache <- tempfile("kegg_cache_")
  dir.create(cache)
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  text <- paste(
    "path:map00500\tStarch and sucrose metabolism",
    "K00700\tstarch synthase",
    "ec:2.4.1.21\tstarch synthase (EC 2.4.1.21)",
    sep = "\n"
  )
  mp <- lazyGas:::.kegg_list_name_map(
    "pathway", cache_dir = cache, text = text, download = FALSE
  )
  expect_equal(unname(mp["map00500"]), "Starch and sucrose metabolism")

  # seed kind maps via text for ko/enzyme too
  lazyGas:::.kegg_list_name_map("ko", cache_dir = cache, text = "K00700\tstarch synthase\n", download = FALSE)
  lazyGas:::.kegg_list_name_map("enzyme", cache_dir = cache, text = "2.4.1.21\tstarch synthase (EC 2.4.1.21)\n", download = FALSE)

  entries <- lazyGas:::.e6_resolve_entries(
    c("map00500", "K00700", "2.4.1.21", "osa:12345", "free name"),
    cache_dir = cache,
    download = FALSE
  )
  kinds <- vapply(entries, function(e) as.character(e$kind %||% "null"), character(1L))
  names_e <- vapply(entries, function(e) e$name, character(1L))
  expect_true("pathway" %in% kinds)
  expect_true("ko" %in% kinds)
  expect_true("ec" %in% kinds)
  expect_true("free name" %in% names_e)
  expect_false(any(grepl("^osa:", names_e)))
})

test_that("Pfam / InterPro list parse and resolve", {
  cache <- tempfile("dom_cache_")
  dir.create(cache)
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  pfam <- file.path(cache, "pfam.tsv")
  writeLines("PF00128\tGlyco_transf_8\tCL0000\tClan\tGlycosyl transferase family 8", pfam)
  ipr <- file.path(cache, "entry.list")
  writeLines(c(
    "ENTRY_AC\tENTRY_TYPE\tENTRY_NAME",
    "IPR000757\tDomain\tGlycosyl transferase, family 1"
  ), ipr)

  doms <- lazyGas:::.e7_resolve_domains(
    c("PF00128", "IPR000757", "unknownPF99999", "Named domain only"),
    cache_dir = cache,
    pfam_path = pfam,
    interpro_path = ipr,
    download = FALSE
  )
  sources <- vapply(doms, `[[`, character(1L), "source")
  names_d <- vapply(doms, `[[`, character(1L), "name")
  expect_true("pfam" %in% sources)
  expect_true("interpro" %in% sources)
  expect_true("Named domain only" %in% names_d)
  expect_false(any(grepl("^PF|^IPR", names_d)))
})

test_that("E2 domain AA overlap same protein ID only", {
  snp <- data.frame(
    Gene_ID = "G1",
    Annotation_Impact = c("MODERATE", "HIGH", "LOW"),
    Annotation = c("missense_variant", "stop_gained", "synonymous"),
    Feature_ID = c("TX1", "TX1", "TX2"),
    `Pos.in.AA` = c("150", "160", "150"),
    `HGVS.p` = c("p.A50V", "p.Q53*", "p.A50A"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  ip <- data.frame(
    protein_id = c("TX1", "TX2"),
    start = c(100, 100),
    stop = c(200, 200),
    description = c("DomA", "DomB"),
    stringsAsFactors = FALSE
  )
  pmap <- data.frame(
    Gene_ID = c("G1", "G1"),
    transcript_id = c("TX1", "TX2"),
    protein_id = c("TX1", "TX2"),
    stringsAsFactors = FALSE
  )
  ov <- lazyGas:::.e2_domain_overlaps(snp, ip, pmap, "G1")
  expect_equal(length(ov$domain_overlaps), 2L)
  # DomA (TX1) should see MODERATE and HIGH → HIGH
  d1 <- ov$domain_overlaps[[1L]]
  expect_identical(d1$worst_overlapping_impact, "HIGH")
  expect_true(length(ov$high_effects) >= 1L)
  expect_true(is.finite(ov$domain_id_match_rate))
  expect_equal(ov$domain_id_match_rate, 1)
})

test_that("E2 does not cross-isoform hit TX2 domain with TX1 variant", {
  snp <- data.frame(
    Gene_ID = "G1",
    Annotation_Impact = "MODERATE",
    Annotation = "missense_variant",
    Feature_ID = "TX1",
    `Pos.in.AA` = "150",
    `HGVS.p` = "p.A50V",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  # Only TX2 has a domain covering AA 150; TX1 variant must not hit it
  ip <- data.frame(
    protein_id = "TX2",
    start = 100,
    stop = 200,
    description = "DomB",
    stringsAsFactors = FALSE
  )
  pmap <- data.frame(
    Gene_ID = c("G1", "G1"),
    transcript_id = c("TX1", "TX2"),
    protein_id = c("TX1", "TX2"),
    stringsAsFactors = FALSE
  )
  # TX1 has aa_pos but no InterPro row → match rate 0 and warn; no false DomB hit
  expect_warning(
    ov <- lazyGas:::.e2_domain_overlaps(snp, ip, pmap, "G1"),
    "InterProScan ID match rate"
  )
  expect_equal(length(ov$domain_overlaps), 0L)
  expect_equal(ov$domain_id_match_rate, 0)
})

test_that("E2 warns when InterProScan ID match rate is below 80%", {
  snp <- data.frame(
    Gene_ID = "G1",
    Annotation_Impact = c("MODERATE", "MODERATE", "MODERATE", "MODERATE", "MODERATE"),
    Annotation = rep("missense_variant", 5L),
    Feature_ID = c("TX1", "TX_MISS", "TX_MISS", "TX_MISS", "TX_MISS"),
    `Pos.in.AA` = c("150", "10", "20", "30", "40"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  ip <- data.frame(
    protein_id = "TX1",
    start = 100,
    stop = 200,
    description = "DomA",
    stringsAsFactors = FALSE
  )
  pmap <- data.frame(
    Gene_ID = "G1",
    transcript_id = "TX1",
    protein_id = "TX1",
    stringsAsFactors = FALSE
  )
  expect_warning(
    ov <- lazyGas:::.e2_domain_overlaps(snp, ip, pmap, "G1"),
    "InterProScan ID match rate"
  )
  expect_true(is.finite(ov$domain_id_match_rate))
  expect_lt(ov$domain_id_match_rate, 0.8)
  expect_match(ov$warning, "80%", fixed = TRUE)
  # Only TX1 variant may overlap DomA
  expect_equal(length(ov$domain_overlaps), 1L)
  expect_identical(ov$domain_overlaps[[1L]]$worst_overlapping_impact, "MODERATE")
})

test_that("simple candidate helper columns", {
  cand <- data.frame(
    peak_ID = 1L,
    Gene_ID = "G1",
    Gene_chr = "chr1",
    Gene_start = 100L,
    dist2peak = 10,
    negLog10P = 5,
    Description = "extra",
    stringsAsFactors = FALSE
  )
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(100, 500),
    strand = "+",
    type = "gene",
    ID = "G1"
  )
  simple <- lazyGas:::.simple_candidate_from_wide(cand, gff = gr)
  expect_true(all(c("Gene_ID", "Gene_chr", "Gene_start", "Gene_end") %in% names(simple)))
  expect_equal(simple$Gene_end[1L], 500L)
  expect_false("Description" %in% names(simple))
})
