test_that("recalc group retain keeps member-block markers under the representative", {
  # Two correlated peaks: rep lead A (taller) absorbs member lead B.
  # Member block carries marker B_cds that must survive recalc tables.
  # 14.2 / 12.4 < 2 → fold gate passes.
  peakcall <- data.frame(
    peak_ID = c(1L, 1L, 1L, 5L, 5L, 5L),
    peak_variant_ID = c("A", "A", "A", "B", "B", "B"),
    variant_ID = c("A", "A_flank", "shared", "B", "B_cds", "shared"),
    dist2peak = c(0, 100, 50, 0, 10, 40),
    LD2peak = c(1, 0.9, 0.8, 1, 0.95, 0.7),
    chr_start = 1L,
    chr_end = 200L,
    FDR = c(1e-10, 1e-8, 1e-7, 1e-9, 1e-8, 1e-7),
    negLog10P = c(14.2, 10, 9, 12.4, 12.0, 8),
    stringsAsFactors = FALSE
  )
  new_peaks <- list(
    newpeaks = "A",
    member = list(c("A", "B"))
  )

  retained <- lazyGas:::.recalc_group_retained_tables(
    peakcall = peakcall,
    new_peaks = new_peaks,
    block_stat_cols = c("FDR", "negLog10P")
  )

  expect_equal(nrow(retained$peaks), 1L)
  expect_equal(as.character(retained$peaks$peak_variant_ID), "A")
  expect_equal(retained$peaks$peak_ID, 1L)

  vids <- as.character(retained$blocks$variant_ID)
  expect_true(all(c("A", "B", "B_cds", "A_flank", "shared") %in% vids))
  expect_true(all(retained$blocks$peak_ID == 1L))
  expect_true(all(as.character(retained$blocks$peak_variant_ID) == "A"))
  # Dedup shared marker: prefer row from representative block
  shared <- retained$blocks[retained$blocks$variant_ID == "shared", , drop = FALSE]
  expect_equal(nrow(shared), 1L)
  expect_equal(shared$LD2peak, 0.8)
})

test_that("recalc group retain drops same-chr member when fold >= 2", {
  peakcall <- data.frame(
    peak_ID = c(1L, 1L, 5L, 5L),
    peak_variant_ID = c("A", "A", "B", "B"),
    variant_ID = c("A", "A_flank", "B", "B_cds"),
    Chr = c("chr04", "chr04", "chr04", "chr04"),
    Pos = c(200L, 100L, 300L, 250L),
    dist2peak = 0,
    LD2peak = 1,
    chr_start = 0L,
    chr_end = 0L,
    FDR = 1e-8,
    negLog10P = c(20, 15, 5, 4),
    stringsAsFactors = FALSE
  )
  new_peaks <- list(newpeaks = "A", member = list(c("A", "B")))
  retained <- lazyGas:::.recalc_group_retained_tables(
    peakcall = peakcall,
    new_peaks = new_peaks,
    block_stat_cols = c("FDR", "negLog10P"),
    group_retain_fold = 2
  )
  expect_equal(as.character(retained$peaks$peak_variant_ID), "A")
  vids <- as.character(retained$blocks$variant_ID)
  expect_true(all(c("A", "A_flank") %in% vids))
  expect_false(any(vids %in% c("B", "B_cds")))
})

test_that("recalc group retain merges same-chr and keeps cross-chr peaks independent", {
  peakcall <- data.frame(
    peak_ID = c(1L, 1L, 5L, 5L, 9L, 9L),
    peak_variant_ID = c("A", "A", "B", "B", "C", "C"),
    variant_ID = c("A", "A_flank", "B", "B_cds", "C", "C_flank"),
    Chr = c("chr04", "chr04", "chr04", "chr04", "chr01", "chr01"),
    Pos = c(200L, 100L, 300L, 250L, 500L, 400L),
    dist2peak = 0,
    LD2peak = 1,
    chr_start = 0L,
    chr_end = 0L,
    FDR = 1e-8,
    negLog10P = 10,
    stringsAsFactors = FALSE
  )
  new_peaks <- list(newpeaks = "A", member = list(c("A", "B", "C")))
  retained <- lazyGas:::.recalc_group_retained_tables(
    peakcall = peakcall,
    new_peaks = new_peaks,
    block_stat_cols = c("FDR", "negLog10P")
  )
  expect_true(all(c(1L, 9L) %in% retained$peaks$peak_ID))
  expect_false(5L %in% retained$peaks$peak_ID) # B merged into A
  a_blocks <- retained$blocks[retained$blocks$peak_ID == 1L, , drop = FALSE]
  c_blocks <- retained$blocks[retained$blocks$peak_ID == 9L, , drop = FALSE]
  expect_true(all(c("A", "A_flank", "B", "B_cds") %in% a_blocks$variant_ID))
  expect_false(any(a_blocks$variant_ID %in% c("C", "C_flank")))
  expect_true(all(c("C", "C_flank") %in% c_blocks$variant_ID))
  expect_true(all(as.character(c_blocks$peak_ID) == "9"))
})

test_that("recalc group retain drops cross-chr independent peak when fold >= 2", {
  peakcall <- data.frame(
    peak_ID = c(1L, 1L, 9L, 9L),
    peak_variant_ID = c("A", "A", "C", "C"),
    variant_ID = c("A", "A_flank", "C", "C_flank"),
    Chr = c("chr04", "chr04", "chr01", "chr01"),
    Pos = c(200L, 100L, 500L, 400L),
    dist2peak = 0,
    LD2peak = 1,
    chr_start = 0L,
    chr_end = 0L,
    FDR = 1e-8,
    negLog10P = c(20, 15, 5, 4),
    stringsAsFactors = FALSE
  )
  new_peaks <- list(newpeaks = "A", member = list(c("A", "C")))
  retained <- lazyGas:::.recalc_group_retained_tables(
    peakcall = peakcall,
    new_peaks = new_peaks,
    block_stat_cols = c("FDR", "negLog10P"),
    group_retain_fold = 2
  )
  expect_equal(retained$peaks$peak_ID, 1L)
  expect_false(9L %in% retained$peaks$peak_ID)
  expect_false(any(retained$blocks$variant_ID %in% c("C", "C_flank")))
})

test_that("recalc groups table maps absorbed leads via grouped_with", {
  peakcall <- data.frame(
    peak_ID = c(1L, 1L, 5L, 5L),
    peak_variant_ID = c("A", "A", "B", "B"),
    variant_ID = c("A", "A_flank", "B", "B_cds"),
    FDR = c(1e-10, 1e-8, 1e-9, 1e-8),
    negLog10P = c(14, 10, 12, 11),
    stringsAsFactors = FALSE
  )
  new_peaks <- list(newpeaks = "A", member = list(c("A", "B")))
  peak_grp <- list(
    pvalue = data.frame(
      variant_ID = c("B"),
      p_values = 0.89,
      stringsAsFactors = FALSE
    )
  )

  groups <- lazyGas:::.recalc_groups_table(
    peakcall = peakcall,
    new_peaks = new_peaks,
    peak_grp = peak_grp,
    from_candidates = c("FDR", "P.model")
  )

  expect_equal(sort(as.character(groups$variant_ID)), c("A", "B"))
  expect_equal(groups$grouped_with[groups$variant_ID == "B"], 1L)
  expect_equal(groups$full_vs_reduce_pvalue[groups$variant_ID == "B"], 0.89)
})

test_that("recalc union member blocks merges into refined peak", {
  refined <- data.frame(
    peak_ID = 1L,
    peak_variant_ID = "A_new",
    source_peak_variant_ID = "A",
    variant_ID = c("A_new", "A_flank"),
    dist2peak = c(0, 20),
    LD2peak = c(1, 0.9),
    chr_start = 1L,
    chr_end = 100L,
    P.model = c(1e-14, 1e-8),
    negLog10P = c(14, 8),
    stringsAsFactors = FALSE
  )
  peakcall <- data.frame(
    peak_ID = c(1L, 1L, 5L, 5L),
    peak_variant_ID = c("A", "A", "B", "B"),
    variant_ID = c("A", "A_flank", "B", "B_cds"),
    dist2peak = c(0, 100, 0, 10),
    LD2peak = c(1, 0.9, 1, 0.95),
    chr_start = 1L,
    chr_end = 200L,
    P.model = c(1e-12, 1e-8, 1e-11, 1e-10),
    negLog10P = c(12, 8, 11, 10),
    stringsAsFactors = FALSE
  )
  new_peaks <- list(newpeaks = "A", member = list(c("A", "B")))

  merged <- lazyGas:::.recalc_union_member_blocks(
    blocks = refined,
    peakcall = peakcall,
    new_peaks = new_peaks,
    block_stat_cols = c("P.model", "negLog10P")
  )

  expect_true(all(c("A_new", "A_flank", "B", "B_cds") %in% merged$variant_ID))
  expect_true(all(merged$peak_ID == 1L))
  expect_equal(nrow(merged), length(unique(merged$variant_ID)))
})

test_that("recalc union member blocks drops member when fold >= 2", {
  refined <- data.frame(
    peak_ID = 1L,
    peak_variant_ID = "A_new",
    source_peak_variant_ID = "A",
    variant_ID = c("A_new", "A_flank"),
    dist2peak = c(0, 20),
    LD2peak = c(1, 0.9),
    chr_start = 1L,
    chr_end = 100L,
    P.model = c(1e-20, 1e-8),
    negLog10P = c(20, 8),
    stringsAsFactors = FALSE
  )
  peakcall <- data.frame(
    peak_ID = c(1L, 1L, 5L, 5L),
    peak_variant_ID = c("A", "A", "B", "B"),
    variant_ID = c("A", "A_flank", "B", "B_cds"),
    dist2peak = c(0, 100, 0, 10),
    LD2peak = c(1, 0.9, 1, 0.95),
    chr_start = 1L,
    chr_end = 200L,
    P.model = c(1e-20, 1e-8, 1e-5, 1e-4),
    negLog10P = c(20, 8, 5, 4),
    stringsAsFactors = FALSE
  )
  new_peaks <- list(newpeaks = "A", member = list(c("A", "B")))
  merged <- lazyGas:::.recalc_union_member_blocks(
    blocks = refined,
    peakcall = peakcall,
    new_peaks = new_peaks,
    block_stat_cols = c("P.model", "negLog10P"),
    group_retain_fold = 2
  )
  expect_true(all(c("A_new", "A_flank") %in% merged$variant_ID))
  expect_false(any(merged$variant_ID %in% c("B", "B_cds")))
})
