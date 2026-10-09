test_that("readAnnColumnMap loads demo map", {
  path <- system.file("extdata", "demo_ann_column_map.csv", package = "lazyGas")
  skip_if_not(nzchar(path) && file.exists(path))
  map <- readAnnColumnMap(path)
  expect_true(all(c("column", "class", "use_for_ranking") %in% names(map)))
  expect_true("Description" %in% map$column)
  expect_error(readAnnColumnMap(tempfile()), "not found")
})

test_that("readExpressionSampleMap loads demo map", {
  path <- system.file("extdata", "demo_expression_sample_map.csv", package = "lazyGas")
  skip_if_not(nzchar(path) && file.exists(path))
  map <- readExpressionSampleMap(path)
  expect_true("sample" %in% names(map))
  expect_true("tissue" %in% names(map))
  expect_true("fruit_mature_1" %in% map$sample)
})

test_that("annotation channel requires ann_column_map", {
  cand <- data.frame(
    Gene_ID = "g1",
    Description = "fruit weight",
    stringsAsFactors = FALSE
  )
  q <- phenotypeQuery("fruit weight", tissues = "fruit")
  expect_error(
    collectGeneEvidence(
      gene_ids = "g1",
      query = q,
      candidate = cand,
      sources = "annotation",
      use_cache = FALSE
    ),
    "ann_column_map"
  )
})

test_that("expression channel requires matrix and sample map", {
  q <- phenotypeQuery("fruit weight", tissues = "fruit")
  expect_error(
    collectGeneEvidence(
      gene_ids = "g1",
      query = q,
      candidate = data.frame(Gene_ID = "g1", stringsAsFactors = FALSE),
      sources = "expression",
      use_cache = FALSE
    ),
    "expression_matrix|expression_sample_map"
  )
})
