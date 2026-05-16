test_that("buildLazyGas documents lazygas_store parameter", {
  expect_true("lazygas_store" %in% names(formals(lazyGas::buildLazyGas)))
})

test_that("importLazyGasResults is exported", {
  expect_true("importLazyGasResults" %in% getNamespaceExports("lazyGas"))
})

test_that("LazyGas has store slot", {
  expect_true("store" %in% methods::slotNames("LazyGas"))
})
