context("Testing corr.normalized function")
library(RSimNorm)
library(testthat)

data(Karen)


test_that("'corr.normalized' function provides expected results", {
  set.seed(1)
  meta = sample_data(Karen)
  meta$pseudo = runif(nrow(meta))
  count = abundances(Karen)
  out = corr.normalized(count_table = count, meta = meta, main_var = 'pseudo')
  expect_equal(sum(out$is.DA==T), 0)
})
