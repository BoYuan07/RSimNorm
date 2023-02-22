context("Testing t.normalized function")
library(RSimNorm)
library(testthat)

data(Karen)


test_that("'t.normalized' function provides expected results", {
  set.seed(1)
  meta = sample_data(Karen)
  meta$pseudo = (runif(nrow(meta))>0.5)
  count = abundances(Karen)
  out = t.normalized(count_table = count, meta = meta, main_var = 'pseudo')
  expect_equal(sum(out$is.DA==T), 0)
})
