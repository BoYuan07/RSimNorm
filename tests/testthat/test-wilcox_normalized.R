context("Testing wilcoxon.normalized function")
library(RSimNorm)
library(testthat)

data(Karen)


test_that("'wilcox.normalized' function provides expected results", {
  set.seed(1)
  meta = sample_data(Karen)
  meta$pseudo = (runif(nrow(meta))>0.5)
  count = abundances(Karen)
  out = wilcox.normalized(count_table = count, meta = meta, main_var = 'pseudo')
  expect_equal(sum(out$is.DA==T), 0)
})
