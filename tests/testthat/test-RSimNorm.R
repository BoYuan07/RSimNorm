context("Testing RSimNorm")
library(RSimNorm)
library(testthat)

data(Karen)

test_that("'RSimNorm' function provides normalized count",{
  set.seed(1)
  out = RSimNorm(Karen)
  P = out$P
  I0 = out$I0
  pi = out$pi0
  expect_equal(ncol(P),ncol(microbiome::abundances(Karen)))
})
