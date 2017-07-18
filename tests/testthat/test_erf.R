context("Error functions")

test_that("classes are correct",{
  x <- 1:10
  xc <- x + 1i
  expect_is(erf_re(x), 'numeric')
  expect_is(erfc_re(x), 'numeric')
  expect_error(erf_re(xc))
  expect_error(erfc_re(xc))
  expect_is(RcppFaddeeva::erf(x), 'complex')
  expect_is(RcppFaddeeva::erf(xc), 'complex')
  expect_is(RcppFaddeeva::erfc(x), 'complex')
  expect_is(RcppFaddeeva::erfc(xc), 'complex')
})

test_that("erf agrees",{
  x <- seq(0,1,by=0.1)
  expect_equal(Re(RcppFaddeeva::erf(x)), erf_re(x))
})

test_that("erfc agrees",{
  x <- seq(0,1,by=0.1)
  expect_equal(Re(RcppFaddeeva::erfc(x)), erfc_re(x))
})

test_that("erfc is consistent",{
  x <- seq(0,1,by=0.1)
  expect_equal(RcppFaddeeva::erfc(x), 1 - RcppFaddeeva::erf(x))
  expect_equal(erfc_re(x), 1 - erf_re(x))
  
  expect_equal(RcppFaddeeva::erfc(-Inf), 2+0i)
  expect_equal(RcppFaddeeva::erfc(0), 1+0i)
  expect_equal(RcppFaddeeva::erfc(Inf), 0+0i)
  
  # special cases
  expect_equal(erfc_re(-Inf), 2)
  expect_equal(erfc_re(0), 1)
  expect_equal(erfc_re(Inf), 0)
  
  # identity
  expect_equal(erfc_re(-x), 2 - erfc_re(x))
  expect_equal(RcppFaddeeva::erfc(-x), 2 - RcppFaddeeva::erfc(x))
})

test_that("ierfc is consistent",{
  x <- seq(0,1,by=0.1)

  # special cases
  expect_equal(ierfc_re(0), Inf)
  expect_equal(ierfc_re(1), 0)
  expect_equal(ierfc_re(2), -Inf)
  
  # identity
  expect_equal(ierfc_re(1-x), ierf_re(x))
})