context("Error functions")

x <- 1:10
xc <- x + 1i
test_that("classes are correct",{
          expect_is(erf_re(x), 'numeric')
          expect_error(erf_re(xc))
          expect_is(RcppFaddeeva::erf(x), 'complex')
          expect_is(RcppFaddeeva::erf(xc), 'complex')
})
