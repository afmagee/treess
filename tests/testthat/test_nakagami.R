test_that("Test Nakagami distribution functions.", {
  nshape <- 0.3123
  nspread <- 1.656
  
  d_nak_jacob <- function(x, shape, spread) {
    scale <- nakagami2gamma(shape, spread)["scale"]
    dgamma(x*x, shape = shape, scale = scale) * 2 * x
  }
  
  testthat::expect_equal(
    d_nak_jacob(1:10, nshape, nspread),
    dnakagami(1:10, nshape, nspread)
  )
  
  f <- function(x) dnakagami(x, nshape, nspread)
  testthat::expect_equal(
    integrate(f, 0, Inf)$value,
    1.0
  )
  
  q <- 0.1231
  testthat::expect_equal(
    pnakagami(q, nshape, nspread),
    integrate(f, 0, q)$value
  )
  
  testthat::expect_equal(
    qnakagami(pnakagami(q, nshape, nspread), nshape, nspread),
    q
  )
})