.rotationMatrix2D <- function(theta) {
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
}

test_that("ProcrustesAlign collapses rotated/scaled/translated copies", {
  set.seed(1)
  base <- matrix(rnorm(10 * 2), 10, 2)
  configs <- list(
    a = base,
    b = (base %*% .rotationMatrix2D(0.7)) * 2.5 + 3,
    c = (base %*% .rotationMatrix2D(-1.2)) * 0.4 - 1
  )
  aligned <- ProcrustesAlign(configs)
  expect_equal(aligned$a, aligned$b, tolerance = 1e-6)
  expect_equal(aligned$a, aligned$c, tolerance = 1e-6)
  expect_named(aligned, names(configs))
})

test_that("RftraAlign recovers an exact similarity transform", {
  set.seed(2)
  base <- matrix(rnorm(10 * 2), 10, 2)
  configs <- list(
    a = base,
    b = sweep((base %*% .rotationMatrix2D(0.4)) * 1.8, 2, c(2, -1), "+")
  )
  aligned <- RftraAlign(configs)
  expect_equal(aligned$a, aligned$b, tolerance = 1e-5)
})

test_that("RftraAlign down-weights a single displaced landmark", {
  set.seed(3)
  base <- matrix(rnorm(12 * 2), 12, 2)
  displaced <- base
  displaced[1, ] <- displaced[1, ] + 10 # one grossly displaced landmark
  configs <- list(a = base, b = displaced)
  aligned <- RftraAlign(configs)
  undisplaced <- setdiff(seq_len(12), 1)
  # the undisplaced landmarks should superimpose far more closely than the
  # displaced one carries the change (the "Pinocchio" landmark)
  errUndisplaced <- max(abs(aligned$a[undisplaced, ] - aligned$b[undisplaced, ]))
  errDisplaced <- max(abs(aligned$a[1, ] - aligned$b[1, ]))
  expect_lt(errUndisplaced, errDisplaced)
})

test_that("ReadShapeConfigs parses the bundled sample shape file", {
  f <- system.file("extdata", "sample_shape_t00_d2_s0.4.tsv", package = "PhyloPCA")
  skip_if(!nzchar(f))
  configs <- ReadShapeConfigs(f, dim = 2)
  expect_true(length(configs) > 0)
  expect_equal(ncol(configs[[1]]), 2)
})
