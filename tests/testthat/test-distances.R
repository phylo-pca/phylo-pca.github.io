test_that("HammingDistance matches manual p-distance", {
  X <- rbind(a = c("0", "1", "1"), b = c("0", "0", "1"), c = c("1", "1", "0"))
  D <- HammingDistance(X)
  expect_equal(D["a", "b"], 1 / 3)
  expect_equal(D["a", "c"], 2 / 3)
  expect_equal(diag(D), c(a = 0, b = 0, c = 0))
})

test_that("CorrectedDistanceMS reduces to the binary JC correction at k=2", {
  set.seed(1)
  X <- matrix(sample(c("0", "1"), 21 * 40, replace = TRUE), nrow = 21,
              dimnames = list(paste0("t", 1:21), NULL))
  D <- CorrectedDistanceMS(X)
  p <- mean(X[1, ] != X[2, ])
  expect_equal(D[1, 2], -0.5 * log(1 - 2 * min(p, 0.5 - 1e-6)), tolerance = 1e-8)
})

test_that("SquaredEuclideanDistance is the elementwise square of dist()", {
  M <- matrix(rnorm(30), 10, 3)
  expect_equal(as.numeric(SquaredEuclideanDistance(M)), as.numeric(stats::dist(M))^2)
})

test_that("Nexus writers round-trip through ape", {
  M <- matrix(rnorm(21 * 5), 21, 5, dimnames = list(paste0("t", 1:21), NULL))
  f <- tempfile(fileext = ".nex")
  WriteContinuousNexus(M, f)
  lines <- readLines(f)
  expect_true(any(grepl("Continuous", lines)))
  expect_true(any(grepl("ntax=21 nchar=5", lines)))

  X <- matrix(sample(0:2, 21 * 8, replace = TRUE), 21, 8,
              dimnames = list(paste0("t", 1:21), NULL))
  g <- tempfile(fileext = ".nex")
  WriteDiscreteNexus(X, g)
  d <- ape::read.nexus.data(g)
  expect_equal(length(d), 21)
  expect_equal(length(d[[1]]), 8)
  expect_true(all(nchar(d[[1]]) == 1)) # one symbol per cell
})
