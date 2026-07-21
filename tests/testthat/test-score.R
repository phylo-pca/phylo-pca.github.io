test_that("ScoreToTruth returns 0 CID/RF for an identical tree", {
  tr <- ape::rtree(10)
  scores <- ScoreToTruth(tr, tr, metrics = c("CID", "RF"))
  expect_equal(unname(scores["CID"]), 0, tolerance = 1e-8)
  expect_equal(unname(scores["RF"]), 0)
})

test_that("ScoreToTruth returns NA SPR for a non-binary tree, not an error", {
  tr <- ape::rtree(8)
  poly <- ape::di2multi(tr, tol = Inf) # collapse to a star (non-binary)
  scores <- ScoreToTruth(poly, tr, metrics = c("CID", "SPR"))
  expect_true(is.na(scores["SPR"]))
  expect_false(is.na(scores["CID"]))
})

test_that("ExactlyRecovered agrees with a CID of exactly 0", {
  tr <- ape::rtree(12)
  expect_true(ExactlyRecovered(tr, tr))
  other <- ape::rtree(12)
  # Vanishingly unlikely to coincide for 12 random tips
  expect_false(ExactlyRecovered(tr, other))
})

test_that("OptimalityRecovery compares within tolerance", {
  expect_true(OptimalityRecovery(10, 10))
  expect_true(OptimalityRecovery(10.00005, 10, tolerance = 1e-4))
  expect_false(OptimalityRecovery(10.1, 10, tolerance = 1e-4))
})

test_that("MongleStateProportions tabulates the bundled Mongle matrix", {
  f <- system.file("extdata", "mongle_2023.nex", package = "PhyloPCA")
  skip_if(!nzchar(f))
  props <- MongleStateProportions(f)
  expect_equal(sum(props), 107)
  expect_true(all(as.integer(names(props)) >= 2))
})

test_that("SimulateContinuous returns a correctly shaped, named matrix", {
  tr <- ape::rtree(9)
  set.seed(42)
  M <- SimulateContinuous(tr, nChar = 15)
  expect_equal(dim(M), c(9, 15))
  expect_equal(rownames(M), tr$tip.label)
})

test_that("SimulateDiscreteMS draws state sizes from the given proportions and stays variable", {
  tr <- ape::rtree(10)
  props <- c(`2` = 34, `3` = 44, `4` = 20, `5` = 8, `6` = 1)
  set.seed(1)
  X <- SimulateDiscreteMS(tr, nChar = 20, stateProportions = props)
  expect_equal(dim(X), c(10, 20))
  expect_equal(rownames(X), tr$tip.label)
  variable <- apply(X, 2, function(col) length(unique(col)) > 1)
  expect_true(all(variable))
})
