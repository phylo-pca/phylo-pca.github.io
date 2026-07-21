test_that(".TntToNewick converts TNT parenthetical notation to parseable Newick", {
  tnt <- "*(0 (1 2) (3 (4 5)))"
  newick <- PhyloPCA:::.TntToNewick(tnt)
  tr <- ape::read.tree(text = newick)
  expect_equal(ape::Ntip(tr), 6)
  expect_true(ape::is.binary(ape::unroot(tr)) || !ape::is.rooted(tr))
})

test_that("ScoreMPTsByMaxEntropy uses every tree when n <= k", {
  skip_if_not_installed("MaxMin")
  set.seed(1)
  truth <- ape::rtree(8)
  trees <- lapply(1:3, function(i) ape::rtree(8, tip.label = truth$tip.label))
  class(trees) <- "multiPhylo"
  result <- ScoreMPTsByMaxEntropy(trees, truth, k = 10)
  expect_equal(result$nTree, 3)
  expect_equal(result$nSel, 3)
  expect_true(is.finite(result$CID))
})

test_that("FindTNT errors informatively when no binary is configured", {
  oldOpt <- getOption("PhyloPCA.tnt")
  oldEnv <- Sys.getenv("TNT_PATH", NA)
  options(PhyloPCA.tnt = "/nonexistent/tnt")
  Sys.setenv(TNT_PATH = "/nonexistent/tnt")
  on.exit({
    options(PhyloPCA.tnt = oldOpt)
    if (is.na(oldEnv)) Sys.unsetenv("TNT_PATH") else Sys.setenv(TNT_PATH = oldEnv)
  })
  # Only meaningful when none of the built-in fallback paths exist either.
  found <- tryCatch(FindTNT(), error = function(e) NULL)
  if (is.null(found)) {
    expect_error(FindTNT(), "Could not locate TNT")
  } else {
    succeed("TNT is installed on this machine; nothing to assert")
  }
})
