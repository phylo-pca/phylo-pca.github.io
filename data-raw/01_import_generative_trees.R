#!/usr/bin/env Rscript
# The 250 generative ("true") trees used throughout Smith (2026), sampled
# by SampleStableTrees() from the RevBayes re-inference of the Mongle
# (2023) matrix (see RunRevBayesReinference(), inst/rev/reinfer_variable.Rev).
# Bundling the sampled trees directly spares most users the multi-hour
# re-inference step; see vignette("01-generative-trees") to reproduce them
# from scratch.
suppressMessages(library(ape))
SCALE <- "C:/Users/pjjg18/GitHub/hominin-pca-scale/out/scale/trees"
CACHE <- file.path("data-raw", "cache", "trees")

if (dir.exists(SCALE)) {
  dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
  invisible(file.copy(list.files(SCALE, full.names = TRUE), CACHE, overwrite = TRUE))
}
files <- sort(list.files(CACHE, pattern = "^tree_[0-9]+\\.nwk$", full.names = TRUE))
stopifnot(length(files) > 0)
generative_trees <- lapply(files, read.tree)
class(generative_trees) <- "multiPhylo"
names(generative_trees) <- sprintf("tree_%03d", seq_along(generative_trees))
save(generative_trees, file = file.path("data", "generative_trees.rda"), compress = "xz")
cat(sprintf("data/generative_trees.rda (%d trees)\n", length(generative_trees)))
