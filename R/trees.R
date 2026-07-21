#' Re-infer a posterior of generative trees
#'
#' Runs the RevBayes re-inference used throughout Smith (2026) to obtain a
#' posterior of plausible generative trees for simulation: the Mk model
#' (Lewis 2001), matching Raskin et al. (2026)'s own model exactly except
#' for the ascertainment-bias correction (`coding = "variable"`, the
#' correct choice for a matrix with no invariant characters, such as
#' Mongle et al. (2023)'s; see [MongleStateProportions()]). Two independent
#' runs are written so that topological convergence can subsequently be
#' assessed with [TreeESS()]. This step is CPU-intensive (hours, even for
#' 21 taxa) - `generative_trees` provides the actual 250-tree sample used
#' in Smith (2026), sparing most users from re-running it.
#'
#' @param nexusFile Path to a discrete character Nexus file.
#' @param outDir Directory for RevBayes' logs and tree traces
#'   (`reinfer.log`, `reinfer_run_1.trees`, `reinfer_run_2.trees`).
#' @param rbPath Path to the RevBayes executable; see [FindRevBayes()].
#' @param nGen Number of MCMC generations per run.
#' @param printGen Sampling frequency, in generations.
#'
#' @return A list with elements `runFiles` (the two `.trees` files
#'   produced, suitable input to [SampleStableTrees()]) and `log`.
#'
#' @export
RunRevBayesReinference <- function(nexusFile, outDir, rbPath = FindRevBayes(),
                                    nGen = 250000, printGen = 50) {
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
  template <- system.file("rev", "reinfer_variable.Rev", package = "phyloPCA")
  header <- tempfile(fileext = ".Rev")
  writeLines(c(
    sprintf('data = "%s"', nexusFile),
    sprintf('outDir = "%s"', outDir),
    sprintf("nGen = %d", nGen),
    sprintf("printGen = %d", printGen),
    sprintf('source("%s")', template)
  ), header)
  log <- file.path(outDir, "reinfer.rb.log")
  status <- system2(rbPath, shQuote(header), stdout = log, stderr = log)
  if (status != 0) warning("RevBayes exited with status ", status, "; see ", log)
  # Return: the two run trace files (input to SampleStableTrees()) + log
  list(runFiles = file.path(outDir, c("reinfer_run_1.trees", "reinfer_run_2.trees")),
       log = log)
}

#' Read a RevBayes tree trace
#'
#' Reads a RevBayes `.trees` log (tab-delimited, one Newick string per
#' generation in the last column, `[&...]` clade annotations stripped) and
#' discards a proportional burn-in.
#'
#' @param file Path to a RevBayes `.trees` file.
#' @param burnin Proportion of samples to discard from the start of the trace.
#'
#' @return A `multiPhylo` of post-burnin trees, or `NULL` if the file
#'   has fewer than four sampled trees.
#'
#' @importFrom ape read.tree
#' @export
ReadRevBayesTrees <- function(file, burnin = 0.25) {
  x <- readLines(file)
  if (length(x) < 3) {
    # Return: too few lines to contain a header and any sampled trees
    return(NULL)
  }
  x <- x[-1] # drop header row
  newick <- vapply(strsplit(x, "\t"), function(p) p[length(p)], character(1))
  newick <- gsub("\\[&[^]]*\\]", "", newick)
  trees <- read.tree(text = newick)
  n <- length(trees)
  if (n < 4) {
    # Return: too few post-header samples to be a usable trace
    return(NULL)
  }
  keep <- seq(max(1, floor(burnin * n) + 1), n)
  # Return: post-burnin trees
  trees[keep]
}

#' Topological effective sample size across RevBayes runs
#'
#' Estimates the effective sample size of tree topology across one or more
#' independent RevBayes runs, using the lower of the Frechet correlation and
#' median pseudo-ESS statistics (the pair used throughout Smith 2026), as
#' recommended by Fabreti & Hohna (2022) and Magee et al. (2024).
#'
#' @param runFiles Character vector of paths to RevBayes `.trees` files, one
#'   per independent run.
#' @param burnin Proportion of samples discarded from each trace.
#' @param cap Maximum number of trees retained per run (traces longer than
#'   this are evenly thinned before the ESS calculation, which scales poorly
#'   with sample size).
#' @param threshold Minimum summed ESS (across runs) for convergence.
#'
#' @return A list with elements `ok` (logical, whether `threshold` is met),
#'   `perRun` (a data frame of per-run ESS estimates), `summed` (the two
#'   ESS statistics summed across runs) and `minEss` (the smaller of the two).
#'
#' @importFrom ape unroot
#' @export
TreeESS <- function(runFiles, burnin = 0.25, cap = 1000, threshold = 200) {
  if (!requireNamespace("treess", quietly = TRUE)) {
    stop("Package 'treess' is required: ",
         "remotes::install_github('ms609/treess')")
  }
  if (!requireNamespace("TreeDist", quietly = TRUE)) {
    stop("Package 'TreeDist' is required for RobinsonFoulds distances.")
  }
  runs <- lapply(runFiles, function(f) {
    tr <- ReadRevBayesTrees(f, burnin)
    if (is.null(tr)) {
      # Return: unusable run, dropped below
      return(NULL)
    }
    if (length(tr) > cap) {
      tr <- tr[unique(round(seq(1, length(tr), length.out = cap)))]
    }
    # Return: thinned post-burnin trees for this run
    tr
  })
  runs <- Filter(Negate(is.null), runs)
  if (!length(runs)) {
    # Return: no usable runs to assess
    return(list(ok = FALSE, reason = "no usable runs"))
  }
  ess <- do.call(rbind, treess::treess(
    runs, TreeDist::RobinsonFoulds,
    methods = c("frechetCorrelationESS", "medianPseudoESS")
  ))
  ess <- as.data.frame(ess)
  cols <- c("frechetCorrelationESS", "medianPseudoESS")
  summed <- colSums(ess[, cols, drop = FALSE])
  # Return: convergence verdict + supporting detail
  list(ok = min(summed) >= threshold, perRun = ess, summed = summed,
       minEss = min(summed), nRun = length(runs),
       nSamp = vapply(runs, length, integer(1)))
}

#' Sample a stable set of generative trees from a posterior
#'
#' Draws a fixed, append-only set of "truth" trees from two independent
#' post-burnin RevBayes tree traces, for use as the generative tree in
#' simulation. Runs are interleaved (so both chains contribute from the
#' first tree onward) and then strided, so that `tree_k` is a fixed function
#' of `k`: raising `nTrees` later only appends new trees and never changes
#' the identity of an existing one.
#'
#' @inheritParams TreeESS
#' @param nTrees Number of trees to sample.
#' @param stride Spacing (in interleaved-pool samples) between sampled trees;
#'   should be odd, and large enough that consecutive draws are
#'   effectively independent (checked via [TreeESS()] when `checkEss = TRUE`).
#' @param checkEss Whether to require the input traces to have converged
#'   (summed topological ESS >= `threshold`) before sampling.
#'
#' @return A `multiPhylo` of `nTrees` sampled trees.
#'
#' @export
SampleStableTrees <- function(runFiles, nTrees, stride = 11L, burnin = 0.25,
                               threshold = 200, checkEss = TRUE) {
  stopifnot(length(runFiles) >= 2)
  minEss <- Inf
  if (checkEss) {
    te <- TreeESS(runFiles, burnin = burnin, threshold = threshold)
    if (!isTRUE(te$ok)) {
      stop(sprintf("tree-ESS %.0f < %d - extend the run before sampling.",
                    te$minEss, threshold))
    }
    minEss <- te$minEss
  }

  pools <- lapply(runFiles, ReadRevBayesTrees, burnin = burnin)
  m <- min(vapply(pools, length, integer(1)))
  # interleave: pool[1], pool[2], pool[1], pool[2], ... so both runs
  # contribute from the very first draw
  ord <- vector("list", 2L * m)
  ord[seq(1, 2L * m, by = 2L)] <- pools[[1]][seq_len(m)]
  ord[seq(2, 2L * m, by = 2L)] <- pools[[2]][seq_len(m)]
  poolLen <- length(ord)

  capacity <- (poolLen - 1L) %/% stride + 1L
  if (nTrees > capacity) {
    stop(sprintf(
      "nTrees=%d exceeds scheme capacity=%d - lower `stride` or extend the run.",
      nTrees, capacity))
  }
  idx <- 1L + (seq_len(nTrees) - 1L) * stride
  trees <- ord[idx]
  class(trees) <- "multiPhylo"
  attr(trees, "poolLen") <- poolLen
  attr(trees, "minEss") <- minEss
  # Return: nTrees stable-stride trees, identity of tree_k independent of nTrees
  trees
}
