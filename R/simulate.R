#' State-space size distribution of the Mongle (2023) character matrix
#'
#' Reads a Nexus character matrix and tabulates the size of each character's
#' state space (`max(observed state) + 1`, following RevBayes'
#' `setNumStatesPartition()` convention; polymorphisms such as `1/2` are
#' split before taking the maximum). Used to draw realistic state-space
#' sizes for [SimulateDiscreteMS()].
#'
#' @param file Path to a Nexus-formatted discrete character matrix.
#'
#' @return A named integer vector (a `table`): counts of characters at each
#'   observed state-space size `k`.
#'
#' @importFrom ape read.nexus.data
#' @export
MongleStateProportions <- function(file) {
  d <- ape::read.nexus.data(file)
  m <- do.call(rbind, d)
  stateCount <- function(column) {
    tokens <- unlist(strsplit(column[!column %in% c("?", "-")], "/"))
    values <- suppressWarnings(as.integer(tokens))
    values <- values[!is.na(values)]
    if (!length(values)) {
      # Return: character carries no scoreable states (all missing)
      return(NA_integer_)
    }
    # Return: observed state-space size for this character
    max(values) + 1L
  }
  k <- apply(m, 2, stateCount)
  k <- k[!is.na(k)]
  # Return: character counts by state-space size
  table(factor(k, levels = 2:max(k)))
}

#' Simulate continuous characters under Brownian motion
#'
#' Thin, documented wrapper around [phytools::fastBM()] matching the
#' simulation used throughout Smith (2026): independent characters, root
#' state 0, rate `sigma^2 = 1`.
#'
#' @param tree A phylogenetic tree (generative/"true" tree).
#' @param nChar Number of independent characters to simulate.
#' @param rate Brownian rate (`sigma^2`).
#'
#' @return A tips x `nChar` numeric matrix, row-named by `tree$tip.label`.
#'
#' @importFrom phytools fastBM
#' @export
SimulateContinuous <- function(tree, nChar, rate = 1) {
  M <- fastBM(tree, sig2 = rate, nsim = nChar)
  rownames(M) <- tree$tip.label
  # Return: tips x nChar continuous character matrix
  M
}

#' Simulate multistate discrete characters under the Mk model
#'
#' Simulates `nChar` unlinked discrete characters, each drawing its
#' state-space size `k` from `stateProportions` (see
#' [MongleStateProportions()]) and evolving under a symmetric `k`-state
#' Markov model (equal rates, flat root frequencies) with no rate variation
#' among characters. Characters that evolve with no variation across tips
#' are rejected and re-drawn (`coding = "variable"`), matching the
#' ascertainment correction used at inference time.
#'
#' @param tree A phylogenetic tree (generative/"true" tree).
#' @param nChar Number of characters to simulate.
#' @param stateProportions A named vector/table of relative frequencies at
#'   which each state-space size `k` should be drawn (as returned by
#'   [MongleStateProportions()]).
#' @param rate Overall character rate (branch lengths are used as-is; this
#'   is a multiplier on them).
#'
#' @return A tips x `nChar` character matrix of state labels `"0".."k-1"`,
#'   row-named by `tree$tip.label`, every column variable across tips.
#'
#' @importFrom TreeTools NTip
#' @importFrom phangorn simSeq
#' @export
SimulateDiscreteMS <- function(tree, nChar, stateProportions, rate = 1) {
  ks <- as.integer(names(stateProportions))
  p <- as.numeric(stateProportions) / sum(stateProportions)
  draw <- as.vector(stats::rmultinom(1, nChar, p))
  blocks <- list()
  for (j in seq_along(ks)) {
    need <- draw[j]
    if (need == 0) next
    k <- ks[j]
    levels <- as.character(0:(k - 1))
    Q <- rep(1, k * (k - 1) / 2)
    bf <- rep(1 / k, k)
    got <- matrix(nrow = NTip(tree), ncol = 0)
    guard <- 0
    while (ncol(got) < need && guard < 500) {
      block <- as.character(simSeq(
        tree, l = max(need * 3L, 100L), type = "USER",
        levels = levels, Q = Q, bf = bf, rate = rate
      ))
      variable <- apply(block, 2, function(col) length(unique(col)) > 1)
      got <- cbind(got, block[, variable, drop = FALSE])
      guard <- guard + 1
    }
    if (ncol(got) < need) {
      stop(sprintf("rejection stalled: k=%d need=%d got=%d", k, need, ncol(got)))
    }
    blocks[[length(blocks) + 1]] <- got[, seq_len(need), drop = FALSE]
  }
  X <- do.call(cbind, blocks)
  # Return: tips x nChar discrete character matrix in canonical tip order
  X[tree$tip.label, , drop = FALSE]
}
