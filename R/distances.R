#' Hamming (p-)distance between discrete character rows
#'
#' Proportion of characters differing between each pair of taxa; valid for
#' any state-space size (unlike a binary Hamming distance, no assumption is
#' made about the number of states a character can take).
#'
#' @param X A taxa x characters matrix of state labels.
#'
#' @return A symmetric taxa x taxa distance matrix.
#'
#' @export
HammingDistance <- function(X) {
  n <- nrow(X)
  D <- matrix(0, n, n, dimnames = list(rownames(X), rownames(X)))
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      D[i, j] <- D[j, i] <- mean(X[i, ] != X[j, ])
    }
  }
  # Return: pairwise proportion-differing distance matrix
  D
}

#' Additive distance for a mixed-state-space discrete matrix
#'
#' A Jukes-Cantor-style corrected distance for a matrix whose characters may
#' have different state-space sizes: within each observed state-space size
#' `k`, the JC-`k` correction is applied to the observed proportion
#' differing, and expected substitutions are summed across state-size blocks
#' and expressed per character. Reduces to the standard binary correction
#' `-0.5 * log(1 - 2p)` when every character has `k = 2` states. This is the
#' additive distance appropriate for neighbour-joining on discrete data,
#' mirroring the squared-Euclidean fix used for continuous data (see
#' [SquaredEuclideanDistance()]).
#'
#' @param X A taxa x characters matrix of state labels (as produced by
#'   [SimulateDiscreteMS()]).
#'
#' @return A symmetric taxa x taxa distance matrix.
#'
#' @export
CorrectedDistanceMS <- function(X) {
  stateSize <- apply(X, 2, function(col) max(as.integer(col)) + 1L)
  stateSize[stateSize < 2] <- 2L
  n <- nrow(X)
  nChar <- ncol(X)
  D <- matrix(0, n, n, dimnames = list(rownames(X), rownames(X)))
  uniqueK <- sort(unique(stateSize))
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      diff <- X[i, ] != X[j, ]
      total <- 0
      for (k in uniqueK) {
        sel <- stateSize == k
        nSel <- sum(sel)
        if (!nSel) next
        pMax <- (k - 1) / k
        p <- min(mean(diff[sel]), pMax - 1e-6) # clamp below saturation
        total <- total + nSel * (-pMax * log(1 - p / pMax)) # JC-k expected subs
      }
      D[i, j] <- D[j, i] <- total / nChar
    }
  }
  # Return: additive (JC-corrected) distance matrix
  D
}

#' Squared-Euclidean distance
#'
#' Neighbour-joining assumes an additive input distance. Ordinary Euclidean
#' distance between PCA scores is not additive, but its *square* is (in
#' expectation, proportional to patristic distance under Brownian motion):
#' this is the "additive fix" used throughout Smith (2026) in place of the
#' raw Euclidean distance used by Raskin et al. (2026).
#'
#' @param X A taxa x variables numeric matrix (e.g. PCA scores).
#'
#' @return A symmetric taxa x taxa distance matrix (a `dist` object squared
#'   element-wise).
#'
#' @export
SquaredEuclideanDistance <- function(X) {
  # Return: element-wise square of the Euclidean distance matrix
  stats::dist(X)^2
}

#' Write a continuous-character Nexus file
#'
#' Writes a `datatype=Continuous` Nexus matrix readable by RevBayes'
#' `readContinuousCharacterData()`.
#'
#' @param M A taxa x characters numeric matrix, row-named by taxon.
#' @param file Output file path.
#'
#' @return `file`, invisibly.
#'
#' @export
WriteContinuousNexus <- function(M, file) {
  con <- file(file, "w")
  on.exit(close(con))
  cat("#NEXUS\n\nBegin data;\n", file = con)
  cat(sprintf("Dimensions ntax=%d nchar=%d;\n", nrow(M), ncol(M)), file = con)
  cat("Format datatype=Continuous missing=?;\nMatrix\n", file = con)
  for (i in seq_len(nrow(M))) {
    cat(sprintf("%s  %s\n", rownames(M)[i],
                paste(formatC(M[i, ], format = "f", digits = 6), collapse = " ")),
        file = con)
  }
  cat(";\nEnd;\n", file = con)
  invisible(file)
}

#' Write a discrete-character Nexus file
#'
#' Writes a `datatype=Standard` Nexus matrix (one symbol per cell) readable
#' by RevBayes' `readDiscreteCharacterData()`.
#'
#' @param X A taxa x characters matrix of single-digit state labels (as
#'   produced by [SimulateDiscreteMS()]).
#' @param file Output file path.
#'
#' @return `file`, invisibly.
#'
#' @export
WriteDiscreteNexus <- function(X, file) {
  maxSymbol <- max(as.integer(X))
  symbols <- paste(0:maxSymbol, collapse = "")
  con <- file(file, "w")
  on.exit(close(con))
  cat("#NEXUS\n\nBegin data;\n", file = con)
  cat(sprintf("Dimensions ntax=%d nchar=%d;\n", nrow(X), ncol(X)), file = con)
  cat(sprintf('Format datatype=Standard missing=? gap=- symbols="%s";\nMatrix\n', symbols),
      file = con)
  for (i in seq_len(nrow(X))) {
    cat(sprintf("%s  %s\n", rownames(X)[i], paste(X[i, ], collapse = "")), file = con)
  }
  cat(";\nEnd;\n", file = con)
  invisible(file)
}
