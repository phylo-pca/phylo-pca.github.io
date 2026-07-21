#' Read landmark configurations from a shape-simulator TSV
#'
#' Reads the long-format landmark table written by the LDDMM shape simulator
#' of Raskin et al. (2026) (one row per landmark per tip, columns `label`,
#' an index, then one column per spatial dimension) and reshapes it into one
#' matrix per tip.
#'
#' @param file Path to a `*_nodeShapes.tsv` file.
#' @param dim Number of spatial dimensions (2 or 3).
#'
#' @return A named list (names = tip labels) of landmarks x `dim` numeric
#'   matrices.
#'
#' @export
ReadShapeConfigs <- function(file, dim) {
  d <- utils::read.delim(file, header = FALSE)
  labels <- unique(d$V1)
  # Return: one landmarks x dim matrix per tip label
  stats::setNames(
    lapply(labels, function(lab) as.matrix(d[d$V1 == lab, 3:(2 + dim), drop = FALSE])),
    labels
  )
}

.CentreScale <- function(X) {
  Xc <- sweep(X, 2, colMeans(X))
  # Return: centred configuration scaled to unit centroid size
  Xc / sqrt(sum(Xc^2))
}

.ProcRotate <- function(A, B) {
  M <- t(A) %*% B
  sv <- svd(M)
  R <- sv$u %*% t(sv$v)
  if (det(R) < 0) { # forbid reflection
    sv$u[, ncol(sv$u)] <- -sv$u[, ncol(sv$u)]
    R <- sv$u %*% t(sv$v)
  }
  # Return: A rotated onto B
  A %*% R
}

#' Generalized Procrustes alignment of landmark configurations
#'
#' A minimal, dependency-free generalized Procrustes analysis (GPA):
#' centres each configuration on its centroid, scales it to unit centroid
#' size, then iteratively rotates every configuration (reflections
#' forbidden) onto the evolving consensus until convergence. Equivalent to
#' `geomorph::gpagen()` (Adams & Otarola-Castillo 2013) for landmark data
#' with no missing values or sliding semilandmarks, used throughout
#' Smith (2026) in preference to `geomorph::gpagen()` so that the exact
#' rotation/scaling steps are visible and independently checkable.
#'
#' @param configs A named list of landmarks x dim numeric matrices, one per
#'   taxon (as returned by [ReadShapeConfigs()]).
#' @param tol Convergence tolerance on the consensus configuration.
#' @param maxit Maximum number of GPA iterations.
#'
#' @return A named list of aligned landmarks x dim matrices, same names and
#'   order as `configs`.
#'
#' @export
ProcrustesAlign <- function(configs, tol = 1e-7, maxit = 200) {
  Z <- lapply(configs, .CentreScale)
  ref <- Z[[1]]
  for (it in seq_len(maxit)) {
    Z <- lapply(Z, .ProcRotate, B = ref)
    consensus <- .CentreScale(Reduce(`+`, Z) / length(Z))
    if (max(abs(consensus - ref)) < tol) {
      ref <- consensus
      break
    }
    ref <- consensus
  }
  # Return: Procrustes-aligned configurations, names preserved
  stats::setNames(Z, names(configs))
}

.RhoResistantScale <- function(X, Y) {
  # Repeated-median interlandmark-distance ratio: a resistant (outlier-safe)
  # scale factor between two configurations, valid in any dimension.
  k <- nrow(X)
  perLandmark <- rep(NA_real_, k)
  for (i in seq_len(k)) {
    dX <- sqrt(rowSums(sweep(X[-i, , drop = FALSE], 2, X[i, ])^2))
    dY <- sqrt(rowSums(sweep(Y[-i, , drop = FALSE], 2, Y[i, ])^2))
    ok <- dY > 1e-9 # skip (near-)coincident landmarks: avoid 1/0
    if (any(ok)) perLandmark[i] <- stats::median(dX[ok] / dY[ok])
  }
  r <- stats::median(perLandmark, na.rm = TRUE)
  # Return: resistant scale factor (falls back to 1 if degenerate)
  if (!is.finite(r) || r <= 0) 1 else r
}

.WeightedProcrustes <- function(X, Y, w) {
  W <- w / sum(w)
  mx <- colSums(W * X)
  my <- colSums(W * Y)
  Xc <- sweep(X, 2, mx)
  Yc <- sweep(Y, 2, my)
  s <- svd(t(Yc * W) %*% Xc) # M = Yc'X; R = UV' aligns Yc onto Xc
  dsign <- sign(det(s$u %*% t(s$v))) # forbid reflection
  R <- s$u %*% diag(c(rep(1, ncol(X) - 1L), dsign), ncol(X)) %*% t(s$v)
  # Return: Y rotated+translated onto X
  sweep(Yc %*% R, 2, -mx)
}

.GrfFit <- function(X, Y, iter = 40, tol = 1e-9, tuning = 4.685) {
  Y <- .RhoResistantScale(X, Y) * Y
  w <- rep(1, nrow(X))
  previous <- NULL
  fitted <- Y
  for (it in seq_len(iter)) {
    fitted <- .WeightedProcrustes(X, Y, w)
    residual <- sqrt(rowSums((fitted - X)^2))
    scale <- stats::median(residual)
    if (scale < 1e-10) {
      # Return: (near-)exact fit
      return(fitted)
    }
    if (!is.null(previous) && max(abs(fitted - previous)) < tol) break
    previous <- fitted
    u <- residual / (tuning * scale)
    w <- ifelse(u < 1, (1 - u^2)^2, 0) # Tukey bisquare weights
    if (all(w == 0)) w <- rep(1, length(w))
  }
  # Return: resistant-fit alignment of Y onto X
  fitted
}

#' Resistant-fit (RFTRA) alignment of landmark configurations
#'
#' Generalized resistant-fit superimposition (Siegel & Benson 1982; Slice
#' 1996): a repeated-median scale plus Tukey-bisquare iteratively reweighted
#' Procrustes rotation, so that a small number of grossly displaced
#' landmarks (the "Pinocchio effect" of Palci & Lee 2019) are down-weighted
#' rather than smeared across the whole configuration, as ordinary
#' least-squares Procrustes ([ProcrustesAlign()]) would. Recovers a pure
#' similarity transform exactly, and matches least-squares Procrustes when
#' no landmark is an outlier.
#'
#' @inheritParams ProcrustesAlign
#' @param iter Maximum number of consensus-update iterations.
#' @param tol Convergence tolerance on the consensus configuration.
#'
#' @return A named list of aligned landmarks x dim matrices, same names and
#'   order as `configs`.
#'
#' @export
RftraAlign <- function(configs, iter = 15, tol = 1e-7) {
  alignTo <- function(ref) lapply(configs, function(Y) .GrfFit(ref, Y))
  ref <- configs[[1]]
  ref <- sweep(ref, 2, colMeans(ref))
  for (it in seq_len(iter)) {
    aligned <- alignTo(ref)
    arr <- simplify2array(aligned) # landmarks x dim x taxa
    consensus <- apply(arr, c(1, 2), stats::median) # resistant consensus
    consensus <- sweep(consensus, 2, colMeans(consensus))
    if (max(abs(consensus - ref)) < tol) {
      ref <- consensus
      break
    }
    ref <- consensus
  }
  # Return: RFTRA-aligned configurations, names preserved
  stats::setNames(alignTo(ref), names(configs))
}
