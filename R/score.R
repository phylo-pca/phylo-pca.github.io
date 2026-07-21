#' Score a reconstructed tree against the true tree
#'
#' Computes one or more tree-distance metrics between a reconstructed tree
#' and the generative ("true") tree it should be compared to. Following
#' Smith (2026), the clustering information distance (Smith 2020) is the
#' primary, robust metric (normalized to 0-1; the expected distance between
#' random 21-tip topologies is ~0.84); Robinson-Foulds distance is reported
#' for comparability with Raskin et al. (2026). The unrooted
#' subtree-prune-and-regraft (SPR) distance is NP-hard to compute exactly
#' and its reference implementation (`TBRDist::USPRDist()`) can hang
#' indefinitely and is uninterruptible from R on a hard tree pair; it is
#' therefore only attempted when both trees are strictly binary, and `NA`
#' is returned for non-binary input (e.g. a maximum-parsimony consensus)
#' rather than risking a hang.
#'
#' @param tree The reconstructed tree.
#' @param truth The true (generative) tree.
#' @param metrics Character vector of metrics to compute: any of `"CID"`
#'   (clustering information distance), `"RF"` (Robinson-Foulds distance)
#'   and `"SPR"` (subtree prune-and-regraft distance).
#'
#' @return A named numeric vector with one element per requested metric.
#'
#' @importFrom ape unroot is.binary
#' @importFrom TreeDist ClusteringInfoDistance RobinsonFoulds
#' @export
ScoreToTruth <- function(tree, truth, metrics = c("CID", "RF", "SPR")) {
  metrics <- match.arg(metrics, several.ok = TRUE)
  tree <- unroot(tree)
  truth <- unroot(truth)
  out <- stats::setNames(rep(NA_real_, length(metrics)), metrics)
  if ("CID" %in% metrics) {
    out["CID"] <- ClusteringInfoDistance(tree, truth, normalize = TRUE)
  }
  if ("RF" %in% metrics) {
    out["RF"] <- as.numeric(RobinsonFoulds(tree, truth, normalize = FALSE))
  }
  if ("SPR" %in% metrics) {
    if (is.binary(tree) && is.binary(truth)) {
      if (!requireNamespace("TBRDist", quietly = TRUE)) {
        stop("Package 'TBRDist' is required for the SPR distance.")
      }
      out["SPR"] <- as.numeric(TBRDist::USPRDist(tree, truth))
    }
    # else: leave NA -- non-binary input is not evaluated, to avoid a hang
  }
  # Return: requested distance metrics, SPR NA on non-binary input
  out
}

#' Exact-match tree recovery
#'
#' Whether a reconstructed tree exactly matches the true generative tree
#' (i.e. a clustering information distance of exactly zero). Reported
#' throughout Smith (2026) alongside continuous accuracy metrics, but
#' interpreted with care: for 21-tip unrooted binary trees there are
#' `(2 x 21 - 5)!! ~= 8.2 x 10^21` possible topologies, so exact-match is a
#' vanishingly strict, non-metric criterion, not a graded measure of
#' accuracy.
#'
#' @param tree The reconstructed tree.
#' @param truth The true (generative) tree.
#'
#' @return A logical.
#'
#' @importFrom TreeDist ClusteringInfoDistance
#' @importFrom ape unroot
#' @export
ExactlyRecovered <- function(tree, truth) {
  # Return: TRUE iff tree and truth have identical clustering information
  ClusteringInfoDistance(unroot(tree), unroot(truth), normalize = TRUE) == 0
}
