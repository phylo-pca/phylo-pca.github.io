# Score a reconstructed tree against the true tree

Computes one or more tree-distance metrics between a reconstructed tree
and the generative ("true") tree it should be compared to. Following
Smith (2026), the clustering information distance (Smith 2020) is the
primary, robust metric (normalized to 0-1; the expected distance between
random 21-tip topologies is ~0.84); Robinson-Foulds distance is reported
for comparability with Raskin et al. (2026). The unrooted
subtree-prune-and-regraft (SPR) distance is NP-hard to compute exactly
and its reference implementation
([`TBRDist::USPRDist()`](https://ms609.github.io/TBRDist/reference/TreeRearrangementDistances.html))
can hang indefinitely and is uninterruptible from R on a hard tree pair;
it is therefore only attempted when both trees are strictly binary, and
`NA` is returned for non-binary input (e.g. a maximum-parsimony
consensus) rather than risking a hang.

## Usage

``` r
ScoreToTruth(tree, truth, metrics = c("CID", "RF", "SPR"))
```

## Arguments

- tree:

  The reconstructed tree.

- truth:

  The true (generative) tree.

- metrics:

  Character vector of metrics to compute: any of `"CID"` (clustering
  information distance), `"RF"` (Robinson-Foulds distance) and `"SPR"`
  (subtree prune-and-regraft distance).

## Value

A named numeric vector with one element per requested metric.
