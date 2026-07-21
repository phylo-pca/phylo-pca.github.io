# Additive distance for a mixed-state-space discrete matrix

A Jukes-Cantor-style corrected distance for a matrix whose characters
may have different state-space sizes: within each observed state-space
size `k`, the JC-`k` correction is applied to the observed proportion
differing, and expected substitutions are summed across state-size
blocks and expressed per character. Reduces to the standard binary
correction `-0.5 * log(1 - 2p)` when every character has `k = 2` states.
This is the additive distance appropriate for neighbour-joining on
discrete data, mirroring the squared-Euclidean fix used for continuous
data (see
[`SquaredEuclideanDistance()`](https://phylo-pca.github.io/reference/SquaredEuclideanDistance.md)).

## Usage

``` r
CorrectedDistanceMS(X)
```

## Arguments

- X:

  A taxa x characters matrix of state labels (as produced by
  [`SimulateDiscreteMS()`](https://phylo-pca.github.io/reference/SimulateDiscreteMS.md)).

## Value

A symmetric taxa x taxa distance matrix.
