# Optimality-based tree recovery

Whether the generative tree is (tied for) most parsimonious, following
the recovery definition used for the parsimony arms in Smith (2026): if
a heuristic search misses the true topology, that is the search
heuristic's fault, not a failure of parsimony as a criterion, so
recovery is judged by comparing the true tree's own parsimony score to
the best score found, rather than by topological identity with a
most-parsimonious tree.

## Usage

``` r
OptimalityRecovery(trueLength, bestLength, tolerance = 1e-04)
```

## Arguments

- trueLength:

  Parsimony score of the generative tree.

- bestLength:

  Best (minimum) parsimony score found by the search.

- tolerance:

  Numeric tolerance for the comparison (guards against floating-point
  score differences for continuous/landmark data).

## Value

A logical: whether the generative tree is tied for most parsimonious.
