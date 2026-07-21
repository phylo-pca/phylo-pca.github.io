# Sample a stable set of generative trees from a posterior

Draws a fixed, append-only set of "truth" trees from two independent
post-burnin RevBayes tree traces, for use as the generative tree in
simulation. Runs are interleaved (so both chains contribute from the
first tree onward) and then strided, so that `tree_k` is a fixed
function of `k`: raising `nTrees` later only appends new trees and never
changes the identity of an existing one.

## Usage

``` r
SampleStableTrees(
  runFiles,
  nTrees,
  stride = 11L,
  burnin = 0.25,
  threshold = 200,
  checkEss = TRUE
)
```

## Arguments

- runFiles:

  Character vector of paths to RevBayes `.trees` files, one per
  independent run.

- nTrees:

  Number of trees to sample.

- stride:

  Spacing (in interleaved-pool samples) between sampled trees; should be
  odd, and large enough that consecutive draws are effectively
  independent (checked via
  [`TreeESS()`](https://phylo-pca.github.io/reference/TreeESS.md) when
  `checkEss = TRUE`).

- burnin:

  Proportion of samples discarded from each trace.

- threshold:

  Minimum summed ESS (across runs) for convergence.

- checkEss:

  Whether to require the input traces to have converged (summed
  topological ESS \>= `threshold`) before sampling.

## Value

A `multiPhylo` of `nTrees` sampled trees.
