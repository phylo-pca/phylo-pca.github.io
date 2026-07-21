# The 250 generative trees used throughout Smith (2026)

Every simulation and reconstruction in Smith (2026) is scored against
one of these 250 trees: a stable-stride sample (see
[`SampleStableTrees()`](https://phylo-pca.github.io/reference/SampleStableTrees.md))
from two independent, converged RevBayes posteriors (see
[`RunRevBayesReinference()`](https://phylo-pca.github.io/reference/RunRevBayesReinference.md))
obtained by re-inferring Mongle et al. (2023)'s discrete character
matrix for 15 hominins and 6 anthropoid outgroups under the Mk model
(Lewis 2001) with the `coding = "variable"` ascertainment-bias
correction. `tree_k` is a fixed function of `k` (raising the tree count
only appends new trees), so trees referenced by index elsewhere in this
package (e.g. `scale$tree`, `shape_combine_nj$t`) always refer to the
same tree.

## Usage

``` r
generative_trees
```

## Format

A `multiPhylo` of 250 trees, named `tree_001`..`tree_250`.

## Source

`hominin-pca-scale/R/20_sample_trees.R`, from
`hominin-pca-scale/out/reinfer/reinfer_run_{1,2}.trees` (see
`inst/rev/reinfer_variable.Rev`).
