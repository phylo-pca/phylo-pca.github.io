# Read a RevBayes maximum *a posteriori* tree

A thin wrapper around
[`ape::read.nexus()`](https://rdrr.io/pkg/ape/man/read.nexus.html) for
the annotated `_map*.tre` files written by RevBayes' `mapTree()` (used
as the single representative tree for each Bayesian analysis in Smith
2026, since it is the topology most sampled in the posterior, and thus
has the highest posterior probability).

## Usage

``` r
ReadMAPTree(file)
```

## Arguments

- file:

  Path to a RevBayes MAP tree file.

## Value

A `phylo` object.
