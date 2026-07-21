# Simulate continuous characters under Brownian motion

Thin, documented wrapper around
[`phytools::fastBM()`](https://rdrr.io/pkg/phytools/man/fastBM.html)
matching the simulation used throughout Smith (2026): independent
characters, root state 0, rate `sigma^2 = 1`.

## Usage

``` r
SimulateContinuous(tree, nChar, rate = 1)
```

## Arguments

- tree:

  A phylogenetic tree (generative/"true" tree).

- nChar:

  Number of independent characters to simulate.

- rate:

  Brownian rate (`sigma^2`).

## Value

A tips x `nChar` numeric matrix, row-named by `tree$tip.label`.
