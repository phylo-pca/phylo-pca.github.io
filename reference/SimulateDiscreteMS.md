# Simulate multistate discrete characters under the Mk model

Simulates `nChar` unlinked discrete characters, each drawing its
state-space size `k` from `stateProportions` (see
[`MongleStateProportions()`](https://phylo-pca.github.io/reference/MongleStateProportions.md))
and evolving under a symmetric `k`-state Markov model (equal rates, flat
root frequencies) with no rate variation among characters. Characters
that evolve with no variation across tips are rejected and re-drawn
(`coding = "variable"`), matching the ascertainment correction used at
inference time.

## Usage

``` r
SimulateDiscreteMS(tree, nChar, stateProportions, rate = 1)
```

## Arguments

- tree:

  A phylogenetic tree (generative/"true" tree).

- nChar:

  Number of characters to simulate.

- stateProportions:

  A named vector/table of relative frequencies at which each state-space
  size `k` should be drawn (as returned by
  [`MongleStateProportions()`](https://phylo-pca.github.io/reference/MongleStateProportions.md)).

- rate:

  Overall character rate (branch lengths are used as-is; this is a
  multiplier on them).

## Value

A tips x `nChar` character matrix of state labels `"0".."k-1"`,
row-named by `tree$tip.label`, every column variable across tips.
