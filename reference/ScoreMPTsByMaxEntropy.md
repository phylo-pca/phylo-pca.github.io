# Score a most-parsimonious tree set by maximum-entropy subsampling

Where a parsimony search finds multiple, equally most-parsimonious
trees, a strict consensus is over-conservative (it penalises a cell for
ambiguity among equally good trees, not for getting the wrong answer).
As in Smith (2026), a density-blind subsample of up to `k` trees is
instead selected to maximize the entropy (spread) within the chosen
subset (via
[`MaxMin::MaxEntropy()`](https://rdrr.io/pkg/MaxMin/man/MaxEntropy.html),
using the clustering information distance as the measure of redundancy
between trees), and their mean distance to the true tree is reported.
When at most `k` trees were found, every tree is used.

## Usage

``` r
ScoreMPTsByMaxEntropy(trees, truth, k = 10L, pool = 200L)
```

## Arguments

- trees:

  A `multiPhylo` of most-parsimonious trees (e.g. from
  [`RunTNTParsimony()`](https://phylo-pca.github.io/reference/RunTNTParsimony.md)).

- truth:

  The true (generative) tree.

- k:

  Maximum subsample size.

- pool:

  If more than `pool` trees are supplied, they are first evenly thinned
  to `pool` trees (most-parsimonious trees are exchangeable, so this
  bounds the O(n^2) distance matrix on pathologically tie-heavy cells
  without materially changing the mean).

## Value

A list with elements `CID` (mean clustering information distance of the
subsample to `truth`), `RF` (mean Robinson-Foulds distance), `nTree`
(the number of most-parsimonious trees supplied) and `nSel` (the
subsample size actually used).
