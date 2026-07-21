# Discrete data maximum parsimony

Maximum-parsimony accuracy (TNT, standard Fitch optimization on
unordered multistate characters; see
[`RunTNTParsimony()`](https://phylo-pca.github.io/reference/RunTNTParsimony.md))
on the identical discrete character matrices used in `scale`, scored by
a maximum-entropy subsample of most-parsimonious trees (see
[`ScoreMPTsByMaxEntropy()`](https://phylo-pca.github.io/reference/ScoreMPTsByMaxEntropy.md))
rather than a strict consensus.

## Usage

``` r
scale_disc_parsimony
```

## Format

A data frame with columns `tree`, `nC`, `CID`, `RF`, `SPR` (`NA`; not
computed for parsimony arms, see
[`ScoreToTruth()`](https://phylo-pca.github.io/reference/ScoreToTruth.md))
and `nmpt` (number of most-parsimonious trees found).

## Source

`hominin-pca-scale/R/34_par_hamilton_gen.R` +
`R/35_par_hamilton_score.R`, run on Hamilton under the unified TNT
search protocol.
