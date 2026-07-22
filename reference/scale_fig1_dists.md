# Figure 1 method-comparison distances (continuous data)

The five method arms compared in Figure 1 of Smith (2026), all on the
identical continuous character matrices used in `scale`: the
PCA-\\neighbour-joining pipeline, all-dimension and additive-distance
neighbour-joining, maximum parsimony (TNT, continuous Wagner; see
[`RunTNTParsimony()`](https://phylo-pca.github.io/reference/RunTNTParsimony.md))
and Bayesian inference, each scored by clustering information distance,
Robinson-Foulds distance and (where both trees are binary) SPR distance.

## Usage

``` r
scale_fig1_dists
```

## Format

A data frame with columns `tree`, `nC`, `arm` (`"pipeline"`, `"allpc"`,
`"additive"`, `"mp"` or `"bayes"`), `CID`, `RF`, `SPR` and `SPR_exact`
(whether the SPR distance could be computed exactly; see
[`ScoreToTruth()`](https://phylo-pca.github.io/reference/ScoreToTruth.md)).

## Source

`hominin-pca-scale/R/15_fig1_rfspr.R` (pipeline/allpc/ additive/bayes
arms), `R/34_par_hamilton_gen.R` + `R/35_par_hamilton_score.R` (the
`"mp"` arm, run on Hamilton under the unified TNT search protocol; see
[`RunTNTParsimony()`](https://phylo-pca.github.io/reference/RunTNTParsimony.md)).

## Note

The `"mp"` arm's `SPR` was recomputed on Hamilton via
[`TBRDist::USPRDist()`](https://ms609.github.io/TBRDist/reference/TreeRearrangementDistances.html)
under a per-cell subprocess timeout (999/1000 cells;
`SPR_exact = TRUE`). The one cell that did not resolve within a 6-hour
budget (a genuinely hard NP-hard instance, not a bug: two markedly
dissimilar 21-tip trees, RF = 20) is filled from
[`TreeDist::SPRDist()`](https://ms609.github.io/TreeDist/reference/SPRDist.html)'s
polynomial-time approximation instead (`SPR_exact = FALSE`), matching
the convention already used for the other arms' occasional approximate
fallbacks.
