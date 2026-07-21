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
