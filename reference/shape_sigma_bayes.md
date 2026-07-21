# Shape data Bayesian inference by dimension and sigma

Bayesian inference accuracy (Brownian motion model on Procrustes-aligned
coordinates; see
[`RunRevBayesBM()`](https://phylo-pca.github.io/reference/RunRevBayesBM.md))
for 2D and 3D shape data separately, across the landmark-integration
parameter `sigma`.

## Usage

``` r
shape_sigma_bayes
```

## Format

A data frame with columns `tnum` (tree index), `dim` (2 or 3), `sigma`,
`Bayes_MAP` (CID of the MAP tree), `MAP_CID` (independent re-scoring
cross-check), `Bayes_pm` (posterior-mean CID), `P_true` (posterior
probability of the exact tree) and `rhat` (PSRF).

## Source

`hominin-pca-scale/R/20b_sigma_bayes_score_maps.R`, using
`inst/rev/bm_infer_body.Rev` on Hamilton.
