# Converged continuous Bayesian re-run

A 2-chain, longer re-run of the continuous Bayesian inference cells in
`scale`, confirming convergence (all cells ASDSF \< 0.01) at `nC` in
`{25, 50, 100, 250}`. This is the dataset `scale`'s continuous arm and
Figures 1 and 2 actually read the Bayesian numbers from.

## Usage

``` r
data(scale_cont_rerun)
```

## Format

A data frame with columns `tree`, `nC`, `Bayes_MAP` (CID of the MAP
tree), `MAP_CID` (an independent re-scoring, kept for a cross-check),
`Bayes_pm` (posterior-mean CID), `P_true` (posterior probability of the
exact tree) and `rhat` (PSRF).

## Source

`hominin-pca-scale/R/14_rerun_cont.R` + `R/16_add_cont25.R`, using
`inst/rev/bm_infer_body.Rev`.
