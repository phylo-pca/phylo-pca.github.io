# Combined 2D+3D shape data Bayesian inference

Bayesian inference accuracy (Brownian motion on the concatenated,
Procrustes-aligned 2D and 3D coordinate matrix) across `sigma`, paired
with the individual-dimension arm in `shape_sigma_bayes` restricted to
the same trees (`unique(shape_both_bayes$t)`).

## Usage

``` r
shape_both_bayes
```

## Format

A data frame with columns `t` (tree index), `sigma`, `Bayes_MAP` (CID of
the MAP tree), `MAP_CID` (independent re-scoring cross-check), `P_true`
and `ASDSF` (average standard deviation of split frequencies, a
convergence diagnostic).

## Source

`hominin-pca-scale/R/24_both_bayes_score.R`, using
`inst/rev/bm_infer_body.Rev` on Hamilton (`rev/both_array.sh`).
