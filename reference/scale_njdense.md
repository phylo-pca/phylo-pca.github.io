# Dense neighbour-joining ladder

Neighbour-joining-only accuracy (no Bayesian inference) across a dense
range of character counts, `nC` in `{10, 25, 50, 100, 250, 500}`, for
both continuous and discrete data, several replicate simulations per
(tree, nC) cell.

## Usage

``` r
data(scale_njdense)
```

## Format

A data frame with columns `tree`, `rep` (replicate index), `nC`, `arm`
(`"cont_raw"`, `"cont_sq"`, `"cont_pc12"` or `"disc_cor"`; see `scale`'s
`NJ_*` columns for the equivalent distances) and `CID`.

## Source

`hominin-pca-scale/R/09_scale_njdense.R`.
