# Shape data Bayesian inference: alignment comparison

As `shape_sigma_nj_align`, but for Bayesian inference (Brownian motion;
see
[`RunRevBayesBM()`](https://phylo-pca.github.io/reference/RunRevBayesBM.md))
rather than neighbour-joining; the `"proc"` rows are shared with
`shape_sigma_bayes`.

## Usage

``` r
data(shape_sigma_bayes_align)
```

## Format

A data frame with columns `t` (tree index), `dim`, `sigma`, `align`
(`"proc"`, `"rftra"` or `"raw"`), `Bayes_MAP` (CID of the MAP tree) and
`MAP_CID` (independent re-scoring cross-check).

## Source

`hominin-pca-scale/R/29_shape_bayes_align_prep.R` +
`R/30_shape_bayes_align_score.R`, using `inst/rev/bm_infer_body.Rev` on
Hamilton (`rev/shape_align_array.sh`). Raw/RFTRA coordinates are
rescaled to match Procrustes' total centred sum-of-squares before
inference, so the Brownian tree-length prior is not railroaded by the
native LDDMM coordinate scale (topology-neutral).
