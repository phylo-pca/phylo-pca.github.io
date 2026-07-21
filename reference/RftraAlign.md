# Resistant-fit (RFTRA) alignment of landmark configurations

Generalized resistant-fit superimposition (Siegel & Benson 1982; Slice
1996): a repeated-median scale plus Tukey-bisquare iteratively
reweighted Procrustes rotation, so that a small number of grossly
displaced landmarks (the "Pinocchio effect" of Palci & Lee 2019) are
down-weighted rather than smeared across the whole configuration, as
ordinary least-squares Procrustes
([`ProcrustesAlign()`](https://phylo-pca.github.io/reference/ProcrustesAlign.md))
would. Recovers a pure similarity transform exactly, and matches
least-squares Procrustes when no landmark is an outlier.

## Usage

``` r
RftraAlign(configs, iter = 15, tol = 1e-07)
```

## Arguments

- configs:

  A named list of landmarks x dim numeric matrices, one per taxon (as
  returned by
  [`ReadShapeConfigs()`](https://phylo-pca.github.io/reference/ReadShapeConfigs.md)).

- iter:

  Maximum number of consensus-update iterations.

- tol:

  Convergence tolerance on the consensus configuration.

## Value

A named list of aligned landmarks x dim matrices, same names and order
as `configs`.
