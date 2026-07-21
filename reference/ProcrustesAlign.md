# Generalized Procrustes alignment of landmark configurations

A minimal, dependency-free generalized Procrustes analysis (GPA):
centres each configuration on its centroid, scales it to unit centroid
size, then iteratively rotates every configuration (reflections
forbidden) onto the evolving consensus until convergence. Equivalent to
`geomorph::gpagen()` (Adams & Otarola-Castillo 2013) for landmark data
with no missing values or sliding semilandmarks, used throughout Smith
(2026) in preference to `geomorph::gpagen()` so that the exact
rotation/scaling steps are visible and independently checkable.

## Usage

``` r
ProcrustesAlign(configs, tol = 1e-07, maxit = 200)
```

## Arguments

- configs:

  A named list of landmarks x dim numeric matrices, one per taxon (as
  returned by
  [`ReadShapeConfigs()`](https://phylo-pca.github.io/reference/ReadShapeConfigs.md)).

- tol:

  Convergence tolerance on the consensus configuration.

- maxit:

  Maximum number of GPA iterations.

## Value

A named list of aligned landmarks x dim matrices, same names and order
as `configs`.
