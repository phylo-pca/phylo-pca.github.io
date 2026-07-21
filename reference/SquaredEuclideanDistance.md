# Squared-Euclidean distance

Neighbour-joining assumes an additive input distance. Ordinary Euclidean
distance between PCA scores is not additive, but its *square* is (in
expectation, proportional to patristic distance under Brownian motion):
this is the "additive fix" used throughout Smith (2026) in place of the
raw Euclidean distance used by Raskin et al. (2026).

## Usage

``` r
SquaredEuclideanDistance(X)
```

## Arguments

- X:

  A taxa x variables numeric matrix (e.g. PCA scores).

## Value

A symmetric taxa x taxa distance matrix (a `dist` object squared
element-wise).
