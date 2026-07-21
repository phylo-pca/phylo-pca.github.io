# Topological effective sample size across RevBayes runs

Estimates the effective sample size of tree topology across one or more
independent RevBayes runs, using the lower of the Frechet correlation
and median pseudo-ESS statistics (the pair used throughout Smith 2026),
as recommended by Fabreti & Hohna (2022) and Magee et al. (2024).

## Usage

``` r
TreeESS(runFiles, burnin = 0.25, cap = 1000, threshold = 200)
```

## Arguments

- runFiles:

  Character vector of paths to RevBayes `.trees` files, one per
  independent run.

- burnin:

  Proportion of samples discarded from each trace.

- cap:

  Maximum number of trees retained per run (traces longer than this are
  evenly thinned before the ESS calculation, which scales poorly with
  sample size).

- threshold:

  Minimum summed ESS (across runs) for convergence.

## Value

A list with elements `ok` (logical, whether `threshold` is met),
`perRun` (a data frame of per-run ESS estimates), `summed` (the two ESS
statistics summed across runs) and `minEss` (the smaller of the two).
