# Bayesian inference of a tree from continuous character data

Runs the two-run, Metropolis-coupled Bayesian phylogenetic analysis of
continuous character data (a Brownian motion model,
`dnPhyloBrownianREML`) used throughout Smith (2026), for both simulated
continuous morphological characters and Procrustes/RFTRA-aligned shape
coordinates (both are continuous data as far as the model is concerned).
Each run stops itself once the Gelman-Rubin PSRF falls below 1.01 and
the parameter ESS exceeds 333 (or `maxHours` elapses).

## Usage

``` r
RunRevBayesBM(
  nexusFile,
  outputPrefix,
  rbPath = FindRevBayes(),
  maxHours = 2,
  seed = NULL
)
```

## Arguments

- nexusFile:

  Path to a continuous-character Nexus file (see
  [`WriteContinuousNexus()`](https://phylo-pca.github.io/reference/WriteContinuousNexus.md)).

- outputPrefix:

  File prefix for RevBayes' outputs (`<prefix>.trees`, `<prefix>.p.log`,
  `<prefix>_map1.tre`, `<prefix>_map2.tre`, ...).

- rbPath:

  Path to the RevBayes executable; see
  [`FindRevBayes()`](https://phylo-pca.github.io/reference/FindRevBayes.md).

- maxHours:

  Wall-clock stopping rule, in hours.

- seed:

  Optional integer random seed, for a reproducible run (RevBayes uses
  its own internal default seeding if omitted).

## Value

A list with elements `mapTrees` (paths to the per-run maximum *a
posteriori* tree files that were produced) and `log` (the RevBayes
console log file).
