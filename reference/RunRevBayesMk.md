# Bayesian inference of a tree from discrete character data

Runs the two-run, Metropolis-coupled Bayesian phylogenetic analysis of
multistate discrete character data used throughout Smith (2026): a
per-state-size Mk model (Lewis 2001; `fnJC(k)` for each observed
state-space size `k`, matching
[`SimulateDiscreteMS()`](https://phylo-pca.github.io/reference/SimulateDiscreteMS.md)),
conditioned on variable characters only (`coding = "variable"`). Each
run stops itself once the Gelman-Rubin PSRF falls below 1.01 and the
parameter ESS exceeds 333 (or `maxHours` elapses).

## Usage

``` r
RunRevBayesMk(
  nexusFile,
  outputPrefix,
  rbPath = FindRevBayes(),
  maxHours = 2,
  seed = NULL
)
```

## Arguments

- nexusFile:

  Path to a discrete-character Nexus file (see
  [`WriteDiscreteNexus()`](https://phylo-pca.github.io/reference/WriteDiscreteNexus.md)).

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
