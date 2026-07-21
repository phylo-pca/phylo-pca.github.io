# Re-infer a posterior of generative trees

Runs the RevBayes re-inference used throughout Smith (2026) to obtain a
posterior of plausible generative trees for simulation: the Mk model
(Lewis 2001), matching Raskin et al. (2026)'s own model exactly except
for the ascertainment-bias correction (`coding = "variable"`, the
correct choice for a matrix with no invariant characters, such as Mongle
et al. (2023)'s; see
[`MongleStateProportions()`](https://phylo-pca.github.io/reference/MongleStateProportions.md)).
Two independent runs are written so that topological convergence can
subsequently be assessed with
[`TreeESS()`](https://phylo-pca.github.io/reference/TreeESS.md). This
step is CPU-intensive (hours, even for 21 taxa) - `generative_trees`
provides the actual 250-tree sample used in Smith (2026), sparing most
users from re-running it.

## Usage

``` r
RunRevBayesReinference(
  nexusFile,
  outDir,
  rbPath = FindRevBayes(),
  nGen = 250000,
  printGen = 50
)
```

## Arguments

- nexusFile:

  Path to a discrete character Nexus file.

- outDir:

  Directory for RevBayes' logs and tree traces (`reinfer.log`,
  `reinfer_run_1.trees`, `reinfer_run_2.trees`).

- rbPath:

  Path to the RevBayes executable; see
  [`FindRevBayes()`](https://phylo-pca.github.io/reference/FindRevBayes.md).

- nGen:

  Number of MCMC generations per run.

- printGen:

  Sampling frequency, in generations.

## Value

A list with elements `runFiles` (the two `.trees` files produced,
suitable input to
[`SampleStableTrees()`](https://phylo-pca.github.io/reference/SampleStableTrees.md))
and `log`.
