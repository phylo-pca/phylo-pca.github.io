# 4. Tree inference: neighbour-joining, parsimony, Bayesian

``` r

library(phyloPCA)
library(ape)
data(generative_trees)
truth <- generative_trees[[1]]
set.seed(1)
M <- SimulateContinuous(truth, nChar = 100)
scores <- prcomp(M)$x
```

Smith (2026) reconstructs a tree from every simulated dataset by three
different methods, so that reconstruction accuracy can be attributed to
the *data* or to the *method* used to analyse it. This vignette covers
all three for continuous data; the same functions apply to discrete and
(for parsimony/Bayesian) shape data, with `type`/template arguments
changed as noted.

## Neighbour-joining

Neighbour-joining ([`ape::nj()`](https://rdrr.io/pkg/ape/man/nj.html))
needs no wrapper – this package’s contribution is providing the correct
*distance* to feed it (see
[`vignette("pca-and-alignment")`](https://phylo-pca.github.io/articles/pca-and-alignment.md)):

``` r

treeNJ <- nj(SquaredEuclideanDistance(scores))
ScoreToTruth(treeNJ, truth, metrics = c("CID", "RF"))
#>        CID         RF 
#> 0.04653366 2.00000000
```

## Maximum parsimony (TNT)

[`RunTNTParsimony()`](https://phylo-pca.github.io/reference/RunTNTParsimony.md)
runs the unified search protocol validated across every parsimony arm in
Smith (2026): 50 random-addition-sequence + TBR starting trees, a New
Technology pass (sectorial search, ratchet, drift, tree fusing –
confirms rather than improves on the RAS+TBR optimum here), then
TBR-swaps every optimal tree to completion and saves all
most-parsimonious trees. `type = "continuous"` treats the data as
Wagner-optimized real-valued characters (Goloboff et al. 2006);
`"discrete"` as standard (Fitch) unordered characters; `"landmark"` as
linear-displacement-optimized landmark configurations (Goloboff &
Catalano’s method, following Palci & Lee 2019).

This needs a local TNT installation (`options(phyloPCA.tnt = ...)` or
the `TNT_PATH` environment variable; see
[`FindTNT()`](https://phylo-pca.github.io/reference/FindTNT.md)) – not
available when this vignette was built, so this chunk is shown but not
run:

``` r

result <- RunTNTParsimony(M, type = "continuous")
result$nTree # essentially always 1 for continuous (real-valued) data
ScoreMPTsByMaxEntropy(result$trees, truth)
```

Where a search finds more than one most-parsimonious tree (routine for
discrete and shape data), a strict consensus is over-conservative –
[`ScoreMPTsByMaxEntropy()`](https://phylo-pca.github.io/reference/ScoreMPTsByMaxEntropy.md)
instead takes a maximum-entropy subsample of up to `k` trees (via
[`MaxMin::MaxEntropy()`](https://rdrr.io/pkg/MaxMin/man/MaxEntropy.html),
using clustering information distance as the measure of redundancy
between trees) and reports their mean distance to the truth.

## Bayesian inference (RevBayes)

[`RunRevBayesBM()`](https://phylo-pca.github.io/reference/RunRevBayesBM.md)
and
[`RunRevBayesMk()`](https://phylo-pca.github.io/reference/RunRevBayesMk.md)
run the two-run, Metropolis-coupled Bayesian analyses used throughout
Smith (2026): a Brownian motion model (`dnPhyloBrownianREML`) for
continuous/shape data, and a per-state-size Mk model for discrete data.
Both stop themselves once the Gelman-Rubin PSRF falls below 1.01 and the
parameter ESS exceeds 333.

This needs a local RevBayes installation (see
[`FindRevBayes()`](https://phylo-pca.github.io/reference/FindRevBayes.md))
and takes on the order of a minute even for this small example – not
available when this vignette was built, so this chunk is shown but not
run:

``` r

nexusFile <- tempfile(fileext = ".nex")
WriteContinuousNexus(M, nexusFile)
outputPrefix <- tempfile("bm")
result <- RunRevBayesBM(nexusFile, outputPrefix, seed = 1)
result$mapTrees

mapTree <- ReadMAPTree(result$mapTrees[1])
ScoreToTruth(mapTree, truth, metrics = c("CID", "RF"))
```

[`RunRevBayesMk()`](https://phylo-pca.github.io/reference/RunRevBayesMk.md)
works the same way, from a discrete-character Nexus file
([`WriteDiscreteNexus()`](https://phylo-pca.github.io/reference/WriteDiscreteNexus.md))
instead.

## Putting a number on convergence

Alongside a point estimate (the MAP tree above), Smith (2026) reports a
Gelman-Rubin PSRF and a topological effective sample size for every
Bayesian analysis; see
[`TreeESS()`](https://phylo-pca.github.io/reference/TreeESS.md) in
[`vignette("generative-trees")`](https://phylo-pca.github.io/articles/generative-trees.md)
for the latter. Because RevBayes’ own trace files
(`<prefix>_run_1.trees`, `<prefix>_run_2.trees`) are exactly the input
[`TreeESS()`](https://phylo-pca.github.io/reference/TreeESS.md) expects,
the same diagnostic used to certify the generative-tree posterior
applies unchanged to every inference run in the pipeline.
