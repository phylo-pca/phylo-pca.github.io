# phyloPCA

<!-- badges: start -->
[![R-CMD-check](https://github.com/phylo-pca/phylo-pca.github.io/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/phylo-pca/phylo-pca.github.io/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Simulation, alignment, tree-inference and scoring functions -- and the
cached results they produced -- behind:

> Smith, M. R. (2026). Recovering phylogenetic structure with principal
> components analysis. A reply to Raskin, Šešelj, Bitarello, Stroustrup, Li
> & Huelsenbeck (2026), *American Journal of Biological Anthropology*.

Raskin et al. (2026) show that neighbour-joining trees built from PCA
scores of simulated continuous morphological data almost never match the
tree the data were generated on. This package provides everything needed
to reproduce our reply's finding that this is a property of *how the data
were analysed*, not of continuous data or PCA in general: an unrooted
21-tip tree topology is drawn from a vanishingly large space (making exact
recovery a poor accuracy measure), neighbour-joining assumes an *additive*
input distance (which raw Euclidean PCA distance is not, but its square
is), and once these are corrected, phylogenetic reconstruction from
continuous, discrete and geometric-morphometric ("shape") data is broadly
comparable in accuracy -- whether reconstructed by neighbour-joining,
maximum parsimony, or Bayesian inference.

## Installation

```r
# install.packages("remotes")
remotes::install_github("phylo-pca/phylo-pca.github.io")
```

Simulating and scoring data needs only R and its dependencies (all on
CRAN, plus [MaxMin](https://github.com/ms609/MaxMin) and
[treess](https://github.com/ms609/treess), installed automatically as
`Remotes`). Reproducing the maximum-parsimony and Bayesian-inference arms
from scratch additionally needs
[TNT](http://www.lillo.org.ar/phylogeny/tnt/) and
[RevBayes](https://revbayes.github.io/) installed locally -- see
`vignette("tree-inference")`. Regenerating the geometric-morphometric
shape data needs Raskin et al. (2026)'s own LDDMM shape simulator, which
this package does not vendor (their repository carries no license) -- see
`vignette("simulating-data")`.

## What's here

- **Simulation**: continuous characters under Brownian motion
  (`SimulateContinuous()`), multistate discrete characters under the Mk
  model (`SimulateDiscreteMS()`), and readers for the externally-simulated
  geometric-morphometric shape data (`ReadShapeConfigs()`).
- **Alignment**: dependency-free generalized Procrustes analysis
  (`ProcrustesAlign()`) and resistant-fit superimposition
  (`RftraAlign()`) for shape data.
- **Distances**: the additive ("squaring") fix for neighbour-joining on
  continuous (`SquaredEuclideanDistance()`) and discrete
  (`CorrectedDistanceMS()`) data.
- **Tree inference**: wrappers for maximum-parsimony search in TNT
  (`RunTNTParsimony()`, continuous/discrete/landmark data, using the
  unified search protocol validated throughout the reply) and Bayesian
  inference in RevBayes (`RunRevBayesBM()`, `RunRevBayesMk()`).
- **Scoring**: tree-distance metrics (`ScoreToTruth()`: clustering
  information distance, Robinson-Foulds, SPR), exact-match and
  optimality-based recovery (`ExactlyRecovered()`, `OptimalityRecovery()`),
  and maximum-entropy subsampling of tied most-parsimonious trees
  (`ScoreMPTsByMaxEntropy()`).
- **Cached results**: the 250 generative trees (`generative_trees`) and
  every dataset behind the reply's figures (`scale`, `scale_fig1_dists`,
  `shape_combine_nj`, and others -- see `vignette("scoring-and-figures")`
  or `help(package = "phyloPCA")` for the full list). Bundling these means
  the reply's results can be inspected and re-plotted without re-running
  the underlying HPC-scale simulation grids.

## Getting started

```r
library(phyloPCA)
data(generative_trees)
truth <- generative_trees[[1]]

M <- SimulateContinuous(truth, nChar = 100)
scores <- prcomp(M)$x

treeNJ <- ape::nj(SquaredEuclideanDistance(scores))
ScoreToTruth(treeNJ, truth, metrics = c("CID", "RF"))
```

See the vignettes for a full walkthrough of each stage of the analysis:

1. [Generative trees](https://phylo-pca.github.io/articles/generative-trees.html)
2. [Simulating character data](https://phylo-pca.github.io/articles/simulating-data.html)
3. [PCA, alignment and distances](https://phylo-pca.github.io/articles/pca-and-alignment.html)
4. [Tree inference](https://phylo-pca.github.io/articles/tree-inference.html)
5. [Scoring reconstructions and reproducing the figures](https://phylo-pca.github.io/articles/scoring-and-figures.html)

## Citation

```r
citation("phyloPCA")
```

## Licence

GPL (>= 3). Mongle et al. (2023)'s character matrix is bundled
(`inst/extdata/mongle_2023.nex`) for convenience, with attribution; it
remains the original authors' data.
