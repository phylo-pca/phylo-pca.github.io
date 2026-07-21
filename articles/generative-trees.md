# 1. Generative trees

Every simulation in Smith (2026) is generated along one of 250 “true”
trees, so that reconstruction accuracy can be judged by comparison back
to a known answer. This vignette explains where those trees came from,
and how to regenerate (or extend) the sample.

## The bundled trees

Regenerating the tree posterior takes hours (see below), so the actual
250-tree sample used throughout the package is bundled directly:

``` r

library(phyloPCA)
#> Registered S3 method overwritten by 'phangorn':
#>   method   from     
#>   [.phyDat TreeTools
data(generative_trees)
generative_trees
#> 250 phylogenetic trees
```

``` r

plot(generative_trees[[1]], cex = 0.6, no.margin = TRUE)
```

![](generative-trees_files/figure-html/unnamed-chunk-3-1.png)

Trees are indexed stably: `generative_trees[[k]]` always refers to the
same tree, however many trees are sampled in total (see “How the sample
is built”, below) – so a tree index used elsewhere in this package (e.g.
`scale$tree`, `shape_combine_nj$t`) unambiguously identifies one of
these trees.

## Where they come from

The trees are a posterior sample from a Bayesian re-inference of Mongle
et al. (2023)’s discrete character matrix (15 hominins, 6 anthropoid
outgroups, 107 characters), under the Mk model (Lewis 2001), matching
Raskin et al. (2026)’s own inference model exactly except for the
ascertainment-bias correction: `coding = "variable"`, the correct choice
for a matrix with no invariant characters. A small excerpt of the matrix
is bundled for illustration:

``` r

mongleFile <- system.file("extdata", "mongle_2023.nex", package = "phyloPCA")
MongleStateProportions(mongleFile)
#> 
#>  2  3  4  5  6 
#> 34 44 20  8  1
```

This tabulates, for each observed number of states `k`, how many of the
107 characters attain it – the distribution
[`SimulateDiscreteMS()`](https://phylo-pca.github.io/reference/SimulateDiscreteMS.md)
draws from when simulating new discrete characters (see
[`vignette("simulating-data")`](https://phylo-pca.github.io/articles/simulating-data.md)).

### Re-running the re-inference (not required)

[`RunRevBayesReinference()`](https://phylo-pca.github.io/reference/RunRevBayesReinference.md)
runs this model via RevBayes. It takes hours even for 21 taxa, so this
chunk does not run automatically:

``` r

result <- RunRevBayesReinference(
  nexusFile = mongleFile,
  outDir = tempfile("reinfer")
)
result$runFiles # the two independent .trees traces
```

### From trace to tree sample

Given two converged tree traces,
[`SampleStableTrees()`](https://phylo-pca.github.io/reference/SampleStableTrees.md)
draws a fixed, append-only sample: trees from both runs are interleaved
(so both chains contribute from the very first draw) and then strided,
so that `tree_k` is a fixed function of `k` – raising the sample size
later only appends new trees, never redraws existing ones. Convergence
is checked first via
[`TreeESS()`](https://phylo-pca.github.io/reference/TreeESS.md), the
topological effective-sample-size diagnostic (the lower of the
Frechet-correlation and median-pseudo ESS statistics; Fabreti & Hohna
2022; Magee et al. 2024) used throughout Smith (2026):

``` r

runFiles <- result$runFiles
trees250 <- SampleStableTrees(runFiles, nTrees = 250, stride = 11)
```

`stride = 11` is large enough that consecutive draws from the
interleaved pool are effectively independent;
[`SampleStableTrees()`](https://phylo-pca.github.io/reference/SampleStableTrees.md)
errors informatively if a larger sample is requested than the trace
supports, rather than silently reusing correlated draws.
