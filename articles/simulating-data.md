# 2. Simulating character data

``` r

library(PhyloPCA)
#> Registered S3 method overwritten by 'phangorn':
#>   method   from     
#>   [.phyDat TreeTools
data(generative_trees)
truth <- generative_trees[[1]]
```

Smith (2026) simulates three kinds of character data along each
generative tree: continuous morphological characters, discrete
morphological characters, and geometric-morphometric (“shape”) landmark
data. All three are compared on an equal footing at the end of the
analysis (see
[`vignette("scoring-and-figures")`](https://phylo-pca.github.io/articles/scoring-and-figures.md)).

## Continuous characters

Continuous characters are simulated under Brownian motion (root state 0,
rate `sigma^2 = 1`, following Raskin et al. 2026), independently for
each character:

``` r

set.seed(1)
continuous <- SimulateContinuous(truth, nChar = 50)
dim(continuous)
#> [1] 21 50
continuous[1:4, 1:4]
#>                                   [,1]        [,2]       [,3]       [,4]
#> Australopithecus_afarensis  0.37596361  0.41737953 -0.3544226 -0.1598788
#> Paranthropus_boisei        -0.28315036 -0.45622196 -0.4200530  0.3664421
#> Paranthropus_robustus      -0.05806199 -0.03238347  0.2002942  0.8465719
#> Paranthropus_aethiopicus   -0.79145422 -0.13586959 -0.9106652  0.4753420
```

## Discrete characters

Discrete characters are simulated under the Mk model (Lewis 2001), with
no rate variation among characters. Each character’s number of possible
states is drawn from the empirical distribution observed in Mongle et
al. (2023)’s matrix (see
[`vignette("generative-trees")`](https://phylo-pca.github.io/articles/generative-trees.md)),
and characters that do not vary across tips are rejected and re-drawn
(`coding = "variable"`), matching the ascertainment correction applied
at inference time:

``` r

mongleFile <- system.file("extdata", "mongle_2023.nex", package = "PhyloPCA")
stateProportions <- MongleStateProportions(mongleFile)
set.seed(1)
discrete <- SimulateDiscreteMS(truth, nChar = 50, stateProportions = stateProportions)
dim(discrete)
#> [1] 21 50
discrete[1:4, 1:8]
#>                            [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> Australopithecus_afarensis "1"  "1"  "0"  "0"  "0"  "1"  "0"  "0" 
#> Paranthropus_boisei        "1"  "0"  "1"  "1"  "0"  "0"  "0"  "1" 
#> Paranthropus_robustus      "1"  "0"  "1"  "1"  "0"  "0"  "0"  "1" 
#> Paranthropus_aethiopicus   "1"  "0"  "1"  "1"  "0"  "0"  "0"  "1"
```

Every column is variable across tips:

``` r

all(apply(discrete, 2, function(col) length(unique(col)) > 1))
#> [1] TRUE
```

## Shape (geometric-morphometric) data

Shape data are simulated with the bespoke LDDMM (large deformation
diffeomorphic metric mapping) shape simulator of Raskin et al. (2026): a
circular (2D) or spherical (3D) template of evenly spaced landmarks
evolves along each branch under an Euler-Maruyama discretisation of a
Matern-like spatial kernel, rejecting any step that would cause
landmarks to cross.

**This package does not vendor that simulator.** Raskin et al. (2026)’s
repository (<https://github.com/Levi-Raskin/PCAPhylogenetics>) carries
no license, so redistributing their C++ source here would not be a
reproducibility aid worth the ambiguity it creates. To regenerate shape
data from scratch: clone their repository, build `ShapeSimulationCode`
(a recompiled, byte-for-byte-faithful harness is described in the
`phaseC/` directory of the analysis repository at
<https://github.com/phylo-pca/phylo-pca.github.io>), and run it on the
trees in `generative_trees`.

Because most users will never need to run the external simulator, a
small sample of its output is bundled so that everything downstream of
simulation – alignment, tree inference, scoring – can be exercised and
verified end-to-end without it:

``` r

sampleFile <- system.file("extdata", "sample_shape_t00_d2_s0.4.tsv", package = "PhyloPCA")
configs <- ReadShapeConfigs(sampleFile, dim = 2)
length(configs) # one landmark configuration per taxon
#> [1] 21
dim(configs[[1]]) # landmarks x dimensions
#> [1] 25  2
```

This particular sample is 2D, 25-landmark data at a shape kernel
`sigma = 0.4`, for the tree `generative_trees[[1]]` (`t00` in the
0-based indexing used by the simulation grid). See
[`vignette("pca-and-alignment")`](https://phylo-pca.github.io/articles/pca-and-alignment.md)
for what to do with it next.
