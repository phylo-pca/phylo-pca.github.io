# Publication-scale grid: continuous and discrete data, NJ and Bayesian inference

The core simulation grid behind Figures 1 and 2 of Smith (2026): 250
generative trees (see
[`SampleStableTrees()`](https://phylo-pca.github.io/reference/SampleStableTrees.md))
x `nC` in `{25, 50, 100, 250}` characters x data type (`continuous`,
simulated by
[`SimulateContinuous()`](https://phylo-pca.github.io/reference/SimulateContinuous.md);
`discrete`, by
[`SimulateDiscreteMS()`](https://phylo-pca.github.io/reference/SimulateDiscreteMS.md))
x method (neighbour-joining, several distances; Bayesian inference,
matched model). Continuous Bayesian cells are superseded by the
converged 2-chain re-run in `scale_cont_rerun`; discrete Bayesian cells
here are the final, 2-chain, converged values.

## Usage

``` r
data(scale)
```

## Format

A data frame with one row per (tree, nC, data type) cell:

- tree, treeIdx:

  Index of the generative tree, `1:250`.

- nC:

  Number of characters simulated.

- data:

  `"continuous"` or `"discrete"`.

- NJ_pipeline:

  CID of the neighbour-joining tree built on Euclidean distance between
  the first two principal components (Raskin et al. 2026's pipeline;
  discrete arm uses Hamming distance).

- NJ_allpc:

  CID of the neighbour-joining tree built on Euclidean distance using
  all principal components.

- NJ_squared:

  CID of the neighbour-joining tree built on the additive
  (squared-Euclidean) distance (continuous), or the JC-corrected
  additive distance
  ([`CorrectedDistanceMS()`](https://phylo-pca.github.io/reference/CorrectedDistanceMS.md);
  discrete).

- NJ_best:

  The better of `NJ_allpc`/`NJ_squared` for that cell.

- Bayes_MAP:

  CID of the Bayesian maximum *a posteriori* tree.

- Bayes_pm:

  Posterior-mean CID (averaged across sampled trees).

- P_true:

  Posterior probability of the exact generative topology.

- rhat:

  Gelman-Rubin PSRF convergence diagnostic for that cell.

## Source

`hominin-pca-scale/R/21_merge_scale.R`, from
`hominin-pca-scale/R/10_scale_run.R` (simulation + neighbour-joining)
and `inst/rev/bm_infer_body.Rev` / `inst/rev/mk_infer_body.Rev`
(Bayesian inference).

## See also

`scale_cont_rerun` for the converged continuous Bayesian re-run this
dataset's continuous `Bayes_*` columns should be read from instead.
