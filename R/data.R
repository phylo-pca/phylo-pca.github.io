#' The 250 generative trees used throughout Smith (2026)
#'
#' Every simulation and reconstruction in Smith (2026) is scored against one
#' of these 250 trees: a stable-stride sample (see [SampleStableTrees()])
#' from two independent, converged RevBayes posteriors (see
#' [RunRevBayesReinference()]) obtained by re-inferring Mongle et al.
#' (2023)'s discrete character matrix for 15 hominins and 6 anthropoid
#' outgroups under the Mk model (Lewis 2001) with the `coding = "variable"`
#' ascertainment-bias correction. `tree_k` is a fixed function of `k`
#' (raising the tree count only appends new trees), so trees referenced by
#' index elsewhere in this package (e.g. `scale$tree`, `shape_combine_nj$t`)
#' always refer to the same tree.
#'
#' @format A `multiPhylo` of 250 trees, named `tree_001`..`tree_250`.
#' @source `hominin-pca-scale/R/20_sample_trees.R`, from
#'   `hominin-pca-scale/out/reinfer/reinfer_run_{1,2}.trees` (see
#'   `inst/rev/reinfer_variable.Rev`).
"generative_trees"

#' Publication-scale grid: continuous and discrete data, NJ and Bayesian inference
#'
#' The core simulation grid behind Figures 1 and 2 of Smith (2026): 250
#' generative trees (see [SampleStableTrees()]) x `nC` in `{25, 50, 100,
#' 250}` characters x data type (`continuous`, simulated by
#' [SimulateContinuous()]; `discrete`, by [SimulateDiscreteMS()]) x method
#' (neighbour-joining, several distances; Bayesian inference, matched
#' model). Continuous Bayesian cells are superseded by the converged 2-chain
#' re-run in `scale_cont_rerun`; discrete Bayesian cells here are the final,
#' 2-chain, converged values.
#'
#' @format A data frame with one row per (tree, nC, data type) cell:
#' \describe{
#'   \item{tree, treeIdx}{Index of the generative tree, `1:250`.}
#'   \item{nC}{Number of characters simulated.}
#'   \item{data}{`"continuous"` or `"discrete"`.}
#'   \item{NJ_pipeline}{CID of the neighbour-joining tree built on Euclidean
#'     distance between the first two principal components (Raskin et al.
#'     2026's pipeline; discrete arm uses Hamming distance).}
#'   \item{NJ_allpc}{CID of the neighbour-joining tree built on Euclidean
#'     distance using all principal components.}
#'   \item{NJ_squared}{CID of the neighbour-joining tree built on the
#'     additive (squared-Euclidean) distance (continuous), or the
#'     JC-corrected additive distance ([CorrectedDistanceMS()]; discrete).}
#'   \item{NJ_best}{The better of `NJ_allpc`/`NJ_squared` for that cell.}
#'   \item{Bayes_MAP}{CID of the Bayesian maximum *a posteriori* tree.}
#'   \item{Bayes_pm}{Posterior-mean CID (averaged across sampled trees).}
#'   \item{P_true}{Posterior probability of the exact generative topology.}
#'   \item{rhat}{Gelman-Rubin PSRF convergence diagnostic for that cell.}
#' }
#' @source `hominin-pca-scale/R/21_merge_scale.R`, from
#'   `hominin-pca-scale/R/10_scale_run.R` (simulation + neighbour-joining)
#'   and `inst/rev/bm_infer_body.Rev` / `inst/rev/mk_infer_body.Rev`
#'   (Bayesian inference).
#' @seealso `scale_cont_rerun` for the converged continuous Bayesian re-run
#'   this dataset's continuous `Bayes_*` columns should be read from instead.
"scale"

#' Converged continuous Bayesian re-run
#'
#' A 2-chain, longer re-run of the continuous Bayesian inference cells in
#' `scale`, confirming convergence (all cells ASDSF < 0.01) at `nC` in
#' `{25, 50, 100, 250}`. This is the dataset `scale`'s continuous arm and
#' Figures 1 and 2 actually read the Bayesian numbers from.
#'
#' @format A data frame with columns `tree`, `nC`, `Bayes_MAP` (CID of the
#'   MAP tree), `MAP_CID` (an independent re-scoring, kept for a
#'   cross-check), `Bayes_pm` (posterior-mean CID), `P_true` (posterior
#'   probability of the exact tree) and `rhat` (PSRF).
#' @source `hominin-pca-scale/R/14_rerun_cont.R` + `R/16_add_cont25.R`,
#'   using `inst/rev/bm_infer_body.Rev`.
"scale_cont_rerun"

#' Dense neighbour-joining ladder
#'
#' Neighbour-joining-only accuracy (no Bayesian inference) across a dense
#' range of character counts, `nC` in `{10, 25, 50, 100, 250, 500}`, for
#' both continuous and discrete data, several replicate simulations per
#' (tree, nC) cell.
#'
#' @format A data frame with columns `tree`, `rep` (replicate index), `nC`,
#'   `arm` (`"cont_raw"`, `"cont_sq"`, `"cont_pc12"` or `"disc_cor"`; see
#'   `scale`'s `NJ_*` columns for the equivalent distances) and `CID`.
#' @source `hominin-pca-scale/R/09_scale_njdense.R`.
"scale_njdense"

#' Figure 1 method-comparison distances (continuous data)
#'
#' The five method arms compared in Figure 1 of Smith (2026), all on the
#' identical continuous character matrices used in `scale`: the
#' PCA-\>neighbour-joining pipeline, all-dimension and additive-distance
#' neighbour-joining, maximum parsimony (TNT, continuous Wagner; see
#' [RunTNTParsimony()]) and Bayesian inference, each scored by clustering
#' information distance, Robinson-Foulds distance and (where both trees are
#' binary) SPR distance.
#'
#' @format A data frame with columns `tree`, `nC`, `arm` (`"pipeline"`,
#'   `"allpc"`, `"additive"`, `"mp"` or `"bayes"`), `CID`, `RF`, `SPR` and
#'   `SPR_exact` (whether the SPR distance could be computed exactly; see
#'   [ScoreToTruth()]).
#' @source `hominin-pca-scale/R/15_fig1_rfspr.R` (pipeline/allpc/
#'   additive/bayes arms), `R/34_par_hamilton_gen.R` + `R/35_par_hamilton_score.R`
#'   (the `"mp"` arm, run on Hamilton under the unified TNT search protocol;
#'   see [RunTNTParsimony()]).
#' @note The `"mp"` arm's `SPR` was recomputed on Hamilton via
#'   `TBRDist::USPRDist()` under a per-cell subprocess timeout (999/1000
#'   cells; `SPR_exact = TRUE`). The one cell that did not resolve within a
#'   6-hour budget (a genuinely hard NP-hard instance, not a bug: two
#'   markedly dissimilar 21-tip trees, RF = 20) is filled from
#'   `TreeDist::SPRDist()`'s polynomial-time approximation instead
#'   (`SPR_exact = FALSE`), matching the convention already used for the
#'   other arms' occasional approximate fallbacks.
"scale_fig1_dists"

#' Combining 2D and 3D shape data (neighbour-joining)
#'
#' Neighbour-joining accuracy from 2D shape data, 3D shape data, and the two
#' combined, across the landmark-integration parameter `sigma`, on a shared
#' set of trees so the three arms are paired. `cid2`/`cid3`/`cidB` are
#' scored on raw Euclidean distance; `cid2sq`/`cid3sq`/`cidBsq` on the
#' additive (squared-Euclidean) distance, with the combined arm computed as
#' `nj(d2^2 + d3^2)` (the correct way two independent additive distances
#' fuse).
#'
#' @format A data frame with columns `t` (tree index), `sigma`,
#'   `cid2`/`cid3`/`cidB` (raw-distance CID for 2D/3D/combined),
#'   `cid2sq`/`cid3sq`/`cidBsq` (additive-distance CID) and `scale_ratio`
#'   (relative scale of the 2D vs 3D coordinate sets before combining).
#' @source `hominin-pca-scale/R/22_shape_combine_nj.R`, from LDDMM shape
#'   simulations aligned by [ProcrustesAlign()].
"shape_combine_nj"

#' Shape data Bayesian inference by dimension and sigma
#'
#' Bayesian inference accuracy (Brownian motion model on Procrustes-aligned
#' coordinates; see [RunRevBayesBM()]) for 2D and 3D shape data separately,
#' across the landmark-integration parameter `sigma`.
#'
#' @format A data frame with columns `tnum` (tree index), `dim` (2 or 3),
#'   `sigma`, `Bayes_MAP` (CID of the MAP tree), `MAP_CID` (independent
#'   re-scoring cross-check), `Bayes_pm` (posterior-mean CID), `P_true`
#'   (posterior probability of the exact tree) and `rhat` (PSRF).
#' @source `hominin-pca-scale/R/20b_sigma_bayes_score_maps.R`, using
#'   `inst/rev/bm_infer_body.Rev` on Hamilton.
"shape_sigma_bayes"

#' Combined 2D+3D shape data Bayesian inference
#'
#' Bayesian inference accuracy (Brownian motion on the concatenated,
#' Procrustes-aligned 2D and 3D coordinate matrix) across `sigma`, paired
#' with the individual-dimension arm in `shape_sigma_bayes` restricted to
#' the same trees (`unique(shape_both_bayes$t)`).
#'
#' @format A data frame with columns `t` (tree index), `sigma`, `Bayes_MAP`
#'   (CID of the MAP tree), `MAP_CID` (independent re-scoring cross-check),
#'   `P_true` and `ASDSF` (average standard deviation of split frequencies,
#'   a convergence diagnostic).
#' @source `hominin-pca-scale/R/24_both_bayes_score.R`, using
#'   `inst/rev/bm_infer_body.Rev` on Hamilton (`rev/both_array.sh`).
"shape_both_bayes"

#' Discrete data maximum parsimony
#'
#' Maximum-parsimony accuracy (TNT, standard Fitch optimization on
#' unordered multistate characters; see [RunTNTParsimony()]) on the
#' identical discrete character matrices used in `scale`, scored by a
#' maximum-entropy subsample of most-parsimonious trees (see
#' [ScoreMPTsByMaxEntropy()]) rather than a strict consensus.
#'
#' @format A data frame with columns `tree`, `nC`, `CID`, `RF`, `SPR`
#'   (`NA`; not computed for parsimony arms, see [ScoreToTruth()]) and
#'   `nmpt` (number of most-parsimonious trees found).
#' @source `hominin-pca-scale/R/34_par_hamilton_gen.R` +
#'   `R/35_par_hamilton_score.R`, run on Hamilton under the unified TNT
#'   search protocol.
"scale_disc_parsimony"

#' Shape data maximum parsimony: native/Procrustes alignment
#'
#' Maximum-parsimony accuracy on shape (landmark) data (TNT's landmark data
#' type; the linear-displacement method of Goloboff & Catalano, following
#' Palci & Lee 2019; see [RunTNTParsimony()]), comparing the "native"
#' (TNT's own internal linear superimposition, no prior alignment) and
#' Procrustes-aligned ([ProcrustesAlign()]) treatments, across dimension and
#' `sigma`.
#'
#' @format A data frame with columns `idx` (job index), `t0` (tree index,
#'   0-based), `dim`, `sigma`, `align` (`"proc"` or `"raw"`), `CID`, `RF`,
#'   `nmpt` (most-parsimonious trees found) and `ksel` (max-entropy
#'   subsample size used).
#' @source `hominin-pca-scale/R/24_shape_par_gen.R` +
#'   `R/26_shape_par_score.R`, run on Hamilton.
#' @note A re-run under the fuller unified search protocol (adding
#'   `lmark xthreads;` and a larger sectorial-search buffer; see
#'   [RunTNTParsimony()]) was in progress on Hamilton at the time this
#'   snapshot was taken. Prior validation showed the search was already at
#'   the global optimum, so the re-run is expected to confirm rather than
#'   change these values, but this file will be refreshed if it does.
"score_native"

#' Shape data maximum parsimony: resistant-fit (RFTRA) alignment
#'
#' As `score_native`, but using resistant-fit superimposition ([RftraAlign()])
#' rather than native or Procrustes alignment.
#'
#' @format As `score_native`, with `align` fixed to `"rftra"`.
#' @source `hominin-pca-scale/R/25_shape_rftra_gen.R` +
#'   `R/26_shape_par_score.R`, run on Hamilton.
#' @note See the note on `score_native`: a confirmatory re-run under the
#'   fuller unified search protocol was in progress when this snapshot was
#'   taken.
"score_rftra"

#' Shape data neighbour-joining: alignment comparison
#'
#' Neighbour-joining accuracy on shape data under three alignments -
#' Procrustes ([ProcrustesAlign()]), RFTRA ([RftraAlign()]) and "native"
#' (no superimposition) - across dimension and `sigma`, the neighbour-
#' joining half of the alignment-comparison Figure (`fig_alignment.png`).
#'
#' @format A data frame with columns `t` (tree index), `dim`, `sigma`,
#'   `align` (`"proc"`, `"rftra"` or `"raw"`) and `CID`.
#' @source `hominin-pca-scale/R/28_shape_nj_align.R`.
"shape_sigma_nj_align"

#' Shape data Bayesian inference: alignment comparison
#'
#' As `shape_sigma_nj_align`, but for Bayesian inference (Brownian motion;
#' see [RunRevBayesBM()]) rather than neighbour-joining; the `"proc"` rows
#' are shared with `shape_sigma_bayes`.
#'
#' @format A data frame with columns `t` (tree index), `dim`, `sigma`,
#'   `align` (`"proc"`, `"rftra"` or `"raw"`), `Bayes_MAP` (CID of the MAP
#'   tree) and `MAP_CID` (independent re-scoring cross-check).
#' @source `hominin-pca-scale/R/29_shape_bayes_align_prep.R` +
#'   `R/30_shape_bayes_align_score.R`, using `inst/rev/bm_infer_body.Rev`
#'   on Hamilton (`rev/shape_align_array.sh`). Raw/RFTRA coordinates are
#'   rescaled to match Procrustes' total centred sum-of-squares before
#'   inference, so the Brownian tree-length prior is not railroaded by the
#'   native LDDMM coordinate scale (topology-neutral).
"shape_sigma_bayes_align"
