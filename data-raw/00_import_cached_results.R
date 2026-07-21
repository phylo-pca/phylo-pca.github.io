#!/usr/bin/env Rscript
# Import the cached results underlying the three figures of Smith (2026), a
# reply to Raskin et al. (2026). These are genuinely cached: the full grids
# were computed on Hamilton (Durham University's HPC service) and locally
# over several days, using the functions in R/ (simulate.R, align.R,
# parsimony.R, bayesian.R, score.R) orchestrated by scripts that are not
# part of this package (see vignette("05-scoring-and-figures") for how each
# object maps onto a figure panel, and the source chain that produced it).
#
# Source, worktree-first: on the author's machine the live analysis
# worktree (hominin-pca-scale) is present and wins, so re-running the
# analysis is never shadowed by a stale copy; everyone else falls back to
# the snapshot vendored under data-raw/cache, which this script also
# refreshes whenever the live worktree *is* available.
SCALE <- "C:/Users/pjjg18/GitHub/hominin-pca-scale/out"
CACHE <- file.path("data-raw", "cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
dir.create("data", showWarnings = FALSE)

files <- c(
  "scale", "scale_cont_rerun", "scale_njdense", "scale_fig1_dists",
  "shape_combine_nj", "shape_sigma_bayes", "shape_both_bayes",
  "scale_disc_parsimony", "score_native", "score_rftra",
  "shape_sigma_nj_align", "shape_sigma_bayes_align"
)

for (nm in files) {
  live <- file.path(SCALE, paste0(nm, ".rds"))
  cached <- file.path(CACHE, paste0(nm, ".rds"))
  src <- if (file.exists(live)) live else cached
  if (file.exists(live)) file.copy(live, cached, overwrite = TRUE)
  stopifnot("cached snapshot is missing and the live worktree is unavailable" =
              file.exists(src))
  assign(nm, readRDS(src))
  save(list = nm, file = file.path("data", paste0(nm, ".rda")), compress = "xz")
  cat(sprintf("data/%s.rda (%d rows)\n", nm, NROW(get(nm))))
}
