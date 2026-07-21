#' Locate the RevBayes executable
#'
#' Looks for a RevBayes installation, in order: the `PhyloPCA.rb` option,
#' the `RB_PATH` environment variable, then `rb`/`rb.exe` on the system
#' `PATH`.
#'
#' @return Path to the RevBayes executable.
#'
#' @export
FindRevBayes <- function() {
  opt <- getOption("PhyloPCA.rb")
  if (!is.null(opt) && file.exists(opt)) {
    # Return: user-configured path
    return(opt)
  }
  env <- Sys.getenv("RB_PATH", NA)
  if (!is.na(env) && file.exists(env)) {
    # Return: environment-configured path
    return(env)
  }
  onPath <- Sys.which(c("rb", "rb.exe"))
  onPath <- onPath[nzchar(onPath)]
  if (length(onPath)) {
    # Return: first match found on PATH
    return(unname(onPath[1]))
  }
  stop("Could not locate RevBayes. Set options(PhyloPCA.rb = '/path/to/rb') ",
       "or the RB_PATH environment variable. RevBayes is available from ",
       "https://revbayes.github.io/.")
}

.RunRevBayes <- function(template, inNex, outPre, rbPath, maxHours, seed = NULL,
                          extraVars = list()) {
  header <- tempfile(fileext = ".Rev")
  vars <- c(list(inNex = inNex, outPre = outPre, maxHours = maxHours), extraVars)
  lines <- vapply(names(vars), function(k) {
    v <- vars[[k]]
    if (is.character(v)) sprintf('%s = "%s"', k, v) else sprintf("%s = %s", k, v)
  }, character(1))
  seedLine <- if (is.null(seed)) character(0) else sprintf("seed(%d)", as.integer(seed))
  writeLines(c(seedLine, lines, sprintf('source("%s")', template)), header)
  log <- paste0(outPre, ".rb.log")
  status <- system2(rbPath, shQuote(header), stdout = log, stderr = log)
  if (status != 0) {
    warning("RevBayes exited with status ", status, "; see ", log)
  }
  mapFiles <- paste0(outPre, c("_map1.tre", "_map2.tre"))
  # Return: per-run MAP tree paths + the log file, whether or not they exist
  list(mapTrees = mapFiles[file.exists(mapFiles)], log = log)
}

#' Bayesian inference of a tree from continuous character data
#'
#' Runs the two-run, Metropolis-coupled Bayesian phylogenetic analysis of
#' continuous character data (a Brownian motion model, `dnPhyloBrownianREML`)
#' used throughout Smith (2026), for both simulated continuous morphological
#' characters and Procrustes/RFTRA-aligned shape coordinates (both are
#' continuous data as far as the model is concerned). Each run stops itself
#' once the Gelman-Rubin PSRF falls below 1.01 and the parameter ESS exceeds
#' 333 (or `maxHours` elapses).
#'
#' @param nexusFile Path to a continuous-character Nexus file (see
#'   [WriteContinuousNexus()]).
#' @param outputPrefix File prefix for RevBayes' outputs (`<prefix>.trees`,
#'   `<prefix>.p.log`, `<prefix>_map1.tre`, `<prefix>_map2.tre`, ...).
#' @param rbPath Path to the RevBayes executable; see [FindRevBayes()].
#' @param maxHours Wall-clock stopping rule, in hours.
#' @param seed Optional integer random seed, for a reproducible run (RevBayes
#'   uses its own internal default seeding if omitted).
#'
#' @return A list with elements `mapTrees` (paths to the per-run maximum
#'   *a posteriori* tree files that were produced) and `log` (the RevBayes
#'   console log file).
#'
#' @export
RunRevBayesBM <- function(nexusFile, outputPrefix, rbPath = FindRevBayes(), maxHours = 2,
                           seed = NULL) {
  template <- system.file("rev", "bm_infer_body.Rev", package = "PhyloPCA")
  # Return: MAP tree file paths + log
  .RunRevBayes(template, nexusFile, outputPrefix, rbPath, maxHours, seed)
}

#' Bayesian inference of a tree from discrete character data
#'
#' Runs the two-run, Metropolis-coupled Bayesian phylogenetic analysis of
#' multistate discrete character data used throughout Smith (2026): a
#' per-state-size Mk model (Lewis 2001; `fnJC(k)` for each observed
#' state-space size `k`, matching [SimulateDiscreteMS()]), conditioned on
#' variable characters only (`coding = "variable"`). Each run stops itself
#' once the Gelman-Rubin PSRF falls below 1.01 and the parameter ESS exceeds
#' 333 (or `maxHours` elapses).
#'
#' @inheritParams RunRevBayesBM
#' @param nexusFile Path to a discrete-character Nexus file (see
#'   [WriteDiscreteNexus()]).
#'
#' @return A list with elements `mapTrees` (paths to the per-run maximum
#'   *a posteriori* tree files that were produced) and `log` (the RevBayes
#'   console log file).
#'
#' @export
RunRevBayesMk <- function(nexusFile, outputPrefix, rbPath = FindRevBayes(), maxHours = 2,
                           seed = NULL) {
  template <- system.file("rev", "mk_infer_body.Rev", package = "PhyloPCA")
  # Return: MAP tree file paths + log
  .RunRevBayes(template, nexusFile, outputPrefix, rbPath, maxHours, seed)
}

#' Read a RevBayes maximum *a posteriori* tree
#'
#' A thin wrapper around [ape::read.nexus()] for the annotated `_map*.tre`
#' files written by RevBayes' `mapTree()` (used as the single representative
#' tree for each Bayesian analysis in Smith 2026, since it is the topology
#' most sampled in the posterior, and thus has the highest posterior
#' probability).
#'
#' @param file Path to a RevBayes MAP tree file.
#'
#' @return A `phylo` object.
#'
#' @importFrom ape read.nexus
#' @export
ReadMAPTree <- function(file) {
  # Return: the MAP tree, annotations stripped
  read.nexus(file)
}
