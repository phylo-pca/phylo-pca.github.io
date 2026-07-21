#' Locate the TNT executable
#'
#' Looks for a TNT installation, in order: the `phyloPCA.tnt` option, the
#' `TNT_PATH` environment variable, then a handful of common install
#' locations.
#'
#' @return Path to the TNT executable.
#'
#' @export
FindTNT <- function() {
  candidates <- c(
    getOption("phyloPCA.tnt"),
    Sys.getenv("TNT_PATH", NA),
    "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe",
    "/usr/local/bin/tnt", "/opt/tnt/tnt"
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  found <- candidates[file.exists(candidates)]
  if (!length(found)) {
    stop("Could not locate TNT. Set options(phyloPCA.tnt = '/path/to/tnt') ",
         "or the TNT_PATH environment variable. TNT is available from ",
         "http://www.lillo.org.ar/phylogeny/tnt/ (free for the Willi Hennig ",
         "Society membership).")
  }
  # Return: first existing candidate path
  found[1]
}

# TNT's parenthetical tree notation differs from Newick only in whitespace
# conventions; this converts one line of `tsave` output to parseable Newick.
.TntToNewick <- function(x) {
  x <- gsub("\\*", "", x)
  x <- gsub(";\\s*$", "", x)
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  x <- gsub(" \\)", ")", x)
  x <- gsub(")(", "),(", x, fixed = TRUE) # TNT omits the separator here
  x <- gsub(" ", ",", x)
  # Return: a parseable Newick string
  paste0(x, ";")
}

# The unified search protocol validated across the continuous, discrete and
# shape parsimony arms (see dev notes `parsimony-comparator-arms`): 50
# RAS+TBR seed trees, then a New Technology pass (sectorial search + ratchet
# + drift + tree fusing, which never beats plain RAS+TBR here but confirms
# the optimum), then TBR-swap every optimal tree to completion and save ALL
# most-parsimonious trees (`save 0.`) into a single large tree buffer (big
# enough that `bbreak` never overflows into an interactive prompt).
.TntSearchLines <- function(type) {
  lmark <- if (type == "landmark") "lmark xthreads;" else NULL
  slack <- if (type == "landmark") 90 else 50 # shape grids tie more heavily
  c(lmark,
    "hold 100000;", "mult=replic 50 tbr;",
    sprintf("sect: slack %d;", slack),
    "xmult=level 10 hits 20 ratchet 15 drift 15 fuse 5;",
    "bbreak;", "tsave *out.tre;", "save 0.;", "tsave/;", "quit;")
}

.TntWriteContinuous <- function(M, path, rseed) {
  # Shift each character to a non-negative minimum (a translation, hence
  # topology-neutral for additive parsimony); down-scale only if any
  # character's range nears TNT's internal ~65-state cap.
  shifted <- sweep(M, 2, apply(M, 2, min), "-")
  range <- max(apply(shifted, 2, function(x) diff(range(x))))
  if (range > 60) shifted <- shifted * (60 / range)
  ntax <- nrow(shifted)
  nChar <- ncol(shifted)
  body <- apply(cbind(sprintf("t%d", seq_len(ntax) - 1L),
                      matrix(sprintf("%.6f", shifted), ntax)), 1, paste, collapse = " ")
  writeLines(c(
    "mxram 500;", "nstates cont;", sprintf("rseed %d;", rseed),
    "xread", "'phyloPCA continuous parsimony'", paste(nChar, ntax), body, ";",
    .TntSearchLines("continuous")
  ), path)
}

.TntWriteDiscrete <- function(X, path, rseed) {
  ntax <- nrow(X)
  nChar <- ncol(X)
  sequences <- apply(X, 1, paste, collapse = "")
  body <- paste(sprintf("t%d", seq_len(ntax) - 1L), sequences)
  writeLines(c(
    "mxram 500;", sprintf("rseed %d;", rseed),
    "xread", "'phyloPCA discrete parsimony'", paste(nChar, ntax), body, ";",
    .TntSearchLines("discrete")
  ), path)
}

.TntWriteLandmark <- function(configs, dim, path, rseed) {
  ntax <- length(configs)
  fmt <- function(config) {
    paste(apply(config, 1, function(row) paste(sprintf("%.5f", row), collapse = ",")),
          collapse = " ")
  }
  body <- paste(sprintf("t%d", seq_len(ntax) - 1L), vapply(configs, fmt, character(1)))
  writeLines(c(
    "mxram 500;", "nstates cont;", sprintf("rseed %d;", rseed),
    "xread", "'phyloPCA shape parsimony'", paste(1, ntax),
    sprintf("&[landmark %dd]", dim), body, ";",
    .TntSearchLines("landmark")
  ), path)
}

#' Maximum-parsimony tree search in TNT
#'
#' Runs a maximum-parsimony search in TNT (Goloboff et al. 2008; Goloboff &
#' Catalano 2016) using the unified search protocol validated in Smith
#' (2026) for continuous (Wagner), discrete (Fitch) and geometric-
#' morphometric landmark (Goloboff-Catalano linear displacement, following
#' Palci & Lee 2019) data, and returns every most-parsimonious tree found.
#'
#' @param data For `type = "continuous"`: a taxa x characters numeric
#'   matrix. For `type = "discrete"`: a taxa x characters matrix of
#'   single-character state labels. For `type = "landmark"`: a named list of
#'   taxa's landmarks x `dim` numeric matrices (e.g. from [ProcrustesAlign()]
#'   or [RftraAlign()]).
#' @param type Data type to search under; see `data`.
#' @param dim Number of spatial dimensions (required, and only used, when
#'   `type = "landmark"`).
#' @param tntPath Path to the TNT executable; see [FindTNT()].
#' @param workDir Working directory for the (purely alphabetic, per TNT's
#'   filename-parsing quirk) job script and tree file; created if absent.
#' @param timeout Maximum seconds to allow the search to run before
#'   aborting (TNT occasionally hangs on a malformed/overflowing buffer).
#' @param retrySeed A second `rseed` to retry with, once, if the first
#'   attempt fails or times out.
#'
#' @return A list with elements `trees` (a `multiPhylo` of every
#'   most-parsimonious tree found, or `NULL` on failure) and `nTree` (the
#'   number of most-parsimonious trees).
#'
#' @importFrom ape read.tree unroot
#' @export
RunTNTParsimony <- function(data, type = c("continuous", "discrete", "landmark"),
                             dim = NULL, tntPath = FindTNT(), workDir = tempfile("tnt"),
                             timeout = 300, retrySeed = 7L) {
  type <- match.arg(type)
  dir.create(workDir, showWarnings = FALSE, recursive = TRUE)
  runFile <- "job.run"
  treeFile <- "out.tre"

  taxa <- if (type == "landmark") names(data) else rownames(data)
  attempt <- function(rseed) {
    path <- file.path(workDir, runFile)
    switch(type,
      continuous = .TntWriteContinuous(data, path, rseed),
      discrete = .TntWriteDiscrete(data, path, rseed),
      landmark = {
        stopifnot("`dim` is required when type = 'landmark'" = !is.null(dim))
        .TntWriteLandmark(data, dim, path, rseed)
      }
    )
    outPath <- file.path(workDir, treeFile)
    if (file.exists(outPath)) file.remove(outPath)
    oldwd <- getwd()
    on.exit(setwd(oldwd))
    setwd(workDir)
    tryCatch(
      system2(tntPath, sprintf('"%s;"', runFile), stdout = FALSE, stderr = FALSE,
              timeout = timeout),
      error = function(e) 124L
    )
    if (!file.exists(outPath) || file.size(outPath) == 0) {
      # Return: no tree output produced
      return(NULL)
    }
    raw <- iconv(readLines(outPath, warn = FALSE), from = "", to = "UTF-8", sub = "")
    treeLines <- grep("^\\s*\\(", raw, value = TRUE)
    if (!length(treeLines)) {
      # Return: output file existed but held no parseable trees
      return(NULL)
    }
    trees <- lapply(treeLines, function(s) {
      t <- read.tree(text = .TntToNewick(s))
      t$tip.label <- taxa[as.integer(t$tip.label) + 1L]
      unroot(t)
    })
    class(trees) <- "multiPhylo"
    # Return: every most-parsimonious tree found
    trees
  }

  trees <- attempt(1L)
  if (is.null(trees) && !is.null(retrySeed)) trees <- attempt(retrySeed)
  # Return: search result (trees = NULL signals search failure/timeout)
  list(trees = trees, nTree = if (is.null(trees)) NA_integer_ else length(trees))
}

#' Score a most-parsimonious tree set by maximum-entropy subsampling
#'
#' Where a parsimony search finds multiple, equally most-parsimonious trees,
#' a strict consensus is over-conservative (it penalises a cell for
#' ambiguity among equally good trees, not for getting the wrong answer). As
#' in Smith (2026), a density-blind subsample of up to `k` trees is instead
#' selected to maximize the entropy (spread) within the chosen subset (via
#' `MaxMin::MaxEntropy()`, using the clustering information distance as the
#' measure of redundancy between trees), and their mean distance to the true
#' tree is reported. When at most `k` trees were found, every tree is used.
#'
#' @param trees A `multiPhylo` of most-parsimonious trees (e.g. from
#'   [RunTNTParsimony()]).
#' @param truth The true (generative) tree.
#' @param k Maximum subsample size.
#' @param pool If more than `pool` trees are supplied, they are first
#'   evenly thinned to `pool` trees (most-parsimonious trees are
#'   exchangeable, so this bounds the O(n^2) distance matrix on
#'   pathologically tie-heavy cells without materially changing the mean).
#'
#' @return A list with elements `CID` (mean clustering information distance
#'   of the subsample to `truth`), `RF` (mean Robinson-Foulds distance),
#'   `nTree` (the number of most-parsimonious trees supplied) and `nSel`
#'   (the subsample size actually used).
#'
#' @importFrom ape unroot
#' @importFrom TreeDist ClusteringInfoDistance RobinsonFoulds
#' @export
ScoreMPTsByMaxEntropy <- function(trees, truth, k = 10L, pool = 200L) {
  if (!requireNamespace("MaxMin", quietly = TRUE)) {
    stop("Package 'MaxMin' is required: remotes::install_github('ms609/MaxMin')")
  }
  nFull <- length(trees)
  if (nFull > pool) {
    trees <- trees[round(seq(1, nFull, length.out = pool))]
  }
  n <- length(trees)
  truth <- unroot(truth)
  sel <- if (n <= k) {
    seq_len(n)
  } else {
    D <- matrix(0, n, n)
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        D[i, j] <- D[j, i] <- ClusteringInfoDistance(unroot(trees[[i]]), unroot(trees[[j]]),
                                                      normalize = TRUE)
      }
    }
    as.integer(MaxMin::MaxEntropy(k, stats::as.dist(D)))
  }
  cidTo <- vapply(trees[sel], function(t) {
    ClusteringInfoDistance(unroot(t), truth, normalize = TRUE)
  }, numeric(1))
  rfTo <- vapply(trees[sel], function(t) {
    as.numeric(RobinsonFoulds(unroot(t), truth, normalize = FALSE))
  }, numeric(1))
  # Return: max-entropy-subsample mean distances + subsample bookkeeping
  list(CID = mean(cidTo), RF = mean(rfTo), nTree = nFull, nSel = length(sel))
}

#' Optimality-based tree recovery
#'
#' Whether the generative tree is (tied for) most parsimonious, following
#' the recovery definition used for the parsimony arms in Smith (2026): if
#' a heuristic search misses the true topology, that is the search
#' heuristic's fault, not a failure of parsimony as a criterion, so recovery
#' is judged by comparing the true tree's own parsimony score to the best
#' score found, rather than by topological identity with a most-parsimonious
#' tree.
#'
#' @param trueLength Parsimony score of the generative tree.
#' @param bestLength Best (minimum) parsimony score found by the search.
#' @param tolerance Numeric tolerance for the comparison (guards against
#'   floating-point score differences for continuous/landmark data).
#'
#' @return A logical: whether the generative tree is tied for most
#'   parsimonious.
#'
#' @export
OptimalityRecovery <- function(trueLength, bestLength, tolerance = 1e-4) {
  # Return: TRUE iff the true tree's score does not exceed the best found
  trueLength <= bestLength + tolerance
}
