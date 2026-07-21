# Maximum-parsimony tree search in TNT

Runs a maximum-parsimony search in TNT (Goloboff et al. 2008; Goloboff &
Catalano 2016) using the unified search protocol validated in Smith
(2026) for continuous (Wagner), discrete (Fitch) and geometric-
morphometric landmark (Goloboff-Catalano linear displacement, following
Palci & Lee 2019) data, and returns every most-parsimonious tree found.

## Usage

``` r
RunTNTParsimony(
  data,
  type = c("continuous", "discrete", "landmark"),
  dim = NULL,
  tntPath = FindTNT(),
  workDir = tempfile("tnt"),
  timeout = 300,
  retrySeed = 7L
)
```

## Arguments

- data:

  For `type = "continuous"`: a taxa x characters numeric matrix. For
  `type = "discrete"`: a taxa x characters matrix of single-character
  state labels. For `type = "landmark"`: a named list of taxa's
  landmarks x `dim` numeric matrices (e.g. from
  [`ProcrustesAlign()`](https://phylo-pca.github.io/reference/ProcrustesAlign.md)
  or
  [`RftraAlign()`](https://phylo-pca.github.io/reference/RftraAlign.md))
  – or, to search two or more shape characters *jointly* (e.g. a
  combined 2D+3D analysis), a list of such lists, one per character,
  with `dim` given as a matching vector.

- type:

  Data type to search under; see `data`.

- dim:

  Number of spatial dimensions (required, and only used, when
  `type = "landmark"`); a vector when `data` holds multiple characters.

- tntPath:

  Path to the TNT executable; see
  [`FindTNT()`](https://phylo-pca.github.io/reference/FindTNT.md).

- workDir:

  Working directory for the (purely alphabetic, per TNT's
  filename-parsing quirk) job script and tree file; created if absent.

- timeout:

  Maximum seconds to allow the search to run before aborting (TNT
  occasionally hangs on a malformed/overflowing buffer).

- retrySeed:

  A second `rseed` to retry with, once, if the first attempt fails or
  times out.

## Value

A list with elements `trees` (a `multiPhylo` of every most-parsimonious
tree found, or `NULL` on failure) and `nTree` (the number of
most-parsimonious trees).
