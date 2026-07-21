# Shape data maximum parsimony: native/Procrustes alignment

Maximum-parsimony accuracy on shape (landmark) data (TNT's landmark data
type; the linear-displacement method of Goloboff & Catalano, following
Palci & Lee 2019; see
[`RunTNTParsimony()`](https://phylo-pca.github.io/reference/RunTNTParsimony.md)),
comparing the "native" (TNT's own internal linear superimposition, no
prior alignment) and Procrustes-aligned
([`ProcrustesAlign()`](https://phylo-pca.github.io/reference/ProcrustesAlign.md))
treatments, across dimension and `sigma`.

## Usage

``` r
data(score_native)
```

## Format

A data frame with columns `idx` (job index), `t0` (tree index, 0-based),
`dim`, `sigma`, `align` (`"proc"` or `"raw"`), `CID`, `RF`, `nmpt`
(most-parsimonious trees found) and `ksel` (max-entropy subsample size
used).

## Source

`hominin-pca-scale/R/24_shape_par_gen.R` + `R/26_shape_par_score.R`, run
on Hamilton.

## Note

A re-run under the fuller unified search protocol (adding
`lmark xthreads;` and a larger sectorial-search buffer; see
[`RunTNTParsimony()`](https://phylo-pca.github.io/reference/RunTNTParsimony.md))
was in progress on Hamilton at the time this snapshot was taken. Prior
validation showed the search was already at the global optimum, so the
re-run is expected to confirm rather than change these values, but this
file will be refreshed if it does.
