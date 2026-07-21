# Shape data maximum parsimony: resistant-fit (RFTRA) alignment

As `score_native`, but using resistant-fit superimposition
([`RftraAlign()`](https://phylo-pca.github.io/reference/RftraAlign.md))
rather than native or Procrustes alignment.

## Usage

``` r
score_rftra
```

## Format

As `score_native`, with `align` fixed to `"rftra"`.

## Source

`hominin-pca-scale/R/25_shape_rftra_gen.R` + `R/26_shape_par_score.R`,
run on Hamilton.

## Note

See the note on `score_native`: a confirmatory re-run under the fuller
unified search protocol was in progress when this snapshot was taken.
