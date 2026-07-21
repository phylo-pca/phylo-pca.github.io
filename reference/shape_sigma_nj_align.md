# Shape data neighbour-joining: alignment comparison

Neighbour-joining accuracy on shape data under three alignments -
Procrustes
([`ProcrustesAlign()`](https://phylo-pca.github.io/reference/ProcrustesAlign.md)),
RFTRA
([`RftraAlign()`](https://phylo-pca.github.io/reference/RftraAlign.md))
and "native" (no superimposition) - across dimension and `sigma`, the
neighbour- joining half of the alignment-comparison Figure
(`fig_alignment.png`).

## Usage

``` r
shape_sigma_nj_align
```

## Format

A data frame with columns `t` (tree index), `dim`, `sigma`, `align`
(`"proc"`, `"rftra"` or `"raw"`) and `CID`.

## Source

`hominin-pca-scale/R/28_shape_nj_align.R`.
