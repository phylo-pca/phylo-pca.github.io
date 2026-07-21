# Combining 2D and 3D shape data (neighbour-joining)

Neighbour-joining accuracy from 2D shape data, 3D shape data, and the
two combined, across the landmark-integration parameter `sigma`, on a
shared set of trees so the three arms are paired. `cid2`/`cid3`/`cidB`
are scored on raw Euclidean distance; `cid2sq`/`cid3sq`/`cidBsq` on the
additive (squared-Euclidean) distance, with the combined arm computed as
`nj(d2^2 + d3^2)` (the correct way two independent additive distances
fuse).

## Usage

``` r
shape_combine_nj
```

## Format

A data frame with columns `t` (tree index), `sigma`,
`cid2`/`cid3`/`cidB` (raw-distance CID for 2D/3D/combined),
`cid2sq`/`cid3sq`/`cidBsq` (additive-distance CID) and `scale_ratio`
(relative scale of the 2D vs 3D coordinate sets before combining).

## Source

`hominin-pca-scale/R/22_shape_combine_nj.R`, from LDDMM shape
simulations aligned by
[`ProcrustesAlign()`](https://phylo-pca.github.io/reference/ProcrustesAlign.md).
