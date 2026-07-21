# Read landmark configurations from a shape-simulator TSV

Reads the long-format landmark table written by the LDDMM shape
simulator of Raskin et al. (2026) (one row per landmark per tip, columns
`label`, an index, then one column per spatial dimension) and reshapes
it into one matrix per tip.

## Usage

``` r
ReadShapeConfigs(file, dim)
```

## Arguments

- file:

  Path to a `*_nodeShapes.tsv` file.

- dim:

  Number of spatial dimensions (2 or 3).

## Value

A named list (names = tip labels) of landmarks x `dim` numeric matrices.
