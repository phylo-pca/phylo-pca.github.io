# Write a discrete-character Nexus file

Writes a `datatype=Standard` Nexus matrix (one symbol per cell) readable
by RevBayes' `readDiscreteCharacterData()`.

## Usage

``` r
WriteDiscreteNexus(X, file)
```

## Arguments

- X:

  A taxa x characters matrix of single-digit state labels (as produced
  by
  [`SimulateDiscreteMS()`](https://phylo-pca.github.io/reference/SimulateDiscreteMS.md)).

- file:

  Output file path.

## Value

`file`, invisibly.
