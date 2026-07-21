# Write a continuous-character Nexus file

Writes a `datatype=Continuous` Nexus matrix readable by RevBayes'
`readContinuousCharacterData()`.

## Usage

``` r
WriteContinuousNexus(M, file)
```

## Arguments

- M:

  A taxa x characters numeric matrix, row-named by taxon.

- file:

  Output file path.

## Value

`file`, invisibly.
