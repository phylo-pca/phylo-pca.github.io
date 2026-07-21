# Read a RevBayes tree trace

Reads a RevBayes `.trees` log (tab-delimited, one Newick string per
generation in the last column, `[&...]` clade annotations stripped) and
discards a proportional burn-in.

## Usage

``` r
ReadRevBayesTrees(file, burnin = 0.25)
```

## Arguments

- file:

  Path to a RevBayes `.trees` file.

- burnin:

  Proportion of samples to discard from the start of the trace.

## Value

A `multiPhylo` of post-burnin trees, or `NULL` if the file has fewer
than four sampled trees.
