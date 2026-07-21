# State-space size distribution of the Mongle (2023) character matrix

Reads a Nexus character matrix and tabulates the size of each
character's state space (`max(observed state) + 1`, following RevBayes'
`setNumStatesPartition()` convention; polymorphisms such as `1/2` are
split before taking the maximum). Used to draw realistic state-space
sizes for
[`SimulateDiscreteMS()`](https://phylo-pca.github.io/reference/SimulateDiscreteMS.md).

## Usage

``` r
MongleStateProportions(file)
```

## Arguments

- file:

  Path to a Nexus-formatted discrete character matrix.

## Value

A named integer vector (a `table`): counts of characters at each
observed state-space size `k`.
