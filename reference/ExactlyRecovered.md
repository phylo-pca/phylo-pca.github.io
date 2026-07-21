# Exact-match tree recovery

Whether a reconstructed tree exactly matches the true generative tree
(i.e. a clustering information distance of exactly zero). Reported
throughout Smith (2026) alongside continuous accuracy metrics, but
interpreted with care: for 21-tip unrooted binary trees there are
`(2 x 21 - 5)!! ~= 8.2 x 10^21` possible topologies, so exact-match is a
vanishingly strict, non-metric criterion, not a graded measure of
accuracy.

## Usage

``` r
ExactlyRecovered(tree, truth)
```

## Arguments

- tree:

  The reconstructed tree.

- truth:

  The true (generative) tree.

## Value

A logical.
