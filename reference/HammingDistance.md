# Hamming (p-)distance between discrete character rows

Proportion of characters differing between each pair of taxa; valid for
any state-space size (unlike a binary Hamming distance, no assumption is
made about the number of states a character can take).

## Usage

``` r
HammingDistance(X)
```

## Arguments

- X:

  A taxa x characters matrix of state labels.

## Value

A symmetric taxa x taxa distance matrix.
