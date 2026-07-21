# Locate the RevBayes executable

Looks for a RevBayes installation, in order: the `phyloPCA.rb` option,
the `RB_PATH` environment variable, then `rb`/`rb.exe` on the system
`PATH`.

## Usage

``` r
FindRevBayes()
```

## Value

Path to the RevBayes executable.
