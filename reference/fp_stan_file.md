# Path to the Stan program

Path to the Stan program

## Usage

``` r
fp_stan_file()
```

## Value

The path to the Stan file shipped with the package, for users who want
to read the model or modify a copy of it. The likelihood of the
observations of each site is in `fp_emissions.stanfunctions`, in the
same folder, which the program includes.

## Examples

``` r
fp_stan_file()
#> [1] "/home/runner/work/_temp/Library/crywolf/bin/stan/fp_occu.stan"
```
