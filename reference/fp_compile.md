# Compile (or load) the Stan model

Loads the false-positive occupancy model compiled when the package was
installed from source. A binary package, or an installation without
CmdStan, comes without it: the model is then compiled the first time it
is needed and kept in the user cache directory
(`tools::R_user_dir("crywolf", "cache")`), so later sessions load it
without recompiling. The cache holds only the Stan program and its
executable, and a new version of the program replaces them. Under
`R CMD check` the cache is in
[`tempdir()`](https://rdrr.io/r/base/tempfile.html), so nothing is
written to the user's filespace.
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) calls
this function itself; call it directly to compile ahead of time, for
example before launching parallel workers in
[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md).

## Usage

``` r
fp_compile(force = FALSE, quiet = TRUE)
```

## Arguments

- force:

  Logical. Recompile even if a current executable exists.

- quiet:

  Logical. Suppress compiler output.

## Value

A `CmdStanModel` object (invisibly).

## Details

The model is compiled with threading enabled (`stan_threads`), so that
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) can
split the likelihood of each chain among several threads (argument
`threads_per_chain`).

## Examples

``` r
# \donttest{
fp_compile()
# }
```
