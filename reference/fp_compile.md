# Compile (or load) the Stan model

Compiles the false-positive occupancy model the first time it is needed
and keeps the executable in the user cache directory
(`tools::R_user_dir("crywolf", "cache")`), so later sessions load it
without recompiling.
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

## Examples

``` r
if (FALSE) { # \dontrun{
fp_compile()
} # }
```
