# Read BirdNET selection tables

Reads the selection tables that BirdNET writes in the Raven format (one
tab-separated `.txt` file per recording or per site), with the column
names standardized as in the birdnet-validation-sampler pipeline: lower
case, without units, and with spaces replaced by underscores
(`Begin Time (s)` becomes `begin_time`).

## Usage

``` r
fp_read_birdnet(path, pattern = "\\.txt$", recursive = TRUE, tz = "UTC")
```

## Arguments

- path:

  Paths of selection tables, or of directories that hold them.

- pattern:

  Regular expression for the table files inside directories.

- recursive:

  Search directories recursively?

- tz:

  Time zone of the recording times in the file names.

## Value

A data frame with one row per detection, with the standardized columns
of the tables plus `table` (the selection table it came from), `file`
(the name of the audio file), `start` (start of the detection in the
audio file, in seconds: `file_offset` when the table has it, otherwise
`begin_time`), and `time` (recording start parsed from a
`YYYYMMDD_HHMMSS` stamp in the file name, plus `start`; `NA` when the
name has none).

## See also

[`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md),
[`fp_read_validation()`](https://provetelab.org/crywolf/reference/fp_read_validation.md)

## Examples

``` r
tab <- tempfile(fileext = ".txt")
writeLines(c(
  paste("Selection", "View", "Channel", "Begin Time (s)", "End Time (s)",
        "Begin Path", "File Offset (s)", "Species Code", "Common Name",
        "Confidence", sep = "\t"),
  paste(1, "Spectrogram 1", 1, 0, 3, "/rec/SITE1_20240101_200000.wav", 0,
        "compau", "Common Pauraque", 0.81, sep = "\t")
), tab)
fp_read_birdnet(tab)
#>   selection          view channel begin_time end_time
#> 1         1 Spectrogram 1       1          0        3
#>                       begin_path file_offset species_code     common_name
#> 1 /rec/SITE1_20240101_200000.wav           0       compau Common Pauraque
#>   confidence                                table                      file
#> 1       0.81 /tmp/Rtmpw6vysN/file1c466200f67c.txt SITE1_20240101_200000.wav
#>   start                time
#> 1     0 2024-01-01 20:00:00
```
