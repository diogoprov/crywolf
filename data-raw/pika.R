# Builds the `pika` dataset from post_process.csv of the Dryad deposit of
# Goldman et al. (2023), doi:10.5061/dryad.ttdz08m2n (CC0). The file is kept here
# because the Dryad API asks for a token to download it. Its SHA-256 is the one
# Dryad publishes, a54f34a12b4efe81abfef4c643d44e75f575da10e0e9e54067524aadd4f12af2.
#
# Rows and columns are those of the deposit; only empty cells become NA.
#
# Run from the package root: Rscript data-raw/pika.R

pika <- read.csv(file.path("data-raw", "post_process.csv"), na.strings = c("NA", ""))
save(pika, file = file.path("data", "pika.rda"), compress = "xz", version = 2)
