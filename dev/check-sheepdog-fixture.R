# Checks that tests/testthat/fixtures/sheepdog-export is still identical to
# the fixture of sheepdog (tests/testthat/fixtures/export-1.2 on its main
# branch), which sheepdog's contract test (test-contract.R) checks against the
# schema. Files are compared by their git blob SHA, as listed by the GitHub API,
# so nothing is downloaded unless the copies differ.
# Run from the package root: Rscript dev/check-sheepdog-fixture.R [--update]
# With --update, the files that differ are downloaded and the ones sheepdog no
# longer has are removed; then rerun the tests of fp_read_sheepdog() and
# fp_draws_data().
update <- "--update" %in% commandArgs(trailingOnly = TRUE)
local <- "tests/testthat/fixtures/sheepdog-export"
api <- paste0("https://api.github.com/repos/diogoprov/sheepdog/contents/",
              "tests/testthat/fixtures/export-1.2?ref=main")
remote <- jsonlite::fromJSON(api)
remote <- remote[remote$type == "file", c("name", "sha", "download_url")]
files <- list.files(local)
sha <- if (length(files)) system2("git", c("hash-object", file.path(local, files)),
                                  stdout = TRUE) else character()
changed <- remote$name[!remote$name %in% files |
                         sha[match(remote$name, files)] != remote$sha]
gone <- setdiff(files, remote$name)
if (!length(changed) && !length(gone)) {
  message("The fixture is identical to sheepdog's (", length(files), " files).")
  quit(status = 0)
}
if (length(changed)) message("Differ from sheepdog: ", toString(changed))
if (length(gone)) message("Not in sheepdog: ", toString(gone))
if (!update) {
  message("Rerun with --update to copy sheepdog's fixture.")
  quit(status = 1)
}
for (f in changed) {
  utils::download.file(remote$download_url[remote$name == f], file.path(local, f),
                       mode = "wb", quiet = TRUE)
}
unlink(file.path(local, gone))
message("Updated. Rerun devtools::test(filter = \"read_sheepdog|draws_data\").")
