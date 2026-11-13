# DAGs of the models for vignettes/crywolf.Rmd, vignettes/dynamic.Rmd and
# vignettes/multispecies.Rmd. The first is adapted from the fpcalib
# compendium (manuscript/figures/fp_dag_v2.dot) to the notation of the
# vignette and the default priors of fp_priors(). Needs Graphviz (`dot`):
# the copy bundled with DiagrammeR misplaces subscripts. Run from the
# package root:
#   Rscript dev/dag/make_dag.R
system2("dot", c("-Tpng", "-Gdpi=120", file.path("dev", "dag", "crywolf_dag.dot"),
                 "-o", file.path("vignettes", "crywolf-dag.png")))
system2("dot", c("-Tpng", "-Gdpi=120", file.path("dev", "dag", "crywolf_dyn_dag.dot"),
                 "-o", file.path("vignettes", "dynamic-dag.png")))
# the multispecies model of vignettes/multispecies.Rmd, in its fs and mixed forms
ms <- readLines(file.path("dev", "dag", "crywolf_ms_dag.dot"))
for (form in c("fs", "mixed")) {
  drop <- if (form == "fs") "// mixed$" else "// fs$"
  lines <- gsub("__T__", if (form == "fs") "" else " + t<SUB>is</SUB>", ms[!grepl(drop, ms)])
  dot <- tempfile(fileext = ".dot")
  writeLines(lines, dot)
  system2("dot", c("-Tpng", "-Gdpi=120", dot,
                   "-o", file.path("vignettes", sprintf("multispecies-dag-%s.png", form))))
}
