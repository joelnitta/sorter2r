# Run from the package root to refresh the pre-computed vignette.
# Requires a working SORTER2 conda environment and the mini dataset.
# First run: ~1 hour (SPAdes + full pipeline).
# Subsequent runs: fast (targets caching).
knitr::knit(
  "vignettes/targets-pipeline.Rmd.orig",
  output = "vignettes/targets-pipeline.Rmd"
)
