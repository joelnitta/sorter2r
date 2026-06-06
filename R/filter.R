#' Copy ingroup assembly folders from Stage 1A output, excluding outgroups
#'
#' SORTER2 Stages 1B onward should be run on ingroup samples only; outgroup
#' assemblies inflate allele counts and ADMIXTURE complexity. This function
#' creates a new directory containing only the ingroup `*_assembly` folders
#' from a Stage 1A output directory, leaving the original untouched.
#'
#' @param input_dir Directory containing Stage 1A output (the `*_assembly`
#'   subdirectories produced by [sorter2_stage1a()]).
#' @param output_dir Directory to create containing only ingroup assemblies.
#'   Must not already exist unless `overwrite = TRUE`.
#' @param exclude_samples Character vector of sample names to exclude
#'   (outgroups). Names should match the `{sample}_assembly` folder prefix —
#'   for example, `"Wade3565_cintermedium"` excludes
#'   `Wade3565_cintermedium_assembly/`.
#' @param overwrite If `TRUE`, delete `output_dir` before running.
#' @return Path to `output_dir` (for use with `tar_file()`).
#' @export
sorter2_filter_assemblies <- function(
  input_dir,
  output_dir,
  exclude_samples = character(),
  overwrite = FALSE
) {
  input_dir  <- sorter2_require_dir(input_dir, "input_dir")
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (!startsWith(output_dir, "/")) {
    output_dir <- file.path(getwd(), output_dir)
  }
  sorter2_check_overwrite(output_dir, overwrite)

  all_dirs     <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
  assembly_dirs <- all_dirs[grepl("_assembly$", basename(all_dirs))]
  sample_names  <- sub("_assembly$", "", basename(assembly_dirs))
  keep          <- !sample_names %in% exclude_samples

  n_excluded <- sum(!keep)
  n_kept     <- sum(keep)
  if (n_excluded == 0L) {
    warning(
      "No assembly folders were excluded. ",
      "Check that exclude_samples names match folder prefixes.",
      call. = FALSE
    )
  }
  message(
    sprintf(
      "Keeping %d ingroup assemblies, excluding %d outgroup(s): %s",
      n_kept,
      n_excluded,
      paste(sample_names[!keep], collapse = ", ")
    )
  )

  dir.create(output_dir, recursive = TRUE)
  for (d in assembly_dirs[keep]) {
    file.copy(d, output_dir, recursive = TRUE)
  }

  output_dir
}
