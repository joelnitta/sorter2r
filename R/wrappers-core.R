#' Run the SORTER2 FormatReads step
#'
#' @param input Input CSV/TSV file mapping original filenames to sample IDs.
#' @param reads_dir Directory containing the FASTQ files to rename.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @param echo If `TRUE`, print the command before running it.
#' @return Path to `reads_dir` (for use with `tar_file()`).
#' @export
sorter2_format_reads <- function(
  input,
  reads_dir = getwd(),
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  dry_run = FALSE,
  echo = TRUE
) {
  input <- sorter2_require_file(input, "input")
  reads_dir <- sorter2_require_dir(reads_dir, "reads_dir")

  args <- c("-i", basename(input))
  sorter2_run(
    script = "SORTER2_FormatReads.py",
    args = args,
    python = python,
    conda_env = conda_env,
    conda = conda,
    script_dir = script_dir,
    working_dir = reads_dir,
    dry_run = dry_run,
    echo = echo
  )
  reads_dir
}

#' Run the SORTER2 Stage1A step
#'
#' @param reads_dir Directory containing the input FASTQ files.
#' @param output_dir Directory where the assembled output will be created.
#' @param spades Logical flag for SPAdes assembly.
#' @param trim Logical flag for TrimGalore trimming.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param overwrite If `TRUE`, delete `output_dir` before running.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @param echo If `TRUE`, print the command before running it.
#' @return Path to `output_dir` (for use with `tar_file()`).
#' @export
sorter2_stage1a <- function(
  reads_dir,
  output_dir,
  spades = TRUE,
  trim = TRUE,
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  overwrite = FALSE,
  dry_run = FALSE,
  echo = TRUE
) {
  reads_dir <- sorter2_require_dir(reads_dir, "reads_dir")
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (!dry_run) sorter2_check_overwrite(output_dir, overwrite)

  args <- c(
    "-o",
    output_dir,
    "-spades",
    sorter2_bool_flag(spades),
    "-trim",
    sorter2_bool_flag(trim)
  )

  sorter2_run(
    script = "SORTER2_Stage1A_TrimSPAdes.py",
    args = args,
    python = python,
    conda_env = conda_env,
    conda = conda,
    script_dir = script_dir,
    working_dir = reads_dir,
    dry_run = dry_run,
    echo = echo
  )
  output_dir
}

#' Run the SORTER2 Stage1B step
#'
#' @param input_dir Directory containing the Stage1A output
#'   (`*_assembly/` subdirs with SPAdes contigs).
#' @param output_dir Directory where `diploids/` and `diploidclusters/` will
#'   be created.
#' @param ref Reference file path.
#' @param loci Number of loci to process.
#' @param clust2id Clustering threshold.
#' @param recluster Logical flag for reclustering.
#' @param contignum Contig count threshold.
#' @param contiglen Contig length threshold.
#' @param aliter Alignment iterations.
#' @param indelrep Indel representation threshold.
#' @param idformat Identifier format.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param overwrite If `TRUE`, delete `output_dir` before running.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @param echo If `TRUE`, print the command before running it.
#' @return Path to `output_dir` (for use with `tar_file()`).
#' @export
sorter2_stage1b <- function(
  input_dir,
  output_dir,
  ref,
  loci,
  clust2id = 0.70,
  recluster = FALSE,
  contignum = 20,
  contiglen = 300,
  aliter = 1000,
  indelrep = 0.1,
  idformat = "onlysample",
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  overwrite = FALSE,
  dry_run = FALSE,
  echo = TRUE
) {
  input_dir <- sorter2_with_trailing_slash(input_dir)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  ref <- sorter2_require_file(ref, "ref")
  if (!dry_run) sorter2_check_overwrite(output_dir, overwrite)

  args <- c(
    "-wd",
    input_dir,
    "-outdir",
    output_dir,
    "-ref",
    ref,
    "-loci",
    as.character(as.integer(loci)),
    "-cl",
    as.character(clust2id),
    "-reclust",
    sorter2_bool_flag(recluster),
    "-csn",
    as.character(as.integer(contignum)),
    "-csl",
    as.character(as.integer(contiglen)),
    "-al",
    as.character(as.integer(aliter)),
    "-indel",
    as.character(indelrep),
    "-idformat",
    as.character(idformat)
  )

  sorter2_run(
    script = "SORTER2_Stage1B_AssembleOrthologs.py",
    args = args,
    python = python,
    conda_env = conda_env,
    conda = conda,
    script_dir = script_dir,
    dry_run = dry_run,
    echo = echo
  )
  output_dir
}

#' Run the SORTER2 Stage2 step
#'
#' @param input_assemblies Directory containing the Stage1A output
#'   (`*_assembly/` subdirs with trimmed FASTQ files).
#' @param input_clusters Directory containing the Stage1B output
#'   (`diploidclusters/` and `diploids/` subdirs).
#' @param output_dir Directory where `diploids_phased/` will be created.
#' @param phasequal Phase quality threshold.
#' @param aliter Alignment iterations.
#' @param indelrep Indel representation threshold.
#' @param idformat Identifier format.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param overwrite If `TRUE`, delete `output_dir` before running.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @param echo If `TRUE`, print the command before running it.
#' @return Path to `output_dir` (for use with `tar_file()`).
#' @export
sorter2_stage2 <- function(
  input_assemblies,
  input_clusters,
  output_dir,
  phasequal = 20,
  aliter = 1000,
  indelrep = 0.1,
  idformat = "phase",
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  overwrite = FALSE,
  dry_run = FALSE,
  echo = TRUE
) {
  input_assemblies <- sorter2_with_trailing_slash(input_assemblies)
  input_clusters <- sorter2_with_trailing_slash(input_clusters)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (!dry_run) sorter2_check_overwrite(output_dir, overwrite)

  args <- c(
    "-wa",
    input_assemblies,
    "-wc",
    input_clusters,
    "-outdir",
    output_dir,
    "-pq",
    as.character(as.integer(phasequal)),
    "-al",
    as.character(as.integer(aliter)),
    "-indel",
    as.character(indelrep),
    "-idformat",
    as.character(idformat)
  )

  sorter2_run(
    script = "SORTER2_Stage2_PhaseOrthologs.py",
    args = args,
    python = python,
    conda_env = conda_env,
    conda = conda,
    script_dir = script_dir,
    dry_run = dry_run,
    echo = echo
  )
  output_dir
}

#' Run the SORTER2 Stage3 step
#'
#' @param input_dir Directory containing the Stage2 output (`diploids_phased/`
#'   and `diploidclusters/`). Also expected to contain a `phaseset/`
#'   subdirectory with per-sample assembly dirs (user-created).
#' @param output_dir Directory where Stage3 results will be written. Currently
#'   must equal `input_dir`; full separation will be supported in a future
#'   update.
#' @param ref Reference file path.
#' @param loci Number of loci to process.
#' @param contigscafnum Contig/scaffold count threshold.
#' @param contigscaflen Contig/scaffold length threshold.
#' @param phasequal Phase quality threshold.
#' @param aliter Alignment iterations.
#' @param indelrep Indel representation threshold.
#' @param filterundiff Logical flag for filtering undifferentiated loci.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @param echo If `TRUE`, print the command before running it.
#' @return Path to `output_dir` (for use with `tar_file()`).
#' @export
sorter2_stage3 <- function(
  input_dir,
  output_dir = input_dir,
  ref,
  loci,
  contigscafnum = 20,
  contigscaflen = 300,
  phasequal = 20,
  aliter = 1000,
  indelrep = 0.1,
  filterundiff = FALSE,
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  dry_run = FALSE,
  echo = TRUE
) {
  input_dir <- sorter2_with_trailing_slash(input_dir)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  ref <- sorter2_require_file(ref, "ref")

  args <- c(
    "-wd",
    input_dir,
    "-outdir",
    output_dir,
    "-ref",
    ref,
    "-loci",
    as.character(as.integer(loci)),
    "-csn",
    as.character(as.integer(contigscafnum)),
    "-csl",
    as.character(as.integer(contigscaflen)),
    "-pq",
    as.character(as.integer(phasequal)),
    "-al",
    as.character(as.integer(aliter)),
    "-indel",
    as.character(indelrep),
    "-fp",
    sorter2_bool_flag(filterundiff)
  )

  sorter2_run(
    script = "SORTER2_Stage3_PhaseHybrids.py",
    args = args,
    python = python,
    conda_env = conda_env,
    conda = conda,
    script_dir = script_dir,
    dry_run = dry_run,
    echo = echo
  )
  output_dir
}

#' Run the SORTER2 Processor step
#'
#' @param input_dir Directory containing the SORTER2 pipeline outputs
#'   (`diploids/`, `diploidclusters/`, and optionally `diploids_phased/` and
#'   `phaseset/`).
#' @param output_dir Directory where processed results will be written.
#'   Currently must equal `input_dir`; full separation will be supported in a
#'   future update.
#' @param repfilt Repeat filter threshold.
#' @param majorclusters Number of major clusters.
#' @param keepal Logical flag for keeping alignments.
#' @param stage2 Logical flag for including Stage2 output.
#' @param stage3 Logical flag for including Stage3 output.
#' @param dovcf Logical flag for writing VCF output.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @param echo If `TRUE`, print the command before running it.
#' @return Path to `output_dir` (for use with `tar_file()`).
#' @export
sorter2_processor <- function(
  input_dir,
  output_dir = input_dir,
  repfilt = 0.50,
  majorclusters = 2,
  keepal = TRUE,
  stage2 = FALSE,
  stage3 = FALSE,
  dovcf = TRUE,
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  dry_run = FALSE,
  echo = TRUE
) {
  input_dir <- sorter2_with_trailing_slash(input_dir)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

  args <- c(
    "-wd",
    input_dir,
    "-outdir",
    output_dir,
    "-rep",
    as.character(repfilt),
    "-mc",
    as.character(as.integer(majorclusters)),
    "-al",
    sorter2_bool_flag(keepal),
    "-st2",
    sorter2_bool_flag(stage2),
    "-st3",
    sorter2_bool_flag(stage3),
    "-dovcf",
    sorter2_bool_flag(dovcf)
  )

  sorter2_run(
    script = "SORTER2_Processor.py",
    args = args,
    python = python,
    conda_env = conda_env,
    conda = conda,
    script_dir = script_dir,
    dry_run = dry_run,
    echo = echo
  )
  output_dir
}
