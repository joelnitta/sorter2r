#' Run the SORTER2 FormatReads step
#'
#' @param input Input file path.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param working_dir Working directory for the command.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @param echo If `TRUE`, print the command before running it.
#' @return A list describing the command execution.
#' @export
sorter2_format_reads <- function(
  input,
  python = Sys.getenv(
    "SORTER2R_PYTHON",
    "python"
  ),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  working_dir = getwd(),
  dry_run = FALSE,
  echo = TRUE
) {
  input <- sorter2_require_file(input, "input")
  working_dir <- sorter2_require_dir(working_dir, "working_dir")

  args <- c("-i", basename(input))
  sorter2_run(
    script = "SORTER2_FormatReads.py",
    args = args,
    python = python,
    conda_env = conda_env,
    conda = conda,
    script_dir = script_dir,
    working_dir = working_dir,
    dry_run = dry_run,
    echo = echo
  )
}

#' Run the SORTER2 Stage1A step
#'
#' @param projname Project name.
#' @param spades Logical flag for SPAdes.
#' @param trim Logical flag for trimming.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param working_dir Working directory for the command.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @param echo If `TRUE`, print the command before running it.
#' @return A list describing the command execution.
#' @export
sorter2_stage1a <- function(
  projname,
  spades = TRUE,
  trim = TRUE,
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  working_dir = getwd(),
  dry_run = FALSE,
  echo = TRUE
) {
  if (!nzchar(trimws(projname))) {
    stop("projname must be non-empty.", call. = FALSE)
  }
  working_dir <- sorter2_require_dir(working_dir, "working_dir")

  args <- c(
    "-n",
    projname,
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
    working_dir = working_dir,
    dry_run = dry_run,
    echo = echo
  )
}

#' Run the SORTER2 Stage1B step
#'
#' @param workingdir Working directory containing the Stage1A output.
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
#' @param dry_run If `TRUE`, return the command without executing it.
#' @param echo If `TRUE`, print the command before running it.
#' @return A list describing the command execution.
#' @export
sorter2_stage1b <- function(
  workingdir,
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
  dry_run = FALSE,
  echo = TRUE
) {
  workingdir <- sorter2_with_trailing_slash(workingdir)
  ref <- sorter2_require_file(ref, "ref")

  args <- c(
    "-wd",
    workingdir,
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
}

#' Run the SORTER2 Stage2 step
#'
#' @param workingdir Working directory containing the Stage1B output.
#' @param phasequal Phase quality threshold.
#' @param aliter Alignment iterations.
#' @param indelrep Indel representation threshold.
#' @param idformat Identifier format.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @param echo If `TRUE`, print the command before running it.
#' @return A list describing the command execution.
#' @export
sorter2_stage2 <- function(
  workingdir,
  phasequal = 20,
  aliter = 1000,
  indelrep = 0.1,
  idformat = "phase",
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  dry_run = FALSE,
  echo = TRUE
) {
  workingdir <- sorter2_with_trailing_slash(workingdir)

  args <- c(
    "-wd",
    workingdir,
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
}

#' Run the SORTER2 Stage3 step
#'
#' @param workingdir Working directory containing the Stage2 output.
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
#' @return A list describing the command execution.
#' @export
sorter2_stage3 <- function(
  workingdir,
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
  workingdir <- sorter2_with_trailing_slash(workingdir)
  ref <- sorter2_require_file(ref, "ref")

  args <- c(
    "-wd",
    workingdir,
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
}

#' Run the SORTER2 Processor step
#'
#' @param workingdir Working directory containing the SORTER2 outputs.
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
#' @return A list describing the command execution.
#' @export
sorter2_processor <- function(
  workingdir,
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
  workingdir <- sorter2_with_trailing_slash(workingdir)

  args <- c(
    "-wd",
    workingdir,
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
}
