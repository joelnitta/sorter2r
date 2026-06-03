sorter2_format_reads <- function(
  input,
  python = Sys.getenv(
    "SORTER2R_PYTHON",
    "python"
  ),
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
    script_dir = script_dir,
    working_dir = working_dir,
    dry_run = dry_run,
    echo = echo
  )
}

sorter2_stage1a <- function(
  projname,
  spades = TRUE,
  trim = TRUE,
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
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
    script_dir = script_dir,
    working_dir = working_dir,
    dry_run = dry_run,
    echo = echo
  )
}

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
    script_dir = script_dir,
    dry_run = dry_run,
    echo = echo
  )
}

sorter2_stage2 <- function(
  workingdir,
  phasequal = 20,
  aliter = 1000,
  indelrep = 0.1,
  idformat = "phase",
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
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
    script_dir = script_dir,
    dry_run = dry_run,
    echo = echo
  )
}

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
    script_dir = script_dir,
    dry_run = dry_run,
    echo = echo
  )
}

sorter2_processor <- function(
  workingdir,
  repfilt = 0.50,
  majorclusters = 2,
  keepal = TRUE,
  stage2 = FALSE,
  stage3 = FALSE,
  dovcf = TRUE,
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
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
    script_dir = script_dir,
    dry_run = dry_run,
    echo = echo
  )
}
