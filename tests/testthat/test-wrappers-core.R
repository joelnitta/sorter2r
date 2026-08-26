test_that("sorter2_bool_flag normalizes logical and T/F inputs", {
  expect_identical(sorter2_bool_flag(TRUE), "T")
  expect_identical(sorter2_bool_flag(FALSE), "F")
  expect_identical(sorter2_bool_flag("t"), "T")
  expect_identical(sorter2_bool_flag(" F "), "F")
  expect_error(sorter2_bool_flag("yes"), "Flag values must be logical")
})

test_that("sorter2_require_file and sorter2_require_dir validate paths", {
  file_path <- tempfile("sorter2r-file-")
  dir_path <- tempfile("sorter2r-dir-")
  file.create(file_path)
  dir.create(dir_path)

  expect_identical(
    sorter2_require_file(file_path),
    normalizePath(file_path, winslash = "/", mustWork = TRUE)
  )
  expect_identical(
    sorter2_require_dir(dir_path),
    normalizePath(dir_path, winslash = "/", mustWork = TRUE)
  )
  expect_error(
    sorter2_require_file(file.path(tempdir(), "missing-file")),
    "does not exist"
  )
  expect_error(
    sorter2_require_dir(file.path(tempdir(), "missing-dir")),
    "does not exist"
  )
})

test_that("sorter2_with_trailing_slash appends one slash", {
  dir_path <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
  expect_identical(sorter2_with_trailing_slash(dir_path), paste0(dir_path, "/"))
})

test_that("sorter2_run serializes dry-run commands", {
  script_dir <- tempfile("sorter2r-scripts-")
  dir.create(script_dir)
  file.create(file.path(script_dir, "SORTER2_Stage1A_TrimSPAdes.py"))

  res <- sorter2_run(
    script = "SORTER2_Stage1A_TrimSPAdes.py",
    args = c("-o", "/tmp/out", "-spades", "T", "-trim", "F"),
    python = "python",
    script_dir = script_dir,
    working_dir = tempdir(),
    dry_run = TRUE,
    echo = FALSE
  )

  expect_true(res$success)
  expect_identical(res$status, 0L)
  expect_true(res$dry_run)
  expect_match(res$command, "python")
  expect_match(res$command, "SORTER2_Stage1A_TrimSPAdes.py")
  expect_match(res$command, "'-o'")
})

test_that("sorter2_run can serialize conda dry-run commands", {
  script_dir <- tempfile("sorter2r-scripts-")
  dir.create(script_dir)
  file.create(file.path(script_dir, "SORTER2_Stage1A_TrimSPAdes.py"))

  res <- sorter2_run(
    script = "SORTER2_Stage1A_TrimSPAdes.py",
    args = c("-o", "/tmp/out"),
    python = "python",
    conda_env = "sorter2",
    conda = "conda",
    script_dir = script_dir,
    working_dir = tempdir(),
    dry_run = TRUE,
    echo = FALSE
  )

  expect_true(res$success)
  expect_identical(res$status, 0L)
  expect_true(res$dry_run)
  expect_match(res$command, "conda")
  expect_match(res$command, "'run' '-n' 'sorter2' 'python'")
  expect_match(res$command, "SORTER2_Stage1A_TrimSPAdes.py")
})

test_that("sorter2_run treats blank conda_env as direct python", {
  script_dir <- tempfile("sorter2r-scripts-")
  dir.create(script_dir)
  file.create(file.path(script_dir, "SORTER2_Stage1A_TrimSPAdes.py"))

  res <- sorter2_run(
    script = "SORTER2_Stage1A_TrimSPAdes.py",
    python = "python",
    conda_env = NULL,
    script_dir = script_dir,
    dry_run = TRUE,
    echo = FALSE
  )

  expect_match(res$command, "python")
  expect_no_match(res$command, "'run' '-n'")
})

test_that("wrapper dry-runs build expected command arguments", {
  script_dir <- tempfile("sorter2r-scripts-")
  dir.create(script_dir)

  scripts <- c(
    "SORTER2_FormatReads.py",
    "SORTER2_Stage1A_TrimSPAdes.py",
    "SORTER2_Stage1B_AssembleOrthologs.py",
    "SORTER2_Stage2_PhaseOrthologs.py",
    "SORTER2_Stage3_PhaseHybrids.py",
    "SORTER2_Processor.py",
    "SORTER2_ProgenitorProcessor.py"
  )
  for (script in scripts) {
    file.create(file.path(script_dir, script))
  }

  input_file <- tempfile("sorter2r-input-", fileext = ".fastq")
  file.create(input_file)
  ref_file <- tempfile("sorter2r-ref-", fileext = ".fasta")
  file.create(ref_file)
  reads_dir <- tempfile("sorter2r-reads-")
  dir.create(reads_dir)
  input_dir <- tempfile("sorter2r-input-")
  dir.create(input_dir)
  output_dir <- tempfile("sorter2r-output-")

  format_res <- sorter2_format_reads(
    input = input_file,
    reads_dir = reads_dir,
    script_dir = script_dir,
    dry_run = TRUE
  )
  expect_match(format_res$command, "SORTER2_FormatReads.py")
  expect_match(format_res$command, basename(input_file))

  stage1a_res <- sorter2_stage1a(
    reads_dir = reads_dir,
    output_dir = output_dir,
    script_dir = script_dir,
    dry_run = TRUE
  )
  expect_match(stage1a_res$command, "SORTER2_Stage1A_TrimSPAdes.py")
  expect_match(stage1a_res$command, "'-o'")
  expect_match(stage1a_res$command, "'-spades' 'T'")
  expect_match(stage1a_res$command, "'-trim' 'T'")

  stage1b_res <- sorter2_stage1b(
    input_dir = input_dir,
    output_dir = output_dir,
    ref = ref_file,
    loci = 10,
    script_dir = script_dir,
    dry_run = TRUE
  )
  expect_match(stage1b_res$command, "SORTER2_Stage1B_AssembleOrthologs.py")
  expect_match(stage1b_res$command, "'-wd'")
  expect_match(stage1b_res$command, "'-outdir'")
  expect_match(stage1b_res$command, "'-ref'")
  expect_match(stage1b_res$command, "'-loci' '10'")

  stage2_res <- sorter2_stage2(
    input_assemblies = input_dir,
    input_clusters = input_dir,
    output_dir = output_dir,
    script_dir = script_dir,
    dry_run = TRUE
  )
  expect_match(stage2_res$command, "SORTER2_Stage2_PhaseOrthologs.py")
  expect_match(stage2_res$command, "'-wa'")
  expect_match(stage2_res$command, "'-wc'")
  expect_match(stage2_res$command, "'-pq' '20'")

  stage3_res <- sorter2_stage3(
    input_phased     = input_dir,
    input_assemblies = input_dir,
    output_dir       = output_dir,
    ref              = ref_file,
    loci             = 10,
    script_dir       = script_dir,
    dry_run          = TRUE
  )
  expect_match(stage3_res$command, "SORTER2_Stage3_PhaseHybrids.py")
  expect_match(stage3_res$command, "'-outdir'")
  expect_match(stage3_res$command, "'-fp' 'F'")
  expect_match(stage3_res$command, "'-specimen' 'F'")

  processor_res <- sorter2_processor(
    input_dir = input_dir,
    output_dir = output_dir,
    script_dir = script_dir,
    dry_run = TRUE
  )
  expect_match(processor_res$command, "SORTER2_Processor.py")
  expect_match(processor_res$command, "'-outdir'")
  expect_match(processor_res$command, "'-dovcf' 'T'")

  reads_dir <- tempfile("sorter2r-reads-")
  dir.create(reads_dir)
  processor_reads_res <- sorter2_processor(
    input_dir = input_dir,
    output_dir = output_dir,
    reads_dir = reads_dir,
    script_dir = script_dir,
    dry_run = TRUE
  )
  expect_match(processor_reads_res$command, "'-reads'")

  progenitor_res <- sorter2_progenitor_processor(
    input_dir  = input_dir,
    output_dir = output_dir,
    mapfile    = ref_file,
    outgroups  = ref_file,
    script_dir = script_dir,
    dry_run    = TRUE
  )
  expect_match(progenitor_res$command, "SORTER2_ProgenitorProcessor.py")
  expect_match(progenitor_res$command, "'-wd'")
  expect_match(progenitor_res$command, "'-indir'")
  expect_match(progenitor_res$command, "'-map'")
  expect_match(progenitor_res$command, "'-dif' 'F'")
})

test_that("sorter2_stage3 tag_specimen dry-run builds expected args", {
  script_dir <- tempfile("sorter2r-scripts-")
  dir.create(script_dir)
  file.create(file.path(script_dir, "SORTER2_Stage3_PhaseHybrids.py"))

  input_dir  <- tempfile("sorter2r-input-")
  output_dir <- tempfile("sorter2r-output-")
  ref_file   <- tempfile("sorter2r-ref-", fileext = ".fasta")
  dir.create(input_dir)
  file.create(ref_file)

  res <- sorter2_stage3(
    input_phased     = input_dir,
    input_assemblies = input_dir,
    output_dir       = output_dir,
    ref              = ref_file,
    loci             = 10,
    tag_specimen     = TRUE,
    script_dir       = script_dir,
    dry_run          = TRUE
  )

  expect_match(res$command, "'-specimen' 'T'")
})

test_that("sorter2_progenitor_processor dry-run builds expected args", {
  script_dir <- tempfile("sorter2r-scripts-")
  dir.create(script_dir)
  file.create(file.path(script_dir, "SORTER2_ProgenitorProcessor.py"))

  input_dir  <- tempfile("sorter2r-input-")
  output_dir <- tempfile("sorter2r-output-")
  dir.create(input_dir)
  mapfile   <- tempfile("map-", fileext = ".csv")
  outgroups <- tempfile("og-",  fileext = ".txt")
  file.create(mapfile)
  file.create(outgroups)

  res <- sorter2_progenitor_processor(
    input_dir    = input_dir,
    output_dir   = output_dir,
    mapfile      = mapfile,
    outgroups    = outgroups,
    minseq       = 20L,
    filterundiff = TRUE,
    script_dir   = script_dir,
    dry_run      = TRUE
  )

  expect_true(res$dry_run)
  expect_match(res$command, "SORTER2_ProgenitorProcessor.py")
  expect_match(res$command, "'-indir'")
  expect_match(res$command, "'-min' '20'")
  expect_match(res$command, "'-dif' 'T'")
})

test_that("sorter2_stage3 warns and returns output_dir when hybrid_samples is empty", {
  script_dir <- tempfile("sorter2r-scripts-")
  dir.create(script_dir)
  file.create(file.path(script_dir, "SORTER2_Stage3_PhaseHybrids.py"))

  ref_file <- tempfile("ref-", fileext = ".fasta")
  file.create(ref_file)
  phaseset_dir <- tempfile("sorter2r-phaseset-")
  dir.create(phaseset_dir)
  output_dir <- tempfile("sorter2r-output-")

  expect_warning(
    result <- sorter2_stage3(
      input_phased     = tempdir(),
      input_assemblies = tempdir(),
      output_dir       = output_dir,
      ref              = ref_file,
      loci             = 10,
      phaseset_dir     = phaseset_dir,
      hybrid_samples   = character(0),
      script_dir       = script_dir
    ),
    "hybrid_samples is empty"
  )
  expect_identical(result, output_dir)
})

test_that("wrapper dry-runs can use a conda environment", {
  script_dir <- tempfile("sorter2r-scripts-")
  dir.create(script_dir)
  file.create(file.path(script_dir, "SORTER2_Stage1A_TrimSPAdes.py"))

  reads_dir <- tempfile("sorter2r-reads-")
  dir.create(reads_dir)
  output_dir <- tempfile("sorter2r-output-")

  stage1a_res <- sorter2_stage1a(
    reads_dir = reads_dir,
    output_dir = output_dir,
    conda_env = "sorter2",
    script_dir = script_dir,
    dry_run = TRUE
  )

  expect_match(stage1a_res$command, "conda")
  expect_match(stage1a_res$command, "'run' '-n' 'sorter2' 'python'")
  expect_match(stage1a_res$command, "'-o'")
})
