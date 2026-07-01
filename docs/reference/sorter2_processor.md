# Run the SORTER2 Processor step

Run the SORTER2 Processor step

## Usage

``` r
sorter2_processor(
  input_dir,
  output_dir = input_dir,
  repfilt = 0.5,
  majorclusters = 2,
  keepal = TRUE,
  stage2 = FALSE,
  stage3 = FALSE,
  dovcf = TRUE,
  reads_dir = NULL,
  success_codes = 0L,
  overwrite = FALSE,
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  dry_run = FALSE,
  verbose = 1L
)
```

## Arguments

- input_dir:

  Directory containing the SORTER2 pipeline outputs (`diploids/`,
  `diploidclusters/`, and optionally `diploids_phased/` and
  `phaseset/`).

- output_dir:

  Directory where processed results will be written.

- repfilt:

  Repeat filter threshold.

- majorclusters:

  Number of major clusters.

- keepal:

  Logical flag for keeping alignments.

- stage2:

  Logical flag for including Stage2 output.

- stage3:

  Logical flag for including Stage3 output.

- dovcf:

  Logical flag for writing VCF output.

- reads_dir:

  Optional directory containing raw FASTQ files (`*_R1.fastq` /
  `*_R2.fastq`). Required when Stage 1A was run with `trim = FALSE`.
  When `NULL` (the default), the Processor expects Trim Galore output
  (`*_R1_val_1.fq`) inside each `*_assembly/` subdirectory.

- success_codes:

  Integer vector of exit codes treated as success. Default `0L`. Pass
  `c(0L, 1L)` if outputs are known to be correct and the script is
  exiting 1 on a re-run (e.g. interrupted pipeline).

- overwrite:

  If `TRUE`, delete `output_dir` before running.

- python:

  Python executable to use.

- conda_env:

  Optional conda environment name.

- conda:

  Conda executable to use when `conda_env` is set.

- script_dir:

  Directory containing the vendored SORTER2 scripts.

- dry_run:

  If `TRUE`, return the command without executing it.

- verbose:

  Integer verbosity level. `1L` (default) shows Python-level section
  headers; `0L` is fully silent; `2L` echoes the command, enables
  per-item debug output (`-v`), and shows subprocess tool output.

## Value

Path to `output_dir` (for use with `tar_file()`).
