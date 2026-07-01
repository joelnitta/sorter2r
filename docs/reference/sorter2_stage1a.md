# Run the SORTER2 Stage1A step

Run the SORTER2 Stage1A step

## Usage

``` r
sorter2_stage1a(
  reads_dir,
  output_dir,
  spades = TRUE,
  trim = TRUE,
  clean_spades_tmp = FALSE,
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  overwrite = FALSE,
  dry_run = FALSE,
  verbose = 1L
)
```

## Arguments

- reads_dir:

  Directory containing the input FASTQ files.

- output_dir:

  Directory where the assembled output will be created.

- spades:

  Logical flag for SPAdes assembly.

- trim:

  Logical flag for TrimGalore trimming.

- clean_spades_tmp:

  If `TRUE`, delete SPAdes intermediate K-mer graph directories (`K21/`,
  `K33/`, etc.) after each sample assembles. Saves substantial disk
  space at the cost of preventing SPAdes restarts. Default `FALSE`.

- python:

  Python executable to use.

- conda_env:

  Optional conda environment name.

- conda:

  Conda executable to use when `conda_env` is set.

- script_dir:

  Directory containing the vendored SORTER2 scripts.

- overwrite:

  If `TRUE`, delete `output_dir` before running.

- dry_run:

  If `TRUE`, return the command without executing it.

- verbose:

  Integer verbosity level. `1L` (default) shows Python-level section
  headers; `0L` is fully silent; `2L` echoes the command, enables
  per-item debug output (`-v`), and shows subprocess tool output.

## Value

Path to `output_dir` (for use with `tar_file()`).
