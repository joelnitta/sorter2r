# Run the SORTER2 FormatReads step

Run the SORTER2 FormatReads step

## Usage

``` r
sorter2_format_reads(
  input,
  reads_dir = getwd(),
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  dry_run = FALSE,
  verbose = 1L
)
```

## Arguments

- input:

  Input CSV/TSV file mapping original filenames to sample IDs.

- reads_dir:

  Directory containing the FASTQ files to rename.

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

Path to `reads_dir` (for use with `tar_file()`).
