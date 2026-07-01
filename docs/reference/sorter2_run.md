# Run a SORTER2 Python script

Run a SORTER2 Python script

## Usage

``` r
sorter2_run(
  script,
  args = character(),
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  working_dir = NULL,
  dry_run = FALSE,
  echo = FALSE,
  show_output = TRUE,
  success_codes = 0L
)
```

## Arguments

- script:

  Python script filename relative to the script directory.

- args:

  Character vector of command-line arguments.

- python:

  Python executable to use.

- conda_env:

  Optional conda environment name. If set, the command is run with
  `conda run -n conda_env`.

- conda:

  Conda executable to use when `conda_env` is set.

- script_dir:

  Directory containing the vendored SORTER2 scripts.

- working_dir:

  Optional working directory for the command.

- dry_run:

  If `TRUE`, return the command without executing it.

- echo:

  If `TRUE`, print the command before running it.

- show_output:

  If `TRUE` (default), stream script stdout/stderr as R
  [`message()`](https://rdrr.io/r/base/message.html) calls. Set to
  `FALSE` for completely silent execution.

- success_codes:

  Integer vector of exit codes treated as success. Defaults to `0L`. Use
  `c(0L, 1L)` for scripts that exit with 1 on success.

## Value

`invisible(NULL)`. Stops with an error on non-zero exit.
