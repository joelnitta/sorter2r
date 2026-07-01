# Process progenitor sequences from Stage 3 output

Runs `SORTER2_ProgenitorProcessor.py`. Collapses progenitor IDs to
user-defined clades, removes outgroup sequences, filters by minimum
sequence count, removes duplicate IDs, and optionally retains only loci
with at least two differentiated progenitor clades.

## Usage

``` r
sorter2_progenitor_processor(
  input_dir,
  output_dir,
  mapfile,
  outgroups,
  minseq = 1L,
  filterundiff = FALSE,
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

  Directory containing Stage 3 FASTA output (read-only).

- output_dir:

  Directory to write processed output into. Created if it does not
  exist.

- mapfile:

  Path to a CSV file with three columns: `hybrid`, `progenitor`,
  `clade`.

- outgroups:

  Path to a plain-text file listing outgroup IDs to remove (one per
  line).

- minseq:

  Minimum number of phased sequence pairs required per progenitor set
  across all loci. The Python script doubles this value internally
  before applying the filter. Default `1L`.

- filterundiff:

  If `TRUE`, discard loci where a hybrid sample has fewer than two
  distinct progenitor clades represented. Default `FALSE`.

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
  headers; `0L` is fully silent; `2L` echoes the command and shows
  subprocess tool output.

## Value

Path to `output_dir` (for use with `tar_file()`).

## Details

The Stage 3 `input_dir` is never modified, so its targets hash remains
stable.
