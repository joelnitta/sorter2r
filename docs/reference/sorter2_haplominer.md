# Run the SORTER2 Hapl-O-Miner step

Mines haploid organellar genomes (e.g. chloroplast) from off-target
by-catch reads in a SORTER2 Stage 1A dataset. Reads from each
`*_assembly/` subdirectory are mapped to a user-supplied haploid
organellar reference; per-sample consensus sequences are generated and
filtered by minimum coverage and read-depth thresholds.

## Usage

``` r
sorter2_haplominer(
  input_dir,
  output_dir,
  organellar_ref,
  coverage = 50,
  depth = 5,
  reads_dir = NULL,
  threads = 1L,
  clean_workfiles = FALSE,
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

  Directory containing Stage 1A `*_assembly/` subdirectories.

- output_dir:

  Directory where all Hapl-O-Miner output will be written.

- organellar_ref:

  Path to the haploid organellar reference FASTA (e.g. a chloroplast
  genome).

- coverage:

  Minimum reference coverage percentage threshold. Samples below this
  value are excluded from the final alignment.

- depth:

  Minimum mean read-depth threshold. Samples below this value are
  excluded from the final alignment.

- reads_dir:

  Optional directory containing raw FASTQ files (`*_R1.fastq` /
  `*_R2.fastq`). Required when Stage 1A was run with `trim = FALSE`.
  When `NULL` (the default), Hapl-O-Miner expects Trim Galore output
  (`*_R1_val_1.fq`) inside each `*_assembly/` subdirectory.

- threads:

  Number of parallel worker processes for per-sample mapping. Default
  `1L` (serial).

- clean_workfiles:

  If `TRUE`, delete the intermediate workfiles directory (BAMs, raw
  consensus FASTAs, per-sample stat files) after the run completes.
  Saves substantial disk space; the compiled `readstats_cp.csv` and
  final consensus FASTAs are unaffected. Default `FALSE`.

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
  headers; `0L` is fully silent; `2L` echoes the command and enables
  per-item debug output (`-v`).

## Value

Path to `output_dir` (for use with `tar_file()`).
