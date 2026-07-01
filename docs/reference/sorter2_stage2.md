# Run the SORTER2 Stage2 step

Run the SORTER2 Stage2 step

## Usage

``` r
sorter2_stage2(
  input_assemblies,
  input_clusters,
  output_dir,
  reads_dir = NULL,
  phasequal = 20,
  aliter = 1000,
  indelrep = 0.1,
  idformat = "phase",
  threads = 1L,
  verbose = 1L,
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  overwrite = FALSE,
  dry_run = FALSE
)
```

## Arguments

- input_assemblies:

  Directory containing the Stage1A output (`*_assembly/` subdirs with
  SPAdes contigs).

- input_clusters:

  Directory containing the Stage1B output (`diploidclusters/` and
  `diploids/` subdirs).

- output_dir:

  Directory where `diploids_phased/` will be created.

- reads_dir:

  Optional directory containing raw FASTQ files (`*_R1.fastq` /
  `*_R2.fastq`). Required when Stage 1A was run with `trim = FALSE`.
  When `NULL` (the default), Stage2 expects Trim Galore output
  (`*_R1_val_1.fq`) inside each `*_assembly/` subdirectory.

- phasequal:

  Phase quality threshold.

- aliter:

  Alignment iterations.

- indelrep:

  Indel representation threshold.

- idformat:

  Identifier format.

- threads:

  Number of parallel worker processes for per-sample phasing and
  per-locus alignment. Default `1L` (serial).

- verbose:

  Integer verbosity level. `1L` (default) shows Python-level section
  headers; `0L` is fully silent; `2L` echoes the command, enables
  per-item debug output (`-v`), and shows subprocess tool output.

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

## Value

Path to `output_dir` (for use with `tar_file()`).
