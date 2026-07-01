# Run the SORTER2 Stage1B step

Run the SORTER2 Stage1B step

## Usage

``` r
sorter2_stage1b(
  input_dir,
  output_dir,
  ref,
  loci,
  clust2id = 0.7,
  recluster = FALSE,
  contignum = 20,
  contiglen = 300,
  aliter = 1000,
  indelrep = 0.1,
  idformat = "onlysample",
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

- input_dir:

  Directory containing the Stage1A output (`*_assembly/` subdirs with
  SPAdes contigs).

- output_dir:

  Directory where `diploids/` and `diploidclusters/` will be created.

- ref:

  Reference file path.

- loci:

  Number of loci to process.

- clust2id:

  Clustering threshold.

- recluster:

  Logical flag for reclustering.

- contignum:

  Contig count threshold.

- contiglen:

  Contig length threshold.

- aliter:

  Alignment iterations.

- indelrep:

  Indel representation threshold.

- idformat:

  Identifier format.

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
