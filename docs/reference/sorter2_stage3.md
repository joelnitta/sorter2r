# Run the SORTER2 Stage3 step

Run the SORTER2 Stage3 step

## Usage

``` r
sorter2_stage3(
  input_phased,
  input_assemblies,
  output_dir,
  phaseset_dir = NULL,
  hybrid_samples = NULL,
  reads_dir = NULL,
  ref,
  loci,
  contigscafnum = 20,
  contigscaflen = 300,
  phasequal = 20,
  aliter = 1000,
  indelrep = 0.1,
  filterundiff = FALSE,
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

- input_phased:

  Directory containing the Stage2 output. Must contain a
  `diploids_phased/` subdirectory with the phased diploid sequences.

- input_assemblies:

  Directory containing the Stage1A output. Used to identify diploid
  sample names from `*_assembly/` subdirectories.

- output_dir:

  Directory where Stage3 results will be written.

- phaseset_dir:

  Optional directory containing hybrid sample `*_assembly/`
  subdirectories (SPAdes output). When provided, the wrapper copies
  those dirs into `output_dir/phaseset/` before running. When `NULL`
  (the default), `output_dir/phaseset/` must already exist and be
  populated by the user before calling this function.

- hybrid_samples:

  Optional character vector of sample base names (e.g.
  `"Iimura46_cgrande"`) used to filter which `*_assembly/` dirs are
  copied from `phaseset_dir`. When `NULL`, all `*_assembly/` dirs are
  copied.

- reads_dir:

  Optional directory containing raw FASTQ files (`*_R1.fastq` /
  `*_R2.fastq`). Required when Stage 1A was run with `trim = FALSE`.
  When `NULL` (the default), Stage3 expects Trim Galore output
  (`*_R1_val_1.fq`) inside each `phaseset/*_assembly/` subdirectory.

- ref:

  Reference file path.

- loci:

  Number of loci to process.

- contigscafnum:

  Contig/scaffold count threshold.

- contigscaflen:

  Contig/scaffold length threshold.

- phasequal:

  Phase quality threshold.

- aliter:

  Alignment iterations.

- indelrep:

  Indel representation threshold.

- filterundiff:

  Logical flag for filtering undifferentiated loci.

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

  If `TRUE`, remove Stage3 outputs (but preserve `phaseset/`) before
  re-running. When `phaseset_dir` is also set, `phaseset/` is fully
  refreshed from the source.

- dry_run:

  If `TRUE`, return the command without executing it.

## Value

Path to `output_dir` (for use with `tar_file()`).
