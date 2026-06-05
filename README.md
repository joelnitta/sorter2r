# sorter2r

`sorter2r` is an R package that wraps core SORTER2 Python scripts.
This initial implementation is thin by design: it validates arguments,
constructs CLI calls, and executes the original Python entry points.

## Current V1 scope

- `sorter2_format_reads()`
- `sorter2_stage1a()`
- `sorter2_stage1b()`
- `sorter2_stage2()`
- `sorter2_stage3()`
- `sorter2_processor()`

## Typical workflow

Each function takes explicit input and output directory arguments, forming
a linear pipeline where each step's output feeds the next. Stage 2 takes
two inputs (Stage 1A and Stage 1B outputs); Stage 3 also takes two inputs
(Stage 2 and Stage 1A outputs).

```r
library(sorter2r)

# Optional: rename raw FASTQ files according to a sample map CSV
sorter2_format_reads(
  input     = "samples.csv",
  reads_dir = "raw_reads/"
)

# Stage 1A — trim (optional) and/or assemble contigs with SPAdes.
# Set trim = FALSE if reads are already cleaned; pass reads_dir to
# sorter2_stage2() in that case so Stage 2 can locate the FASTQ files.
sorter2_stage1a(
  reads_dir  = "raw_reads/",
  output_dir = "results/stage1a/",
  trim       = TRUE,
  spades     = TRUE,
  conda_env  = "SORTER2"
)

# Stage 1B — cluster assembled contigs into ortholog groups
sorter2_stage1b(
  input_dir  = "results/stage1a/",
  output_dir = "results/stage1b/",
  ref        = "reference/baits.fasta",
  loci       = 100,
  conda_env  = "SORTER2"
)

# Stage 2 — phase bi-allelic haplotypes
# input_assemblies: Stage 1A output (*_assembly/ dirs)
# input_clusters:   Stage 1B output (diploidclusters/ and diploids/)
# reads_dir:        directory with raw FASTQ files (*_R1.fastq / *_R2.fastq).
#   Required when Stage 1A was run with trim = FALSE (pre-cleaned reads).
#   Omit when trim = TRUE — reads are already inside the *_assembly/ dirs.
sorter2_stage2(
  input_assemblies = "results/stage1a/",
  input_clusters   = "results/stage1b/",
  output_dir       = "results/stage2/",
  reads_dir        = "raw_reads/",
  conda_env        = "SORTER2"
)

# Stage 3 — phase hybrid loci
# phaseset_dir: where to find hybrid *_assembly/ dirs (usually Stage 1A output).
#   The wrapper copies them into output_dir/phaseset/ automatically.
# hybrid_samples: names of hybrid samples to copy; omit to copy all *_assembly/.
sorter2_stage3(
  input_phased     = "results/stage2/",
  input_assemblies = "results/stage1a/",
  output_dir       = "results/stage3/",
  phaseset_dir     = "results/stage1a/",
  hybrid_samples   = c("Sample1_species", "Sample2_species"),
  ref              = "reference/baits.fasta",
  loci             = 100,
  conda_env        = "SORTER2"
)

# Processor — filter and summarise final outputs
sorter2_processor(
  input_dir  = "results/stage3/",
  output_dir = "results/processed/",
  conda_env  = "SORTER2"
)
```

Set `dry_run = TRUE` on any call to preview the command without running it.
Set `SORTER2R_CONDA_ENV` once to avoid repeating `conda_env` on every call:

```r
Sys.setenv(SORTER2R_CONDA_ENV = "SORTER2")
```

## Development usage from source

1. Vendor scripts to `inst/python/` (see sync script below).
2. Pass `script_dir = "inst/python"` to any wrapper (or set no `script_dir`
   once the package is installed — it resolves automatically).
3. Use a configured Python/conda environment that contains SORTER2 tools.

## Upstream sync workflow

Run from package root:

```bash
bash tools/sync_upstream.sh
```

This script copies core Python scripts plus `SORTER2.yml` from upstream into:

- `inst/python/` (runtime scripts)
- `inst/upstream/` (tracking copy)

Override upstream branch if needed:

```bash
UPSTREAM_REF=main bash tools/sync_upstream.sh
```
