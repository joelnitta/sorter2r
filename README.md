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

## Development usage from source

1. Vendor scripts to `inst/python/` (see sync script below).
2. Point wrappers to that folder with `script_dir`.
3. Use a configured Python/conda environment that contains SORTER2 tools.

Example:

```r
res <- sorter2_stage1a(
  projname = "ProjectName",
  trim = TRUE,
  spades = TRUE,
  script_dir = "inst/python",
  conda_env = "SORTER2",
  working_dir = ".",
  dry_run = TRUE
)

res$command
```

You can also set `SORTER2R_CONDA_ENV` once instead of passing
`conda_env` to each wrapper. If needed, set `SORTER2R_CONDA` to the path
of a specific conda executable.

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
