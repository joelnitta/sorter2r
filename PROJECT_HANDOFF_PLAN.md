# sorter2r project handoff plan

## Why this exists

This file captures current context so work can resume quickly after switching to
another repository.

## Current objective

Create and maintain an R package that wraps core SORTER2 Python scripts while
staying in sync with upstream changes from JonasMendez/SORTER2.

## Decisions already made

- Package location model: separate repository for sorter2r.
- V1 implementation strategy: thin wrappers that call Python CLI scripts.
- Upstream sync strategy: vendor core scripts and refresh with a sync script.
- V1 scope: core pipeline only.
  FormatReads, Stage1A, Stage1B, Stage2, Stage3, Processor.

## Current implementation status in this working tree

- Package scaffold exists in this directory.
- Core wrapper functions exist in R/wrappers-core.R.
- Shared execution helpers exist in R/aaa-utils.R.
- Vendored scripts exist in inst/python and tracking copies in inst/upstream.
- Upstream sync script exists at tools/sync_upstream.sh.

## Immediate next steps after switching repos

1. Create the separate sorter2r repository and copy this directory into it.
2. Run basic checks: source R files and run dry-run wrapper calls.
3. Add testthat tests for argument validation and command serialization.
4. Add a small smoke test path for Stage1A -> Stage1B -> Stage2 -> Processor.
5. Add CI for R CMD check and wrapper tests.
6. Add scheduled upstream drift check and document release checklist.

## Sync workflow reference

From package root:

```bash
bash tools/sync_upstream.sh
```

Optional override for branch/ref:

```bash
UPSTREAM_REF=main bash tools/sync_upstream.sh
```

## Compatibility notes

- Wrappers assume external bioinformatics tools are available in the selected
  Python/conda environment.
- Boolean flags are normalized to T/F to match SORTER2 CLI behavior.
- Paths for -wd arguments should include a trailing slash; wrappers normalize
  this when possible.

## Suggested first commit in new repo

- Add scaffold and wrappers.
- Add vendored core scripts.
- Add sync script.
- Add this handoff file.
