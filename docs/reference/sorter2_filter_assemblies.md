# Copy ingroup assembly folders from Stage 1A output, excluding outgroups

SORTER2 Stages 1B onward should be run on ingroup samples only; outgroup
assemblies inflate allele counts and ADMIXTURE complexity. This function
creates a new directory containing only the ingroup `*_assembly` folders
from a Stage 1A output directory, leaving the original untouched.

## Usage

``` r
sorter2_filter_assemblies(
  input_dir,
  output_dir,
  exclude_samples = character(),
  overwrite = FALSE
)
```

## Arguments

- input_dir:

  Directory containing Stage 1A output (the `*_assembly` subdirectories
  produced by [`sorter2_stage1a()`](sorter2_stage1a.md)).

- output_dir:

  Directory to create containing only ingroup assemblies. Must not
  already exist unless `overwrite = TRUE`.

- exclude_samples:

  Character vector of sample names to exclude (outgroups). Names should
  match the `{sample}_assembly` folder prefix — for example,
  `"Wade3565_cintermedium"` excludes `Wade3565_cintermedium_assembly/`.

- overwrite:

  If `TRUE`, delete `output_dir` before running.

## Value

Path to `output_dir` (for use with `tar_file()`).
