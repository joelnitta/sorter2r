# Reformat a reference FASTA for SORTER2

Rewrites sequence IDs to the `>L{n}_{Family-Genus-species}` format
required by Stage 1B. The locus number `n` is taken from the `L{n}`
substring after the last `-` in the original ID; sequences whose
post-hyphen tag does not already carry an `L{n}_` prefix are assigned
sequential numbers starting after the highest embedded one. The
taxonomic label (fields 2–4 of the original `_`-delimited name) is
appended with `-` as a separator, ensuring every output ID is unique.

## Usage

``` r
sorter2_format_ref(input, output = NULL)
```

## Arguments

- input:

  Path to the input reference FASTA file.

- output:

  Path for the reformatted output file. Defaults to the input path with
  `_sorter2` appended before the extension.

## Value

Path to the output file.
