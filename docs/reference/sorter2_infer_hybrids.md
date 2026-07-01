# Identify candidate hybrid or polyploid samples from Stage 1B output

Reads the consensus-allele summary table produced by
[`sorter2_stage1b()`](sorter2_stage1b.md) and computes per-sample
allele-count statistics useful for distinguishing diploids from hybrids
or polyploids. Polyploids carry extra genome copies, so their total
consensus-allele count across loci is proportionally elevated relative
to diploids.

## Usage

``` r
sorter2_infer_hybrids(input_dir, min_mean = 2)
```

## Arguments

- input_dir:

  Directory containing Stage 1B output. Must contain a `diploids/`
  subdirectory with a file matching
  `ALLsamples_consensusallele_*_SUMMARY_TABLE.csv`.

- min_mean:

  Minimum mean consensus alleles per locus required to flag a sample as
  a hybrid/polyploid candidate. Default `2.0`: samples with on average
  more than two allele copies per locus are flagged.

## Value

A data frame with one row per sample, sorted by descending
`mean_per_locus`, with columns:

- sample:

  Sample name (trailing `_` stripped; matches the names accepted by the
  `hybrid_samples` argument of [`sorter2_stage3()`](sorter2_stage3.md)).

- total_alleles:

  Sum of consensus alleles across all loci.

- n_loci:

  Number of loci with at least one consensus allele.

- mean_per_locus:

  `total_alleles / n_loci`. Primary ploidy proxy: diploids typically
  show 1–2; values above 2 suggest polyploidy.

- candidate:

  `TRUE` if `mean_per_locus >= min_mean`.

## Details

Interpretation guidance (from SORTER2 documentation, full-dataset
scale):

- Diploids: ~1–2 consensus alleles per locus on average

- Tetraploids: ~2–3 per locus

- Hexaploids: ~4+ per locus

Note that absolute counts depend on the number of loci analysed; use
`mean_per_locus` for cross-dataset comparisons rather than
`total_alleles`.

For a complete hybrid assessment, combine these allele-count statistics
with ADMIXTURE ancestry proportions and phylogenetic placement.
