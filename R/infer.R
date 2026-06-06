#' Identify candidate hybrid or polyploid samples from Stage 1B output
#'
#' Reads the consensus-allele summary table produced by [sorter2_stage1b()]
#' and computes per-sample allele-count statistics useful for distinguishing
#' diploids from hybrids or polyploids. Polyploids carry extra genome copies,
#' so their total consensus-allele count across loci is proportionally
#' elevated relative to diploids.
#'
#' Interpretation guidance (from SORTER2 documentation, full-dataset scale):
#' - Diploids: ~1–2 consensus alleles per locus on average
#' - Tetraploids: ~2–3 per locus
#' - Hexaploids: ~4+ per locus
#'
#' Note that absolute counts depend on the number of loci analysed; use
#' `mean_per_locus` for cross-dataset comparisons rather than `total_alleles`.
#'
#' For a complete hybrid assessment, combine these allele-count statistics
#' with ADMIXTURE ancestry proportions and phylogenetic placement.
#'
#' @param input_dir Directory containing Stage 1B output. Must contain a
#'   `diploids/` subdirectory with a file matching
#'   `ALLsamples_consensusallele_*_SUMMARY_TABLE.csv`.
#' @param min_mean Minimum mean consensus alleles per locus required to flag
#'   a sample as a hybrid/polyploid candidate. Default `2.0`: samples with on
#'   average more than two allele copies per locus are flagged.
#' @return A data frame with one row per sample, sorted by descending
#'   `mean_per_locus`, with columns:
#'   \describe{
#'     \item{sample}{Sample name (trailing `_` stripped; matches the names
#'       accepted by the `hybrid_samples` argument of [sorter2_stage3()]).}
#'     \item{total_alleles}{Sum of consensus alleles across all loci.}
#'     \item{n_loci}{Number of loci with at least one consensus allele.}
#'     \item{mean_per_locus}{`total_alleles / n_loci`. Primary ploidy proxy:
#'       diploids typically show 1–2; values above 2 suggest polyploidy.}
#'     \item{candidate}{`TRUE` if `mean_per_locus >= min_mean`.}
#'   }
#' @importFrom utils read.csv
#' @export
sorter2_infer_hybrids <- function(input_dir, min_mean = 2.0) {
  diploids_dir <- file.path(input_dir, "diploids")
  csvs <- list.files(
    diploids_dir,
    pattern = "ALLsamples_consensusallele_.*_SUMMARY_TABLE\\.csv$",
    full.names = TRUE
  )
  if (length(csvs) == 0L) {
    stop(
      "No consensus-allele summary table found in: ", diploids_dir, "\n",
      "Expected a file matching ",
      "ALLsamples_consensusallele_*_SUMMARY_TABLE.csv",
      call. = FALSE
    )
  }
  tbl <- read.csv(csvs[[1L]], check.names = FALSE)

  sample_data <- tbl[, -1L, drop = FALSE]

  rows <- lapply(names(sample_data), function(s) {
    counts       <- sample_data[[s]]
    total        <- sum(counts)
    n_loci       <- sum(counts > 0L)
    mean_per_loc <- if (n_loci > 0L) total / n_loci else 0
    data.frame(
      sample         = sub("_$", "", s),
      total_alleles  = total,
      n_loci         = n_loci,
      mean_per_locus = round(mean_per_loc, 2),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out$candidate <- out$mean_per_locus >= min_mean

  if (max(out$n_loci) < 50L) {
    warning(
      "Fewer than 50 loci per sample; mean_per_locus estimates may be ",
      "unstable. Results are indicative only.",
      call. = FALSE
    )
  }

  out[order(-out$mean_per_locus), ]
}
