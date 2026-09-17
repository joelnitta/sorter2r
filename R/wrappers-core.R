#' Run the SORTER2 FormatReads step
#'
#' @param input Input CSV/TSV file mapping original filenames to sample IDs.
#' @param reads_dir Directory containing the FASTQ files to rename.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @param verbose Integer verbosity level. `1L` (default) shows Python-level
#'   section headers; `0L` is fully silent; `2L` echoes the command, enables
#'   per-item debug output (`-v`), and shows subprocess tool output.
#' @return Path to `reads_dir` (for use with `tar_file()`).
#' @export
sorter2_format_reads <- function(
  input,
  reads_dir = getwd(),
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  dry_run = FALSE,
  verbose = 1L
) {
  input <- sorter2_require_file(input, "input")
  reads_dir <- sorter2_require_dir(reads_dir, "reads_dir")

  args <- c("-i", basename(input))
  if (verbose >= 2L) args <- c(args, "-v")
  res <- sorter2_run(
    script = "SORTER2_FormatReads.py",
    args = args,
    python = python,
    conda_env = conda_env,
    conda = conda,
    script_dir = script_dir,
    working_dir = reads_dir,
    dry_run = dry_run,
    echo = (verbose >= 2L),
    show_output = (verbose >= 1L)
  )
  if (dry_run) return(invisible(res))
  reads_dir
}

#' Run the SORTER2 Stage1A step
#'
#' @param reads_dir Directory containing the input FASTQ files.
#' @param output_dir Directory where the assembled output will be created.
#' @param spades Logical flag for SPAdes assembly.
#' @param trim Logical flag for TrimGalore trimming.
#' @param clean_spades_tmp If `TRUE`, delete SPAdes intermediate K-mer graph
#'   directories (`K21/`, `K33/`, etc.) after each sample assembles. Saves
#'   substantial disk space at the cost of preventing SPAdes restarts.
#'   Default `FALSE`.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param overwrite If `TRUE`, delete `output_dir` before running.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @param verbose Integer verbosity level. `1L` (default) shows Python-level
#'   section headers; `0L` is fully silent; `2L` echoes the command, enables
#'   per-item debug output (`-v`), and shows subprocess tool output.
#' @return Path to `output_dir` (for use with `tar_file()`).
#' @export
sorter2_stage1a <- function(
  reads_dir,
  output_dir,
  spades = TRUE,
  trim = TRUE,
  clean_spades_tmp = FALSE,
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  overwrite = FALSE,
  dry_run = FALSE,
  verbose = 1L
) {
  reads_dir <- sorter2_require_dir(reads_dir, "reads_dir")
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  # normalizePath may return a relative path for non-existent dirs on some
  # Linux platforms; the Python script runs with working_dir = reads_dir, so
  # a relative output_dir would be resolved inside the reads directory.
  if (!startsWith(output_dir, "/")) {
    output_dir <- file.path(getwd(), output_dir)
  }
  if (!dry_run) sorter2_check_overwrite(output_dir, overwrite)

  args <- c(
    "-o",
    output_dir,
    "-spades",
    sorter2_bool_flag(spades),
    "-trim",
    sorter2_bool_flag(trim),
    "-clean_tmp",
    sorter2_bool_flag(clean_spades_tmp)
  )

  res <- sorter2_run(
    script = "SORTER2_Stage1A_TrimSPAdes.py",
    args = args,
    python = python,
    conda_env = conda_env,
    conda = conda,
    script_dir = script_dir,
    working_dir = reads_dir,
    dry_run = dry_run,
    echo = (verbose >= 2L),
    show_output = (verbose >= 1L)
  )
  if (dry_run) return(invisible(res))
  output_dir
}

#' Run the SORTER2 Stage1B step
#'
#' @param input_dir Directory containing the Stage1A output
#'   (`*_assembly/` subdirs with SPAdes contigs).
#' @param output_dir Directory where `diploids/` and `diploidclusters/` will
#'   be created.
#' @param ref Reference file path.
#' @param loci Number of loci to process.
#' @param clust2id Clustering threshold.
#' @param recluster Logical flag for reclustering. When `TRUE`, reuses
#'   `assembly_workfiles/` (assembled contigs) from a previous Stage1B run
#'   in `output_dir` instead of re-assembling.
#' @param contignum Contig count threshold.
#' @param contiglen Contig length threshold.
#' @param aliter Alignment iterations.
#' @param indelrep Indel representation threshold.
#' @param idformat Identifier format.
#' @param verbose Integer verbosity level. `1L` (default) shows Python-level
#'   section headers; `0L` is fully silent; `2L` echoes the command, enables
#'   per-item debug output (`-v`), and shows subprocess tool output.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param overwrite If `TRUE`, remove existing output before running. When
#'   `recluster = TRUE`, only prior Stage1B output is removed and
#'   `assembly_workfiles/` is preserved so `-reclust` can reuse it.
#'   Otherwise `output_dir` is deleted entirely.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @return Path to `output_dir` (for use with `tar_file()`).
#' @export
sorter2_stage1b <- function(
  input_dir,
  output_dir,
  ref,
  loci,
  clust2id = 0.70,
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
) {
  input_dir <- sorter2_with_trailing_slash(input_dir)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (!startsWith(output_dir, "/")) {
    output_dir <- file.path(getwd(), output_dir)
  }
  ref <- sorter2_require_file(ref, "ref")

  workfiles_path <- file.path(output_dir, "assembly_workfiles")

  if (!dry_run) {
    if (dir.exists(output_dir) && overwrite) {
      if (recluster) {
        # -reclust depends on assembly_workfiles/ from the previous run;
        # only clear prior Stage1B output, not the workfiles it needs.
        entries <- list.files(output_dir, full.names = TRUE)
        for (e in setdiff(entries, workfiles_path)) unlink(e, recursive = TRUE)
      } else {
        message("Removing existing output: ", output_dir)
        unlink(output_dir, recursive = TRUE)
      }
    } else if (dir.exists(output_dir) && !overwrite) {
      stop(
        sprintf(
          "Output already exists: %s\nSet overwrite = TRUE to remove it and re-run.",
          output_dir
        ),
        call. = FALSE
      )
    }
  }

  args <- c(
    "-wd",
    input_dir,
    "-outdir",
    output_dir,
    "-ref",
    ref,
    "-loci",
    as.character(as.integer(loci)),
    "-cl",
    as.character(clust2id),
    "-reclust",
    sorter2_bool_flag(recluster),
    "-csn",
    as.character(as.integer(contignum)),
    "-csl",
    as.character(as.integer(contiglen)),
    "-al",
    as.character(as.integer(aliter)),
    "-indel",
    as.character(indelrep),
    "-idformat",
    as.character(idformat)
  )
  if (verbose >= 2L) args <- c(args, "-v")

  res <- sorter2_run(
    script = "SORTER2_Stage1B_AssembleOrthologs.py",
    args = args,
    python = python,
    conda_env = conda_env,
    conda = conda,
    script_dir = script_dir,
    dry_run = dry_run,
    echo = (verbose >= 2L),
    show_output = (verbose >= 1L)
  )
  if (dry_run) return(invisible(res))
  output_dir
}

#' Run the SORTER2 Stage2 step
#'
#' @param input_assemblies Directory containing the Stage1A output
#'   (`*_assembly/` subdirs with SPAdes contigs).
#' @param input_clusters Directory containing the Stage1B output
#'   (`diploidclusters/` and `diploids/` subdirs).
#' @param output_dir Directory where `diploids_phased/` will be created.
#' @param reads_dir Optional directory containing raw FASTQ files
#'   (`*_R1.fastq` / `*_R2.fastq`). Required when Stage 1A was run with
#'   `trim = FALSE`. When `NULL` (the default), Stage2 expects Trim Galore
#'   output (`*_R1_val_1.fq`) inside each `*_assembly/` subdirectory.
#' @param phasequal Phase quality threshold.
#' @param aliter Alignment iterations.
#' @param indelrep Indel representation threshold.
#' @param idformat Identifier format.
#' @param threads Number of parallel worker processes for per-sample phasing
#'   and per-locus alignment. Default `1L` (serial).
#' @param verbose Integer verbosity level. `1L` (default) shows Python-level
#'   section headers; `0L` is fully silent; `2L` echoes the command, enables
#'   per-item debug output (`-v`), and shows subprocess tool output.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param overwrite If `TRUE`, delete `output_dir` before running.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @return Path to `output_dir` (for use with `tar_file()`).
#' @export
sorter2_stage2 <- function(
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
) {
  input_assemblies <- sorter2_with_trailing_slash(input_assemblies)
  input_clusters <- sorter2_with_trailing_slash(input_clusters)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (!startsWith(output_dir, "/")) {
    output_dir <- file.path(getwd(), output_dir)
  }
  if (!dry_run) sorter2_check_overwrite(output_dir, overwrite)

  args <- c(
    "-wa",
    input_assemblies,
    "-wc",
    input_clusters,
    "-outdir",
    output_dir,
    "-pq",
    as.character(as.integer(phasequal)),
    "-al",
    as.character(as.integer(aliter)),
    "-indel",
    as.character(indelrep),
    "-idformat",
    as.character(idformat),
    "-t",
    as.character(as.integer(threads))
  )

  if (!is.null(reads_dir)) {
    reads_dir <- sorter2_with_trailing_slash(reads_dir)
    args <- c(args, "-reads", reads_dir)
  }
  if (verbose >= 2L) args <- c(args, "-v")

  res <- sorter2_run(
    script = "SORTER2_Stage2_PhaseOrthologs.py",
    args = args,
    python = python,
    conda_env = conda_env,
    conda = conda,
    script_dir = script_dir,
    dry_run = dry_run,
    echo = (verbose >= 2L),
    show_output = (verbose >= 1L)
  )
  if (dry_run) return(invisible(res))
  output_dir
}

#' Run the SORTER2 Stage3 step
#'
#' @param input_phased Directory containing the Stage2 output. Must contain a
#'   `diploids_phased/` subdirectory with the phased diploid sequences.
#' @param input_assemblies Directory containing the Stage1A output. Used to
#'   identify diploid sample names from `*_assembly/` subdirectories.
#' @param output_dir Directory where Stage3 results will be written.
#' @param phaseset_dir Optional directory containing hybrid sample
#'   `*_assembly/` subdirectories (SPAdes output). When provided, the wrapper
#'   copies those dirs into `output_dir/phaseset/` before running. When `NULL`
#'   (the default), `output_dir/phaseset/` must already exist and be populated
#'   by the user before calling this function.
#' @param hybrid_samples Optional character vector of sample base names (e.g.
#'   `"Iimura46_cgrande"`) used to filter which `*_assembly/` dirs are copied
#'   from `phaseset_dir`. When `NULL`, all `*_assembly/` dirs are copied.
#' @param reads_dir Optional directory containing raw FASTQ files
#'   (`*_R1.fastq` / `*_R2.fastq`). Required when Stage 1A was run with
#'   `trim = FALSE`. When `NULL` (the default), Stage3 expects Trim Galore
#'   output (`*_R1_val_1.fq`) inside each `phaseset/*_assembly/` subdirectory.
#' @param ref Reference file path.
#' @param loci Number of loci to process.
#' @param contigscafnum Contig/scaffold count threshold.
#' @param contigscaflen Contig/scaffold length threshold.
#' @param phasequal Phase quality threshold.
#' @param aliter Alignment iterations.
#' @param indelrep Indel representation threshold.
#' @param filterundiff Logical flag for filtering undifferentiated loci.
#' @param tag_specimen Logical flag. When `TRUE`, a hybrid's alt-lineage
#'   haplotype header also records the specific diploid specimen (voucher)
#'   it matched, not just the species code, so a progenitor mapfile can
#'   target a specific voucher. Default `FALSE` preserves the
#'   original species-only header.
#' @param verbose Integer verbosity level. `1L` (default) shows Python-level
#'   section headers; `0L` is fully silent; `2L` echoes the command, enables
#'   per-item debug output (`-v`), and shows subprocess tool output.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param overwrite If `TRUE`, remove Stage3 outputs (but preserve
#'   `phaseset/`) before re-running. When `phaseset_dir` is also set,
#'   `phaseset/` is fully refreshed from the source.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @return Path to `output_dir` (for use with `tar_file()`).
#' @export
sorter2_stage3 <- function(
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
  tag_specimen = FALSE,
  verbose = 1L,
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  overwrite = FALSE,
  dry_run = FALSE
) {
  input_phased <- sorter2_with_trailing_slash(input_phased)
  input_assemblies <- sorter2_with_trailing_slash(input_assemblies)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (!startsWith(output_dir, "/")) {
    output_dir <- file.path(getwd(), output_dir)
  }
  ref <- sorter2_require_file(ref, "ref")

  phaseset_path <- file.path(output_dir, "phaseset")

  if (!dry_run) {
    if (dir.exists(output_dir) && overwrite) {
      # Remove Stage3 outputs but keep phaseset/ unless phaseset_dir refreshes it
      entries <- list.files(output_dir, full.names = TRUE)
      keep <- if (is.null(phaseset_dir)) phaseset_path else character(0)
      for (e in setdiff(entries, keep)) unlink(e, recursive = TRUE)
    } else if (dir.exists(output_dir) && !overwrite) {
      stop(
        sprintf(
          "Output already exists: %s\nSet overwrite = TRUE to re-run.",
          output_dir
        ),
        call. = FALSE
      )
    }

    if (!is.null(phaseset_dir)) {
      phaseset_dir <- sorter2_with_trailing_slash(phaseset_dir)
      dir.create(phaseset_path, recursive = TRUE, showWarnings = FALSE)
      asm_dirs <- list.dirs(phaseset_dir, recursive = FALSE, full.names = TRUE)
      asm_dirs <- asm_dirs[grepl("_assembly$", basename(asm_dirs))]
      if (!is.null(hybrid_samples)) {
        asm_dirs <- asm_dirs[
          sub("_assembly$", "", basename(asm_dirs)) %in% hybrid_samples
        ]
      }
      if (length(asm_dirs) == 0L) {
        if (!is.null(hybrid_samples) && length(hybrid_samples) == 0L) {
          warning(
            "hybrid_samples is empty; no hybrid assemblies to process. ",
            "Stage 3 skipped.",
            call. = FALSE
          )
          return(invisible(output_dir))
        }
        stop(
          "No *_assembly/ directories found in phaseset_dir: ", phaseset_dir,
          call. = FALSE
        )
      }
      for (src in asm_dirs) {
        dst <- file.path(phaseset_path, basename(src))
        if (dir.exists(dst)) unlink(dst, recursive = TRUE)
        file.copy(src, phaseset_path, recursive = TRUE)
      }
      message(
        "Populated phaseset/ with ", length(asm_dirs),
        " assembly dir(s) from ", phaseset_dir
      )
    }

    if (!dir.exists(phaseset_path)) {
      stop(
        "output_dir/phaseset/ does not exist: ", phaseset_path,
        "\nEither set phaseset_dir to auto-populate it, or create it manually ",
        "with hybrid sample *_assembly/ subdirectories.",
        call. = FALSE
      )
    }
  }

  args <- c(
    "-wp",
    input_phased,
    "-wa",
    input_assemblies,
    "-outdir",
    output_dir,
    "-ref",
    ref,
    "-loci",
    as.character(as.integer(loci)),
    "-csn",
    as.character(as.integer(contigscafnum)),
    "-csl",
    as.character(as.integer(contigscaflen)),
    "-pq",
    as.character(as.integer(phasequal)),
    "-al",
    as.character(as.integer(aliter)),
    "-indel",
    as.character(indelrep),
    "-fp",
    sorter2_bool_flag(filterundiff),
    "-specimen",
    sorter2_bool_flag(tag_specimen)
  )
  if (!is.null(reads_dir)) {
    reads_dir <- sorter2_with_trailing_slash(reads_dir)
    args <- c(args, "-reads", reads_dir)
  }
  if (verbose >= 2L) args <- c(args, "-v")

  res <- sorter2_run(
    script = "SORTER2_Stage3_PhaseHybrids.py",
    args = args,
    python = python,
    conda_env = conda_env,
    conda = conda,
    script_dir = script_dir,
    dry_run = dry_run,
    echo = (verbose >= 2L),
    show_output = (verbose >= 1L)
  )
  if (dry_run) return(invisible(res))
  output_dir
}

#' Run the SORTER2 Hapl-O-Miner step
#'
#' Mines haploid organellar genomes (e.g. chloroplast) from off-target
#' by-catch reads in a SORTER2 Stage 1A dataset.  Reads from each
#' `*_assembly/` subdirectory are mapped to a user-supplied haploid
#' organellar reference; per-sample consensus sequences are generated
#' and filtered by minimum coverage and read-depth thresholds.
#'
#' @param input_dir Directory containing Stage 1A `*_assembly/`
#'   subdirectories.
#' @param output_dir Directory where all Hapl-O-Miner output will be
#'   written.
#' @param organellar_ref Path to the haploid organellar reference FASTA
#'   (e.g. a chloroplast genome).
#' @param coverage Minimum reference coverage percentage threshold.
#'   Samples below this value are excluded from the final alignment.
#' @param depth Minimum mean read-depth threshold.  Samples below this
#'   value are excluded from the final alignment.
#' @param reads_dir Optional directory containing raw FASTQ files
#'   (`*_R1.fastq` / `*_R2.fastq`). Required when Stage 1A was run with
#'   `trim = FALSE`. When `NULL` (the default), Hapl-O-Miner expects
#'   Trim Galore output (`*_R1_val_1.fq`) inside each `*_assembly/`
#'   subdirectory.
#' @param threads Number of parallel worker processes for per-sample
#'   mapping. Default `1L` (serial).
#' @param bwa_threads Threads passed to each `bwa mem` call (`-t`). Runs
#'   inside each of the `threads` worker processes, so total CPU usage is
#'   `threads * bwa_threads` - keep the product within your core budget.
#'   Default `1L` (matches the previous single-threaded behavior).
#' @param clean_workfiles If `TRUE`, delete the intermediate workfiles
#'   directory (BAMs, raw consensus FASTAs, per-sample stat files) after
#'   the run completes. Saves substantial disk space; the compiled
#'   `readstats_cp.csv` and final consensus FASTAs are unaffected.
#'   Default `FALSE`.
#' @param align If `TRUE` (default), align the filtered whole-consensus
#'   sequences with mafft (fast single-pass mode) into
#'   `HaplOMiner_*_filtered_al.fasta`. This is separate from, and not
#'   required for, the per-sample consensus FASTAs written to
#'   `all_chloroplasts/`. Callers that re-split each consensus per-locus
#'   themselves (e.g. with BLAT) and never read this alignment should
#'   pass `FALSE` - even fast mode does not scale well to large sample
#'   groups (all-pairs distance calc on whole-organellar-length
#'   sequences with FFT disabled), and can run for hours once a group
#'   reaches dozens of samples.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @param verbose Integer verbosity level. `1L` (default) shows
#'   Python-level section headers; `0L` is fully silent; `2L` echoes
#'   the command and enables per-item debug output (`-v`).
#' @return Path to `output_dir` (for use with `tar_file()`).
#' @export
sorter2_haplominer <- function(
  input_dir,
  output_dir,
  organellar_ref,
  coverage = 50,
  depth = 5,
  reads_dir = NULL,
  threads = 1L,
  bwa_threads = 1L,
  clean_workfiles = FALSE,
  align = TRUE,
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  dry_run = FALSE,
  verbose = 1L
) {
  input_dir <- sorter2_with_trailing_slash(input_dir)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (!startsWith(output_dir, "/")) {
    output_dir <- file.path(getwd(), output_dir)
  }
  organellar_ref <- sorter2_require_file(organellar_ref, "organellar_ref")

  args <- c(
    "-wd",
    input_dir,
    "-outdir",
    output_dir,
    "-cpref",
    organellar_ref,
    "-c",
    as.character(as.numeric(coverage)),
    "-d",
    as.character(as.numeric(depth)),
    "-t",
    as.character(as.integer(threads)),
    "-bwa_t",
    as.character(as.integer(bwa_threads)),
    "-clean_workfiles",
    sorter2_bool_flag(clean_workfiles),
    "-align",
    sorter2_bool_flag(align)
  )

  if (!is.null(reads_dir)) {
    reads_dir <- sorter2_with_trailing_slash(reads_dir)
    args <- c(args, "-reads", reads_dir)
  }
  if (verbose >= 2L) args <- c(args, "-v")

  res <- sorter2_run(
    script = "SORTER2_HaplOMiner.py",
    args = args,
    python = python,
    conda_env = conda_env,
    conda = conda,
    script_dir = script_dir,
    dry_run = dry_run,
    echo = (verbose >= 2L),
    show_output = (verbose >= 1L)
  )
  if (dry_run) return(invisible(res))
  output_dir
}

#' Run the SORTER2 Processor step
#'
#' @param input_dir Directory containing the SORTER2 pipeline outputs
#'   (`diploids/`, `diploidclusters/`, and optionally `diploids_phased/` and
#'   `phaseset/`).
#' @param output_dir Directory where processed results will be written.
#' @param repfilt Repeat filter threshold.
#' @param majorclusters Number of major clusters.
#' @param keepal Logical flag for keeping alignments.
#' @param stage2 Logical flag for including Stage2 output.
#' @param stage3 Logical flag for including Stage3 output.
#' @param dovcf Logical flag for writing VCF output.
#' @param reads_dir Optional directory containing raw FASTQ files
#'   (`*_R1.fastq` / `*_R2.fastq`). Required when Stage 1A was run
#'   with `trim = FALSE`. When `NULL` (the default), the Processor
#'   expects Trim Galore output (`*_R1_val_1.fq`) inside each
#'   `*_assembly/` subdirectory.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param success_codes Integer vector of exit codes treated as success.
#'   Default `0L`. Pass `c(0L, 1L)` if outputs are known to be correct
#'   and the script is exiting 1 on a re-run (e.g. interrupted pipeline).
#' @param overwrite If `TRUE`, delete `output_dir` before running.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @param verbose Integer verbosity level. `1L` (default) shows Python-level
#'   section headers; `0L` is fully silent; `2L` echoes the command, enables
#'   per-item debug output (`-v`), and shows subprocess tool output.
#' @return Path to `output_dir` (for use with `tar_file()`).
#' @export
sorter2_processor <- function(
  input_dir,
  output_dir = input_dir,
  repfilt = 0.50,
  majorclusters = 2,
  keepal = TRUE,
  stage2 = FALSE,
  stage3 = FALSE,
  dovcf = TRUE,
  reads_dir = NULL,
  success_codes = 0L,
  overwrite = FALSE,
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  dry_run = FALSE,
  verbose = 1L
) {
  input_dir <- sorter2_with_trailing_slash(input_dir)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (!dry_run) sorter2_check_overwrite(output_dir, overwrite)

  args <- c(
    "-wd",
    input_dir,
    "-outdir",
    output_dir,
    "-rep",
    as.character(repfilt),
    "-mc",
    as.character(as.integer(majorclusters)),
    "-al",
    sorter2_bool_flag(keepal),
    "-st2",
    sorter2_bool_flag(stage2),
    "-st3",
    sorter2_bool_flag(stage3),
    "-dovcf",
    sorter2_bool_flag(dovcf)
  )

  if (!is.null(reads_dir)) {
    reads_dir <- sorter2_with_trailing_slash(reads_dir)
    args <- c(args, "-reads", reads_dir)
  }
  if (verbose >= 2L) args <- c(args, "-v")

  res <- sorter2_run(
    script = "SORTER2_Processor.py",
    args = args,
    python = python,
    conda_env = conda_env,
    conda = conda,
    script_dir = script_dir,
    dry_run = dry_run,
    echo = (verbose >= 2L),
    show_output = (verbose >= 1L),
    success_codes = success_codes
  )
  if (dry_run) return(invisible(res))
  output_dir
}

#' Process progenitor sequences from Stage 3 output
#'
#' Runs \code{SORTER2_ProgenitorProcessor.py}. Collapses progenitor IDs
#' to user-defined clades, removes outgroup sequences, filters by minimum
#' sequence count, removes duplicate IDs, and optionally retains only loci
#' with at least two differentiated progenitor clades.
#'
#' The Stage 3 \code{input_dir} is never modified, so its targets hash
#' remains stable.
#'
#' @param input_dir Directory containing Stage 3 FASTA output (read-only).
#' @param output_dir Directory to write processed output into. Created if
#'   it does not exist.
#' @param mapfile Path to a CSV file with three columns:
#'   \code{hybrid}, \code{progenitor}, \code{clade}. \code{hybrid} may
#'   be either a bare species code, applying to every voucher with
#'   that code, or a \code{voucher_species} compound key (e.g.
#'   \code{"Nitta4471_ointermedia"}), applying to that one voucher
#'   only. A compound row takes precedence over a species-level row
#'   for the same \code{progenitor} value. Use compound rows when a
#'   species is not monophyletic and different vouchers need
#'   different progenitor clades; a species-level match prints a
#'   monophyly-reminder warning.
#' @param outgroups Path to a plain-text file listing outgroup IDs to
#'   remove (one per line).
#' @param minseq Minimum number of phased sequence pairs required per
#'   progenitor set across all loci. The Python script doubles this
#'   value internally before applying the filter. Default \code{1L}.
#' @param filterundiff If \code{TRUE}, discard loci where a hybrid
#'   sample has fewer than two distinct progenitor clades represented.
#'   Default \code{FALSE}.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name.
#' @param conda Conda executable to use when \code{conda_env} is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param dry_run If \code{TRUE}, return the command without executing it.
#' @param verbose Integer verbosity level. \code{1L} (default) shows
#'   Python-level section headers; \code{0L} is fully silent; \code{2L}
#'   echoes the command and shows subprocess tool output.
#' @return Path to \code{output_dir} (for use with \code{tar_file()}).
#' @export
sorter2_progenitor_processor <- function(
  input_dir,
  output_dir,
  mapfile,
  outgroups,
  minseq = 1L,
  filterundiff = FALSE,
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  dry_run = FALSE,
  verbose = 1L
) {
  input_dir  <- sorter2_with_trailing_slash(input_dir)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  mapfile    <- sorter2_require_file(mapfile, "mapfile")
  outgroups  <- sorter2_require_file(outgroups, "outgroups")

  if (!dry_run) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  args <- c(
    "-wd",    paste0(output_dir, "/"),
    "-indir", input_dir,
    "-map",   mapfile,
    "-og",    outgroups,
    "-min",   as.character(as.integer(minseq)),
    "-dif",   sorter2_bool_flag(filterundiff)
  )
  if (verbose >= 2L) args <- c(args, "-v")

  res <- sorter2_run(
    script      = "SORTER2_ProgenitorProcessor.py",
    args        = args,
    python      = python,
    conda_env   = conda_env,
    conda       = conda,
    script_dir  = script_dir,
    dry_run     = dry_run,
    echo        = (verbose >= 2L),
    show_output = (verbose >= 1L)
  )
  if (dry_run) return(invisible(res))
  output_dir
}
