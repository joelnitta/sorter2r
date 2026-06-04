sorter2_bool_flag <- function(x) {
  if (is.logical(x)) {
    return(ifelse(x, "T", "F"))
  }

  x <- as.character(x)
  x <- toupper(trimws(x))
  if (!x %in% c("T", "F")) {
    stop("Flag values must be logical or one of 'T'/'F'.", call. = FALSE)
  }
  x
}

sorter2_require_file <- function(path, label = "file") {
  if (!file.exists(path)) {
    stop(sprintf("%s does not exist: %s", label, path), call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

sorter2_require_dir <- function(path, label = "directory") {
  if (!dir.exists(path)) {
    stop(sprintf("%s does not exist: %s", label, path), call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

sorter2_with_trailing_slash <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (grepl("/$", path)) {
    return(path)
  }
  paste0(path, "/")
}

sorter2_default_script_dir <- function() {
  installed <- system.file("python", package = "sorter2r")
  if (!nzchar(installed)) {
    stop(
      "Could not find vendored scripts. Set script_dir explicitly ",
      "when running from source.",
      call. = FALSE
    )
  }
  installed
}

#' Reformat a reference FASTA for SORTER2
#'
#' Rewrites sequence IDs to the `>L{n}_{Family-Genus-species}` format
#' required by Stage 1B. The locus number `n` is taken from the `L{n}`
#' substring after the last `-` in the original ID; sequences whose
#' post-hyphen tag does not already carry an `L{n}_` prefix are assigned
#' sequential numbers starting after the highest embedded one. The
#' taxonomic label (fields 2–4 of the original `_`-delimited name) is
#' appended with `-` as a separator, ensuring every output ID is unique.
#'
#' @param input Path to the input reference FASTA file.
#' @param output Path for the reformatted output file. Defaults to the input
#'   path with `_sorter2` appended before the extension.
#' @return Path to the output file.
#' @export
sorter2_format_ref <- function(input, output = NULL) {
  input <- sorter2_require_file(input, "input")

  if (is.null(output)) {
    output <- paste0(
      tools::file_path_sans_ext(input), "_sorter2.fasta"
    )
  }

  lines <- readLines(input)
  header_idx <- grep("^>", lines)

  raw_names <- sub("^>", "", lines[header_idx])

  # Post-hyphen substring identifies the locus group
  locus_tags <- sub(".*-", "", raw_names)

  # Assign locus numbers: preserve embedded L{n} where present, else sequential
  unique_tags <- unique(locus_tags)
  has_ln <- grepl("^L[0-9]+_", unique_tags)

  existing_nums <- as.integer(
    sub("^L([0-9]+)_.*", "\\1", unique_tags[has_ln])
  )
  next_n <- if (length(existing_nums)) max(existing_nums) else 0L

  locus_num <- integer(length(unique_tags))
  names(locus_num) <- unique_tags
  locus_num[has_ln] <- as.integer(
    sub("^L([0-9]+)_.*", "\\1", unique_tags[has_ln])
  )
  for (tag in unique_tags[!has_ln]) {
    next_n <- next_n + 1L
    locus_num[[tag]] <- next_n
  }

  # Taxonomic label: fields 2-4 of the original name (family, genus, species)
  taxon <- vapply(strsplit(raw_names, "_"), function(parts) {
    paste(parts[2:4], collapse = "-")
  }, character(1))

  lines[header_idx] <- paste0(
    ">L", locus_num[locus_tags], "_", taxon
  )

  writeLines(lines, output)
  normalizePath(output, winslash = "/", mustWork = TRUE)
}

#' Run a SORTER2 Python script
#'
#' @param script Python script filename relative to the script directory.
#' @param args Character vector of command-line arguments.
#' @param python Python executable to use.
#' @param conda_env Optional conda environment name. If set, the command is run
#'   with `conda run -n conda_env`.
#' @param conda Conda executable to use when `conda_env` is set.
#' @param script_dir Directory containing the vendored SORTER2 scripts.
#' @param working_dir Optional working directory for the command.
#' @param dry_run If `TRUE`, return the command without executing it.
#' @param echo If `TRUE`, print the command before running it.
#' @param success_codes Integer vector of exit codes treated as success.
#'   Defaults to `0L`. Use `c(0L, 1L)` for scripts that exit with 1 on success.
#' @return A list describing the command execution.
#' @export
sorter2_run <- function(
  script,
  args = character(),
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  conda_env = Sys.getenv("SORTER2R_CONDA_ENV", ""),
  conda = Sys.getenv("SORTER2R_CONDA", "conda"),
  script_dir = NULL,
  working_dir = NULL,
  dry_run = FALSE,
  echo = TRUE,
  success_codes = 0L
) {
  if (is.null(script_dir)) {
    script_dir <- sorter2_default_script_dir()
  }
  script_dir <- sorter2_require_dir(script_dir, "script_dir")

  script_path <- file.path(script_dir, script)
  sorter2_require_file(script_path, "script")

  if (!is.null(working_dir)) {
    working_dir <- sorter2_require_dir(working_dir, "working_dir")
  }

  cmd_args <- c(script_path, args)
  if (is.null(conda_env)) {
    conda_env <- ""
  } else {
    conda_env <- as.character(conda_env)
    if (length(conda_env) != 1L) {
      stop("conda_env must be a single environment name.", call. = FALSE)
    }
    conda_env <- trimws(conda_env)
  }
  use_conda <- nzchar(conda_env)

  if (use_conda) {
    conda <- as.character(conda)
    if (length(conda) != 1L || !nzchar(trimws(conda))) {
      stop("conda must be a single executable path or command.", call. = FALSE)
    }
    command <- conda
    command_args <- c("run", "-n", conda_env, python, cmd_args)
  } else {
    command <- python
    command_args <- cmd_args
  }

  command_string <- paste(c(shQuote(command), shQuote(command_args)), collapse = " ")

  if (echo) {
    message("Running: ", command_string)
  }

  if (dry_run) {
    return(list(
      success = TRUE,
      status = 0L,
      command = command_string,
      stdout = character(),
      stderr = character(),
      dry_run = TRUE
    ))
  }

  stdout_file <- tempfile("sorter2r-stdout-")
  stderr_file <- tempfile("sorter2r-stderr-")
  on.exit(unlink(c(stdout_file, stderr_file), force = TRUE), add = TRUE)

  old_working_dir <- getwd()
  on.exit(setwd(old_working_dir), add = TRUE)
  if (!is.null(working_dir)) {
    setwd(working_dir)
  }

  exit_code <- system2(
    command = command,
    args = command_args,
    stdout = stdout_file,
    stderr = stderr_file
  )

  stdout <- if (file.exists(stdout_file)) {
    readLines(stdout_file)
  } else {
    character()
  }
  stderr <- if (file.exists(stderr_file)) {
    readLines(stderr_file)
  } else {
    character()
  }

  if (echo) {
    if (length(stdout) > 0) message(paste(stdout, collapse = "\n"))
    if (length(stderr) > 0) message(paste(stderr, collapse = "\n"))
  }

  list(
    success = as.integer(exit_code) %in% as.integer(success_codes),
    status = as.integer(exit_code),
    command = command_string,
    stdout = stdout,
    stderr = stderr,
    dry_run = FALSE
  )
}
