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

sorter2_run <- function(
  script,
  args = character(),
  python = Sys.getenv("SORTER2R_PYTHON", "python"),
  script_dir = NULL,
  working_dir = NULL,
  dry_run = FALSE,
  echo = TRUE
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
  command_string <- paste(c(shQuote(python), shQuote(cmd_args)), collapse = " ")

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

  exit_code <- system2(
    command = python,
    args = cmd_args,
    stdout = stdout_file,
    stderr = stderr_file,
    wd = working_dir
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

  list(
    success = identical(exit_code, 0L),
    status = as.integer(exit_code),
    command = command_string,
    stdout = stdout,
    stderr = stderr,
    dry_run = FALSE
  )
}
