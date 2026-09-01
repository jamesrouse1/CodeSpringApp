# Durable Slurm submission helpers for confirmed nf-core/sarek runs.
# The Shiny module owns presentation; this file owns filesystem and scheduler work.

SAREK_SUBMISSION_VERSION <- "0.5.0"

sarek_submission_require_helpers <- function() {
  required <- c(
    "sarek_text",
    "sarek_identifier",
    "sarek_is_absolute_path",
    "sarek_normalize_sex",
    "sarek_write_manifest",
    "sarek_write_nextflow_samplesheet"
  )
  missing <- required[!vapply(required, exists, logical(1), mode = "function")]
  if (length(missing)) {
    stop(
      "Source R/sarek_manifest.R and R/sarek_nextflow_input.R before R/sarek_submission.R. Missing helpers: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

sarek_submission_requested_tools <- function(analysis_mode, preset = "core") {
  analysis_mode <- sarek_text(analysis_mode)
  preset <- sarek_text(preset, "core")
  if (!identical(preset, "core")) {
    stop("Unsupported Sarek analysis preset: ", preset, call. = FALSE)
  }
  switch(
    analysis_mode,
    germline = "haplotypecaller",
    tumor_only = "mutect2,controlfreec",
    matched_tumor_normal = "mutect2,manta,ascat",
    annotation_only = "vep",
    stop("Unsupported Sarek analysis mode: ", analysis_mode, call. = FALSE)
  )
}

sarek_submission_missing_sex_patients <- function(manifest) {
  if (is.null(manifest) || !is.list(manifest) || !length(manifest$patients)) return(character(0))
  missing <- vapply(manifest$patients, function(patient) {
    !identical(sarek_normalize_sex(patient$sex, "NA"), "XX") &&
      !identical(sarek_normalize_sex(patient$sex, "NA"), "XY")
  }, logical(1))
  patient_ids <- vapply(manifest$patients, function(patient) {
    sarek_text(patient$patient_id, "unnamed patient")
  }, character(1))
  unique(patient_ids[missing])
}

sarek_submission_tool_resolution <- function(analysis_mode, preset = "core", manifest = NULL) {
  analysis_mode <- sarek_text(analysis_mode)
  requested <- sarek_submission_requested_tools(analysis_mode, preset)
  selected <- strsplit(requested, ",", fixed = TRUE)[[1]]
  missing_sex_patients <- if (is.null(manifest)) {
    character(0)
  } else {
    sarek_submission_missing_sex_patients(manifest)
  }
  sex_dependent <- switch(
    analysis_mode,
    tumor_only = "controlfreec",
    matched_tumor_normal = "ascat",
    character(0)
  )
  skipped <- if (length(missing_sex_patients)) intersect(selected, sex_dependent) else character(0)
  selected <- setdiff(selected, skipped)
  warnings <- if (length(skipped)) {
    paste0(
      "Sex chromosomes are not provided for patient(s): ",
      paste(missing_sex_patients, collapse = ", "),
      ". ", paste(toupper(skipped), collapse = ", "),
      " will be skipped for this run because Sarek requires sex for that caller. ",
      "The run will continue with: ", paste(selected, collapse = ", "), "."
    )
  } else character(0)
  list(
    requested_tools = requested,
    tools = paste(selected, collapse = ","),
    skipped_tools = skipped,
    missing_sex_patients = missing_sex_patients,
    warnings = warnings
  )
}

sarek_submission_tools <- function(analysis_mode, preset = "core", manifest = NULL) {
  sarek_submission_tool_resolution(analysis_mode, preset, manifest)$tools
}

sarek_submission_paths <- function(manifest) {
  manifest_id <- sarek_identifier(manifest$manifest_id, "sarek_analysis")
  results_root <- sarek_text(manifest$storage$results_root)
  work_root <- sarek_text(manifest$storage$work_root)
  if (!sarek_is_absolute_path(results_root) || !sarek_is_absolute_path(work_root)) {
    stop("Confirmed Sarek results and work roots must be absolute paths.", call. = FALSE)
  }
  run_dir <- normalizePath(file.path(results_root, manifest_id), winslash = "/", mustWork = FALSE)
  internal_dir <- file.path(run_dir, ".codespring")
  list(
    manifest_id = manifest_id,
    run_dir = run_dir,
    output_dir = file.path(run_dir, "results"),
    input_dir = file.path(internal_dir, "input"),
    internal_dir = internal_dir,
    log_dir = file.path(internal_dir, "logs"),
    work_dir = normalizePath(file.path(work_root, manifest_id), winslash = "/", mustWork = FALSE),
    manifest_path = file.path(internal_dir, "manifest.json"),
    samplesheet_path = file.path(internal_dir, "input", "samplesheet.csv"),
    params_path = file.path(internal_dir, "params.json"),
    run_config = file.path(internal_dir, "nextflow.config"),
    launch_script = file.path(internal_dir, "launch.sh"),
    nextflow_log = file.path(internal_dir, "logs", "nextflow.log"),
    trace_path = file.path(internal_dir, "logs", "trace.tsv"),
    stdout = file.path(internal_dir, "logs", "controller.out"),
    stderr = file.path(internal_dir, "logs", "controller.err"),
    submission_record = file.path(internal_dir, "submission.tsv"),
    controller_attempts = file.path(internal_dir, "controller_attempts.tsv"),
    runtime_status = file.path(internal_dir, "runtime_status.tsv"),
    active_children = file.path(internal_dir, "active_children.tsv"),
    controller_info = file.path(internal_dir, "controller.tsv")
  )
}

sarek_submission_params <- function(manifest, nextflow_input, paths) {
  assay <- sarek_text(manifest$assay$type)
  if (!assay %in% c("WGS", "annotation_only")) {
    stop(
      "Automated submission currently supports WGS and annotation-only runs. ",
      "WES and targeted runs require an intervals field before they can be submitted safely.",
      call. = FALSE
    )
  }
  list(
    input = paths$samplesheet_path,
    outdir = paths$output_dir,
    step = sarek_text(nextflow_input$step),
    genome = sarek_text(manifest$reference$sarek_genome),
    tools = sarek_submission_tools(manifest$analysis$mode, manifest$analysis$preset, manifest)
  )
}

sarek_submission_validate_runtime <- function(launcher, config, nxf_home, singularity_cache, sbatch) {
  values <- list(
    launcher = sarek_text(launcher),
    config = sarek_text(config),
    nxf_home = sarek_text(nxf_home),
    singularity_cache = sarek_text(singularity_cache),
    sbatch = sarek_text(sbatch, "sbatch")
  )
  if (!sarek_is_absolute_path(values$launcher) || !file.exists(values$launcher) || dir.exists(values$launcher)) {
    stop("The configured Sarek Nextflow executable is missing: ", values$launcher, call. = FALSE)
  }
  if (file.access(values$launcher, mode = 1) != 0) {
    stop("The configured Sarek Nextflow executable is not executable: ", values$launcher, call. = FALSE)
  }
  if (!sarek_is_absolute_path(values$config) || !file.exists(values$config) || dir.exists(values$config)) {
    stop("The configured CSHL Sarek profile is missing: ", values$config, call. = FALSE)
  }
  if (file.access(values$config, mode = 4) != 0) {
    stop("The configured CSHL Sarek profile is not readable: ", values$config, call. = FALSE)
  }
  for (field in c("nxf_home", "singularity_cache")) {
    if (!sarek_is_absolute_path(values[[field]])) {
      stop(field, " must be an absolute path.", call. = FALSE)
    }
  }
  sbatch_path <- if (sarek_is_absolute_path(values$sbatch)) values$sbatch else Sys.which(values$sbatch)
  if (!nzchar(sbatch_path) || !file.exists(sbatch_path) || file.access(sbatch_path, mode = 1) != 0) {
    stop("The Slurm sbatch command is unavailable: ", values$sbatch, call. = FALSE)
  }
  values$sbatch <- normalizePath(sbatch_path, winslash = "/", mustWork = TRUE)
  values
}

sarek_write_text_atomic <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".sarek_write_", tmpdir = dirname(path))
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  writeLines(as.character(lines), temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) stop("Could not atomically write: ", path, call. = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

sarek_submission_launch_script <- function(paths, runtime, pipeline, pipeline_version, nextflow_version) {
  exports <- c(
    paste("export NXF_HOME=", shQuote(runtime$nxf_home), sep = ""),
    paste("export NXF_SINGULARITY_CACHEDIR=", shQuote(runtime$singularity_cache), sep = ""),
    paste("export NXF_VER=", shQuote(nextflow_version), sep = "")
  )

  command <- c(
    shQuote(runtime$launcher),
    "-log", shQuote(paths$nextflow_log),
    "-c", shQuote(paths$run_config),
    "run", shQuote(pipeline),
    "-ansi-log false",
    "-r", shQuote(pipeline_version),
    "-profile", "singularity",
    "-params-file", shQuote(paths$params_path),
    "-work-dir", shQuote(paths$work_dir),
    "-with-trace", shQuote(paths$trace_path)
  )

  c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    exports,
    paste("mkdir -p", shQuote(paths$work_dir)),
    paste("cd", shQuote(paths$run_dir)),
    paste("status_file=", shQuote(paths$runtime_status), sep = ""),
    paste("active_children_file=", shQuote(paths$active_children), sep = ""),
    paste("controller_info_file=", shQuote(paths$controller_info), sep = ""),
    paste("child_tag=", shQuote(paths$child_tag), sep = ""),
    'controller_mode="${CSL_SAREK_CONTROLLER_MODE:-initial}"',
    "started_at=$(date -u +%Y-%m-%dT%H:%M:%SZ)",
    "controller_host=$(uname -n)",
    "{",
    "  printf 'field\\tvalue\\n'",
    "  printf 'pid\\t%s\\n' \"$$\"",
    "  printf 'host\\t%s\\n' \"$controller_host\"",
    "  printf 'started_at\\t%s\\n' \"$started_at\"",
    "} > \"${controller_info_file}.tmp\"",
    "mv -f \"${controller_info_file}.tmp\" \"$controller_info_file\"",
    "",
    'if [[ "${CSL_SAREK_RESUME:-false}" == "true" ]]; then',
    "  resume_args=(-resume)",
    "else",
    paste0(
      "  resume_args=(-name ",
      shQuote(paste0("codespring_", paths$manifest_id)),
      ")"
    ),
    "fi",
    "",
    "write_runtime_status() {",
    '  local state="$1"',
    '  local exit_code="$2"',
    '  local ended_at="$3"',
    '  local active_children="$4"',
    "  {",
    "    printf 'field\\tvalue\\n'",
    '    printf \'state\\t%s\\n\' "$state"',
    '    printf \'mode\\t%s\\n\' "$controller_mode"',
    '    printf \'started_at\\t%s\\n\' "$started_at"',
    '    printf \'ended_at\\t%s\\n\' "$ended_at"',
    '    printf \'exit_code\\t%s\\n\' "$exit_code"',
    '    printf \'active_children\\t%s\\n\' "$active_children"',
    '    printf \'child_tag\\t%s\\n\' "$child_tag"',
    '  } > "${status_file}.tmp"',
    '  mv -f "${status_file}.tmp" "$status_file"',
    "}",
    "",
    "snapshot_active_children() {",
    '  local raw="${active_children_file}.raw.$$"',
    '  local tmp="${active_children_file}.tmp.$$"',
    "",
    '  if ! squeue -h -u "$USER" -o "%i|%T|%M|%l|%k|%j" > "$raw" 2>/dev/null; then',
    '    rm -f "$raw" "$tmp"',
    "    return 1",
    "  fi",
    "",
    "  {",
    "    printf 'job_id\\tstate\\telapsed\\ttime_limit\\tcomment\\tjob_name\\n'",
    '    awk -F"|" -v tag="$child_tag" \'$5 == tag {',
    '      printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n", $1, $2, $3, $4, $5, $6',
    "    }' \"$raw\"",
    '  } > "$tmp"',
    "",
    '  mv -f "$tmp" "$active_children_file"',
    '  rm -f "$raw"',
    "  return 0",
    "}",
    "",
    "finish_controller() {",
    "  local exit_code=$?",
    "  trap - EXIT",
    "",
    "  local state=FAILED",
    '  local active_count=""',
    "  local child_query_ok=0",
    "",
    "  if snapshot_active_children; then",
    "    child_query_ok=1",
    '    active_count=$(awk \'NR > 1 { n++ } END { print n + 0 }\' "$active_children_file")',
    "  fi",
    "",
    '  if [ "$exit_code" -eq 0 ]; then',
    '    if [ "$child_query_ok" -eq 1 ] && [ "${active_count:-0}" -eq 0 ]; then',
    "      state=COMPLETED",
    "    else",
    "      state=INCOMPLETE",
    "    fi",
    "  fi",
    "",
    '  write_runtime_status "$state" "$exit_code" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${active_count:-unknown}"',
    '  exit "$exit_code"',
    "}",
    "",
    "trap finish_controller EXIT",
    "write_runtime_status RUNNING '' '' ''",
    paste(
      paste(command, collapse = " "),
      '"${resume_args[@]}"'
    )
  )
}

sarek_read_key_value_file <- function(path) {
  if (!nzchar(sarek_text(path)) || !file.exists(path) || dir.exists(path)) return(character(0))
  lines <- tryCatch(readLines(path, warn = FALSE), error = function(error) character(0))
  if (!length(lines)) return(character(0))
  pieces <- strsplit(lines, "\t", fixed = TRUE)
  keys <- vapply(pieces, function(parts) if (length(parts)) trimws(parts[[1]]) else "", character(1))
  values <- vapply(pieces, function(parts) {
    if (length(parts) < 2L) return("")
    paste(parts[-1L], collapse = "\t")
  }, character(1))
  keep <- nzchar(keys) & keys != "field"
  if (!any(keep)) return(character(0))
  stats::setNames(values[keep], keys[keep])
}

sarek_tail_text_lines <- function(path, max_lines = 200L, max_bytes = 262144L) {
  path <- sarek_text(path)
  max_lines <- suppressWarnings(as.integer(max_lines)[1])
  max_bytes <- suppressWarnings(as.numeric(max_bytes)[1])
  if (!length(max_lines) || is.na(max_lines) || max_lines < 1L) max_lines <- 200L
  if (!length(max_bytes) || is.na(max_bytes) || max_bytes < 1024) max_bytes <- 262144
  if (!nzchar(path) || !file.exists(path) || dir.exists(path) || file.access(path, mode = 4) != 0) {
    return(character(0))
  }
  size <- suppressWarnings(as.numeric(file.info(path)$size[[1]]))
  if (!is.finite(size) || size <= 0) return(character(0))
  start <- max(0, size - max_bytes)
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  if (start > 0) seek(connection, where = start, origin = "start")
  raw <- readBin(connection, what = "raw", n = as.integer(min(size, max_bytes)))
  if (!length(raw)) return(character(0))
  text <- tryCatch(rawToChar(raw), error = function(error) "")
  if (!nzchar(text)) return(character(0))
  text <- iconv(text, from = "UTF-8", to = "UTF-8", sub = "")
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  if (start > 0 && length(lines)) lines <- lines[-1L]
  lines <- sub("\r$", "", lines)
  utils::tail(lines[nzchar(trimws(lines))], max_lines)
}

sarek_read_nextflow_trace <- function(path, max_rows = 20000L) {
  empty <- data.frame(stringsAsFactors = FALSE)
  path <- sarek_text(path)
  if (!nzchar(path) || !file.exists(path) || dir.exists(path) || file.access(path, mode = 4) != 0) return(empty)
  trace <- tryCatch(
    utils::read.delim(
      path,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      fill = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    error = function(error) empty
  )
  if (!NROW(trace)) return(empty)
  names(trace) <- tolower(gsub("[^A-Za-z0-9]+", "_", names(trace)))
  max_rows <- suppressWarnings(as.integer(max_rows)[1])
  if (!length(max_rows) || is.na(max_rows) || max_rows < 1L) max_rows <- 20000L
  if (NROW(trace) > max_rows) trace <- utils::tail(trace, max_rows)
  trace
}

sarek_duration_seconds <- function(value) {
  values <- as.character(value)
  vapply(values, function(item) {
    item <- trimws(tolower(item))
    if (!nzchar(item) || item %in% c("na", "-") ) return(NA_real_)
    if (grepl("^[0-9]+(?:[.][0-9]+)?$", item, perl = TRUE)) return(suppressWarnings(as.numeric(item)))
    if (grepl("^[0-9]+:[0-5]?[0-9]:[0-5]?[0-9](?:[.][0-9]+)?$", item, perl = TRUE)) {
      parts <- suppressWarnings(as.numeric(strsplit(item, ":", fixed = TRUE)[[1]]))
      return(parts[[1]] * 3600 + parts[[2]] * 60 + parts[[3]])
    }
    matches <- gregexpr("[0-9]+(?:[.][0-9]+)?[[:space:]]*(?:ms|d|h|m|s)", item, perl = TRUE)
    tokens <- regmatches(item, matches)[[1]]
    if (!length(tokens) || identical(tokens, "")) return(NA_real_)
    numbers <- suppressWarnings(as.numeric(sub("[[:space:]]*(ms|d|h|m|s)$", "", tokens, perl = TRUE)))
    units <- sub("^.*[0-9][[:space:]]*", "", tokens, perl = TRUE)
    multipliers <- c(ms = 0.001, s = 1, m = 60, h = 3600, d = 86400)
    if (anyNA(numbers) || any(!units %in% names(multipliers))) return(NA_real_)
    sum(numbers * unname(multipliers[units]))
  }, numeric(1))
}

sarek_format_duration <- function(seconds) {
  seconds <- suppressWarnings(as.numeric(seconds)[1])
  if (!length(seconds) || !is.finite(seconds) || seconds < 0) return("")
  seconds <- round(seconds)
  if (seconds < 60) return(paste0(seconds, " sec"))
  if (seconds < 3600) return(paste0(round(seconds / 60), " min"))
  if (seconds < 86400) {
    hours <- floor(seconds / 3600)
    minutes <- round((seconds %% 3600) / 60)
    return(paste0(hours, " hr", if (hours == 1L) "" else "s", if (minutes > 0) paste0(" ", minutes, " min") else ""))
  }
  days <- floor(seconds / 86400)
  hours <- round((seconds %% 86400) / 3600)
  paste0(days, " day", if (days == 1L) "" else "s", if (hours > 0) paste0(" ", hours, " hr") else "")
}

sarek_trace_column <- function(trace, candidates, default = "") {
  selected <- candidates[candidates %in% names(trace)]
  if (!length(selected)) return(rep(default, NROW(trace)))
  as.character(trace[[selected[[1]]]])
}

sarek_process_group <- function(name) {
  trimws(sub("[[:space:]]*[(][^()]*[)][[:space:]]*$", "", as.character(name), perl = TRUE))
}

sarek_process_label <- function(name) {
  name <- trimws(as.character(name))
  sub("^.*:", "", name)
}

sarek_nextflow_log_events <- function(lines, max_events = 12L) {
  if (!length(lines)) return(list(events = character(0), submitted = character(0)))
  events <- character(0)
  submitted <- character(0)
  for (line in lines) {
    process_match <- regexec("(Submitted|Cached) process >[[:space:]]*(.+)$", line, perl = TRUE)
    process_parts <- regmatches(line, process_match)[[1]]
    if (length(process_parts) >= 3L) {
      kind <- process_parts[[2]]
      process <- trimws(process_parts[[3]])
      events <- c(events, paste0(kind, ": ", sarek_process_label(process)))
      if (identical(kind, "Submitted")) submitted <- c(submitted, process)
      next
    }
    if (grepl("ERROR[[:space:]]*~|Session aborted|Execution complete|Workflow completed|WARN[[:space:]]*~", line, perl = TRUE)) {
      message <- sub("^.* - ", "", line)
      events <- c(events, trimws(message))
    }
  }
  list(events = utils::tail(events, max_events), submitted = submitted)
}

sarek_multiset_difference <- function(values, remove) {
  remaining <- as.character(values)
  for (item in as.character(remove)) {
    position <- match(item, remaining)
    if (!is.na(position)) remaining <- remaining[-position]
  }
  remaining
}

sarek_task_time_estimate <- function(trace, active_names, now = Sys.time(), min_completed = 2L) {
  unavailable <- list(available = FALSE, label = "Not enough completed tasks of this step to estimate its duration.")
  if (!NROW(trace) || !length(active_names)) return(unavailable)
  names_value <- sarek_trace_column(trace, c("name", "process"))
  status <- toupper(sarek_trace_column(trace, "status"))
  duration_value <- sarek_trace_column(trace, c("realtime", "duration"), default = NA_character_)
  duration <- sarek_duration_seconds(duration_value)
  groups <- sarek_process_group(names_value)
  active <- as.character(active_names[[length(active_names)]])
  active_group <- sarek_process_group(active)
  completed <- status %in% c("COMPLETED", "CACHED") & groups == active_group & is.finite(duration) & duration > 0
  completed_count <- sum(completed)
  min_completed <- suppressWarnings(as.integer(min_completed)[1])
  if (!length(min_completed) || is.na(min_completed) || min_completed < 1L) min_completed <- 2L
  if (completed_count < min_completed) return(unavailable)
  expected <- stats::median(duration[completed])

  active_rows <- which(status %in% c("RUNNING", "SUBMITTED", "NEW") & groups == active_group)
  elapsed <- NA_real_
  if (length(active_rows)) {
    start_value <- sarek_trace_column(trace[active_rows, , drop = FALSE], c("start", "submit"))
    starts <- suppressWarnings(as.POSIXct(
      start_value,
      tz = "UTC",
      tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%OS")
    ))
    starts <- starts[!is.na(starts)]
    if (length(starts)) elapsed <- max(0, as.numeric(difftime(now, utils::tail(starts, 1L), units = "secs")))
  }
  remaining <- if (is.finite(elapsed)) max(0, expected - elapsed) else NA_real_
  detail <- paste0("median of ", completed_count, " completed task", if (completed_count == 1L) "" else "s", " of this step")
  label <- if (is.finite(remaining)) {
    paste0("Estimated time remaining for the current task: about ", sarek_format_duration(remaining), " (", detail, ").")
  } else {
    paste0("Typical duration for this task: about ", sarek_format_duration(expected), " (", detail, ").")
  }
  list(
    available = TRUE,
    label = label,
    expected_seconds = expected,
    elapsed_seconds = elapsed,
    remaining_seconds = remaining,
    completed_examples = completed_count
  )
}

sarek_run_activity <- function(run, state = "", max_events = 12L) {
  scalar <- function(name, default = "") {
    value <- tryCatch(run[[name]], error = function(error) NULL)
    sarek_text(value, default)
  }
  run_dir <- scalar("run_dir")
  internal_dir <- if (nzchar(run_dir)) file.path(run_dir, ".codespring") else ""
  log_path <- scalar("nextflow_log", if (nzchar(internal_dir)) file.path(internal_dir, "logs", "nextflow.log") else "")
  trace_path <- scalar("trace_path", if (nzchar(internal_dir)) file.path(internal_dir, "logs", "trace.tsv") else "")
  trace <- sarek_read_nextflow_trace(trace_path)
  log_lines <- sarek_tail_text_lines(log_path, max_lines = 2000L)
  log <- sarek_nextflow_log_events(log_lines, max_events = max_events)
  trace_names <- sarek_trace_column(trace, c("name", "process"))
  trace_status <- toupper(sarek_trace_column(trace, "status"))
  terminal <- trace_status %in% c("COMPLETED", "CACHED", "FAILED", "ABORTED", "CANCELLED")
  trace_active <- trace_names[trace_status %in% c("RUNNING", "SUBMITTED", "NEW") & nzchar(trace_names)]
  inferred_active <- if (NROW(trace)) {
    sarek_multiset_difference(log$submitted, trace_names[terminal])
  } else {
    utils::tail(log$submitted, 1L)
  }
  active_names <- unique(c(trace_active, inferred_active))
  normalized_state <- sarek_normalize_slurm_state(state)
  if (nzchar(normalized_state) && !normalized_state %in% c("RUNNING", "COMPLETING", "SUSPENDED")) {
    active_names <- character(0)
  }
  estimate <- sarek_task_time_estimate(trace, active_names)
  completed <- sum(trace_status %in% c("COMPLETED", "CACHED"))
  failed <- sum(trace_status %in% c("FAILED", "ABORTED", "CANCELLED"))
  running <- max(sum(trace_status %in% c("RUNNING", "SUBMITTED", "NEW")), length(active_names))
  list(
    available = length(log_lines) > 0L || NROW(trace) > 0L,
    current_steps = unname(utils::tail(active_names, 5L)),
    current_labels = unname(vapply(utils::tail(active_names, 5L), sarek_process_label, character(1))),
    completed = completed,
    failed = failed,
    running = running,
    observed = completed + failed + running,
    recent_events = log$events,
    estimate = estimate,
    log_path = log_path,
    trace_path = trace_path,
    trace_available = NROW(trace) > 0L
  )
}

sarek_submission_catalog <- function(results_root) {
  empty <- data.frame(
    run_id = character(0), status = character(0), job_id = character(0),
    controller_pid = character(0), controller_host = character(0),
    child_tag = character(0),
    submitted_at = character(0), step = character(0), tools = character(0),
    run_dir = character(0), output_dir = character(0), work_dir = character(0),
    runtime_status = character(0), updated_at = as.POSIXct(character(0)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  results_root <- sarek_text(results_root)
  if (!sarek_is_absolute_path(results_root) || !dir.exists(results_root) || file.access(results_root, mode = 5) != 0) {
    return(empty)
  }
  run_dirs <- list.dirs(results_root, recursive = FALSE, full.names = TRUE)
  records <- file.path(run_dirs, ".codespring", "submission.tsv")
  keep <- file.exists(records)
  run_dirs <- run_dirs[keep]
  records <- records[keep]
  if (!length(run_dirs)) return(empty)

  rows <- lapply(seq_along(run_dirs), function(index) {
    run_dir <- normalizePath(run_dirs[[index]], winslash = "/", mustWork = FALSE)
    values <- sarek_read_key_value_file(records[[index]])
    controller_info <- sarek_read_key_value_file(
      file.path(run_dir, ".codespring", "controller.tsv")
    )
    params_path <- file.path(run_dir, ".codespring", "params.json")
    params <- if (file.exists(params_path) && requireNamespace("jsonlite", quietly = TRUE)) {
      tryCatch(jsonlite::read_json(params_path, simplifyVector = TRUE), error = function(error) list())
    } else {
      list()
    }
    record_info <- file.info(records[[index]])
    data.frame(
      run_id = basename(run_dir),
      status = sarek_text(values["status"], "submitted"),
      job_id = sarek_text(values["job_id"]),
      controller_pid = sarek_text(
        controller_info["pid"],
        sarek_text(values["controller_pid"])
      ),
      controller_host = sarek_text(
        controller_info["host"],
        sarek_text(values["controller_host"])
      ),
      child_tag = sarek_text(values["child_tag"]),
      submitted_at = sarek_text(values["submitted_at"]),
      step = sarek_text(values["step"], sarek_text(params$step)),
      tools = sarek_text(values["tools"], sarek_text(params$tools)),
      run_dir = sarek_text(values["run_dir"], run_dir),
      output_dir = sarek_text(values["output_dir"], file.path(run_dir, "results")),
      work_dir = sarek_text(values["work_dir"]),
      runtime_status = file.path(run_dir, ".codespring", "runtime_status.tsv"),
      updated_at = record_info$mtime[[1]],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  result <- do.call(rbind, rows)
  order_value <- suppressWarnings(as.POSIXct(result$submitted_at, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  missing_time <- is.na(order_value)
  order_value[missing_time] <- result$updated_at[missing_time]
  result[order(order_value, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
}

sarek_result_file_type <- function(path) {
  name <- tolower(basename(sarek_text(path)))
  if (grepl("multiqc|report|timeline|dag", name) && grepl("\\.(html?|pdf|png|svg)$", name)) return("Report")
  if (grepl("\\.(vcf|vcf\\.gz|bcf)$", name)) return("Variants")
  if (grepl("\\.(bam|cram|sam)$", name)) return("Alignment")
  if (grepl("\\.(bai|crai|csi|tbi)$", name)) return("Index")
  if (grepl("\\.(bed|bed\\.gz|interval_list)$", name)) return("Regions")
  if (grepl("\\.(tsv|csv|txt|json|yaml|yml)$", name)) return("Table / metadata")
  if (grepl("\\.(html?|pdf|png|jpg|jpeg|svg)$", name)) return("Report")
  "Other"
}

sarek_format_file_size <- function(bytes) {
  bytes <- suppressWarnings(as.numeric(bytes))
  labels <- c("B", "KB", "MB", "GB", "TB")
  vapply(bytes, function(value) {
    if (is.na(value) || value < 0) return("")
    level <- if (value <= 0) 1L else min(length(labels), floor(log(value, 1024)) + 1L)
    scaled <- value / (1024 ^ (level - 1L))
    paste0(format(round(scaled, if (level <= 2L) 0L else 1L), trim = TRUE, scientific = FALSE), " ", labels[[level]])
  }, character(1))
}

sarek_result_file_catalog <- function(output_dir, max_files = 5000L) {
  empty <- data.frame(
    file = character(0), folder = character(0), type = character(0),
    size = character(0), modified = character(0), path = character(0),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  output_dir <- sarek_text(output_dir)
  max_files <- suppressWarnings(as.integer(max_files)[1])
  if (!length(max_files) || is.na(max_files) || max_files < 1L) max_files <- 5000L
  if (!sarek_is_absolute_path(output_dir) || !dir.exists(output_dir) || file.access(output_dir, mode = 5) != 0) {
    return(empty)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  queue <- output_dir
  files <- character(0)
  scanned_dirs <- 0L
  truncated <- FALSE
  while (length(queue) && length(files) < max_files && scanned_dirs < max_files) {
    directory <- queue[[1]]
    queue <- queue[-1L]
    scanned_dirs <- scanned_dirs + 1L
    entries <- tryCatch(
      list.files(directory, full.names = TRUE, recursive = FALSE, all.files = FALSE, no.. = TRUE),
      error = function(error) character(0)
    )
    if (!length(entries)) next
    readable_dirs <- entries[dir.exists(entries) & file.access(entries, mode = 5) == 0]
    queue <- c(queue, sort(readable_dirs))
    regular <- entries[!dir.exists(entries) & file.exists(entries) & file.access(entries, mode = 4) == 0]
    remaining <- max_files - length(files)
    if (length(regular) > remaining) truncated <- TRUE
    files <- c(files, utils::head(sort(regular), remaining))
  }
  if (length(queue)) truncated <- TRUE
  if (!length(files)) {
    attr(empty, "truncated") <- truncated
    return(empty)
  }
  info <- file.info(files)
  relative <- substring(dirname(files), nchar(output_dir) + 1L)
  relative <- sub("^/", "", relative)
  relative[!nzchar(relative)] <- "."
  result <- data.frame(
    file = basename(files),
    folder = relative,
    type = vapply(files, sarek_result_file_type, character(1)),
    size = sarek_format_file_size(info$size),
    modified = format(info$mtime, "%Y-%m-%d %H:%M"),
    path = normalizePath(files, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  type_order <- match(result$type, c("Report", "Variants", "Alignment", "Index", "Regions", "Table / metadata", "Other"))
  result <- result[order(type_order, result$folder, result$file), , drop = FALSE]
  rownames(result) <- NULL
  attr(result, "truncated") <- truncated
  result
}

sarek_normalize_slurm_state <- function(state) {
  state <- toupper(trimws(sarek_text(state)))
  state <- sub("[+ ].*$", "", state)
  aliases <- c(
    PD = "PENDING", R = "RUNNING", CG = "COMPLETING", CD = "COMPLETED",
    F = "FAILED", CA = "CANCELLED", TO = "TIMEOUT", OOM = "OUT_OF_MEMORY",
    ERROR = "FAILED"
  )
  if (state %in% names(aliases)) unname(aliases[[state]]) else state
}

sarek_run_command <- function(command, args, runner = NULL) {
  if (is.function(runner)) return(runner(command, args))
  executable <- if (sarek_is_absolute_path(command)) command else Sys.which(command)
  if (!nzchar(executable) || !file.exists(executable)) return(character(0))
  output <- tryCatch(
    system2(executable, vapply(args, shQuote, character(1)), stdout = TRUE, stderr = TRUE),
    error = function(error) character(0)
  )
  status <- attr(output, "status")
  if (!is.null(status) && !identical(as.integer(status), 0L)) character(0) else output
}

sarek_query_slurm_job <- function(job_id, runner = NULL, squeue = "squeue", sacct = "sacct") {
  job_id <- sarek_text(job_id)
  if (!grepl("^[0-9]+$", job_id)) return(list(state = "", elapsed = "", source = "record"))
  queue_output <- sarek_run_command(
    squeue,
    c("-h", "-j", job_id, "-o", "%T|%M"),
    runner
  )
  queue_lines <- trimws(as.character(queue_output))
  queue_lines <- queue_lines[nzchar(queue_lines) & !startsWith(queue_lines, "slurm_")]
  if (length(queue_lines)) {
    fields <- strsplit(queue_lines[[1]], "|", fixed = TRUE)[[1]]
    return(list(
      state = sarek_normalize_slurm_state(fields[[1]]),
      elapsed = if (length(fields) >= 2L) trimws(fields[[2]]) else "",
      source = "squeue"
    ))
  }
  account_output <- sarek_run_command(
    sacct,
    c("-n", "-P", "-j", job_id, "--format=State,Elapsed"),
    runner
  )
  account_lines <- trimws(as.character(account_output))
  account_lines <- account_lines[nzchar(account_lines) & !startsWith(account_lines, "slurm_")]
  if (length(account_lines)) {
    fields <- strsplit(account_lines[[1]], "|", fixed = TRUE)[[1]]
    return(list(
      state = sarek_normalize_slurm_state(fields[[1]]),
      elapsed = if (length(fields) >= 2L) trimws(fields[[2]]) else "",
      source = "sacct"
    ))
  }
  list(state = "", elapsed = "", source = "record")
}

sarek_pid_alive <- function(pid, host = "") {
  pid <- sarek_text(pid)

  if (!grepl("^[0-9]+$", pid)) return(FALSE)

  host <- sarek_text(host)
  local_host <- sarek_text(Sys.info()[["nodename"]])

  if (nzchar(host) && nzchar(local_host) && !identical(host, local_host)) {
    return(NA)
  }

  file.exists(file.path("/proc", pid))
}

sarek_query_active_children <- function(
  child_tag,
  user = Sys.info()[["user"]],
  runner = NULL,
  squeue = "squeue"
) {
  empty <- data.frame(
    job_id = character(0),
    state = character(0),
    elapsed = character(0),
    time_limit = character(0),
    comment = character(0),
    job_name = character(0),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  child_tag <- sarek_text(child_tag)
  user <- sarek_text(user)

  if (!nzchar(child_tag) || !nzchar(user)) return(empty)

  output <- sarek_run_command(
    squeue,
    c("-h", "-u", user, "-o", "%i|%T|%M|%l|%k|%j"),
    runner = runner
  )

  lines <- trimws(as.character(output))
  lines <- lines[nzchar(lines)]

  if (!length(lines)) return(empty)

  fields <- strsplit(lines, "|", fixed = TRUE)

  rows <- lapply(fields, function(parts) {
    if (length(parts) < 6L || trimws(parts[[5]]) != child_tag) return(NULL)

    data.frame(
      job_id = trimws(parts[[1]]),
      state = sarek_normalize_slurm_state(parts[[2]]),
      elapsed = trimws(parts[[3]]),
      time_limit = trimws(parts[[4]]),
      comment = trimws(parts[[5]]),
      job_name = trimws(paste(parts[6:length(parts)], collapse = "|")),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(empty)

  do.call(rbind, rows)
}

sarek_run_status <- function(run, runner = NULL, squeue = "squeue", sacct = "sacct") {
  scalar <- function(name, default = "") {
    value <- tryCatch(run[[name]], error = function(error) NULL)
    sarek_text(value, default)
  }

  runtime <- sarek_read_key_value_file(scalar("runtime_status"))
  runtime_state <- toupper(sarek_text(runtime["state"]))

  run_dir <- scalar("run_dir")
  controller_info <- if (nzchar(run_dir)) {
    sarek_read_key_value_file(
      file.path(run_dir, ".codespring", "controller.tsv")
    )
  } else {
    character(0)
  }

  controller_pid <- sarek_text(
    controller_info["pid"],
    scalar("controller_pid")
  )
  controller_host <- sarek_text(
    controller_info["host"],
    scalar("controller_host")
  )
  child_tag <- scalar("child_tag", sarek_text(runtime["child_tag"]))

  children <- sarek_query_active_children(
    child_tag,
    runner = runner,
    squeue = squeue
  )

  active_children <- NROW(children)

  terminal <- c(
    "COMPLETED",
    "FAILED",
    "INCOMPLETE",
    "CANCELLED",
    "TIMEOUT",
    "OUT_OF_MEMORY"
  )

  controller_alive <- sarek_pid_alive(
    controller_pid,
    controller_host
  )

  state <- ""

  if (runtime_state %in% terminal) {
    state <- runtime_state
  } else if (isTRUE(controller_alive)) {
    state <- "RUNNING"
  } else if (identical(runtime_state, "RUNNING")) {
    state <- "INCOMPLETE"
  } else if (active_children > 0L) {
    state <- "INCOMPLETE"
  } else {
    state <- toupper(scalar("status", "SUBMITTED"))
  }

  list(
    state = state,
    elapsed = "",
    source = if (runtime_state %in% terminal) "runtime" else "controller",
    exit_code = sarek_text(runtime["exit_code"]),
    started_at = sarek_text(runtime["started_at"]),
    ended_at = sarek_text(runtime["ended_at"]),
    controller_pid = controller_pid,
    controller_alive = controller_alive,
    active_children = active_children,
    child_jobs = children
  )
}

sarek_run_progress <- function(state) {
  state <- sarek_normalize_slurm_state(state)
  if (state %in% c("PENDING", "CONFIGURING", "SUBMITTED")) {
    return(list(percent = 10L, label = "Queued", kind = "queued", active = TRUE))
  }
  if (state %in% c("RUNNING", "SUSPENDED")) {
    return(list(percent = 55L, label = "Running", kind = "running", active = TRUE))
  }
  if (identical(state, "COMPLETING")) {
    return(list(percent = 90L, label = "Finishing", kind = "running", active = TRUE))
  }
  if (identical(state, "COMPLETED")) {
    return(list(percent = 100L, label = "Completed", kind = "success", active = FALSE))
  }
  if (state %in% c("FAILED", "INCOMPLETE", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY")) {
    return(list(percent = 100L, label = gsub("_", " ", tools::toTitleCase(tolower(state))), kind = "error", active = FALSE))
  }
  list(percent = 15L, label = if (nzchar(state)) gsub("_", " ", state) else "Status unavailable", kind = "unknown", active = FALSE)
}

sarek_build_submission_bundle <- function(
  manifest,
  nextflow_input,
  runtime,
  nextflow_version = "25.10.2"
) {
  sarek_submission_require_helpers()
  if (!is.list(manifest) || !identical(sarek_text(manifest$status), "confirmed")) {
    stop("Only a confirmed Sarek manifest can be submitted.", call. = FALSE)
  }
  if (!is.list(nextflow_input) || !is.data.frame(nextflow_input$samplesheet) || !NROW(nextflow_input$samplesheet)) {
    stop("A valid generated Sarek samplesheet is required for submission.", call. = FALSE)
  }
  if (!identical(sarek_text(nextflow_input$pipeline), SAREK_NEXTFLOW_PIPELINE) ||
      !identical(sarek_text(nextflow_input$pipeline_version), SAREK_NEXTFLOW_PIPELINE_VERSION)) {
    stop("The generated input does not match the configured Sarek pipeline version.", call. = FALSE)
  }

  paths <- sarek_submission_paths(manifest)
  if (file.exists(paths$run_dir)) {
    stop(
      "A Sarek run named '", paths$manifest_id,
      "' already exists at ", paths$run_dir,
      ". Choose a new manifest ID; resume support will be added separately.",
      call. = FALSE
    )
  }
  paths$child_tag <- substr(
    paste0(
      "codespring_sarek_",
      paths$manifest_id,
      "_",
      format(Sys.time(), "%Y%m%d%H%M%S", tz = "UTC"),
      "_",
      Sys.getpid()
    ),
    1L,
    220L
  )

  params <- sarek_submission_params(manifest, nextflow_input, paths)
  tool_resolution <- sarek_submission_tool_resolution(
    manifest$analysis$mode,
    manifest$analysis$preset,
    manifest
  )
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The jsonlite R package is required to prepare a Sarek run.", call. = FALSE)
  }

  for (directory in c(paths$input_dir, paths$log_dir, paths$output_dir)) {
    if (!dir.create(directory, recursive = TRUE, showWarnings = FALSE) && !dir.exists(directory)) {
      stop("Could not create the Sarek run directory: ", directory, call. = FALSE)
    }
  }
  sarek_write_manifest(manifest, paths$manifest_path)
  sarek_write_nextflow_samplesheet(nextflow_input, paths$samplesheet_path)
  jsonlite::write_json(params, paths$params_path, auto_unbox = TRUE, pretty = TRUE, null = "null")

  sarek_write_text_atomic(
    c(
      paste("includeConfig", shQuote(runtime$config)),
      "",
      "process {",
      paste0(
        "  clusterOptions = '--comment=",
        paths$child_tag,
        "'"
      ),
      "}"
    ),
    paths$run_config
  )

  sarek_write_text_atomic(
    sarek_submission_launch_script(
      paths,
      runtime,
      nextflow_input$pipeline,
      nextflow_input$pipeline_version,
      nextflow_version
    ),
    paths$launch_script
  )
  Sys.chmod(paths$launch_script, mode = "0700")
  list(paths = paths, params = params, tool_resolution = tool_resolution)
}

sarek_parse_slurm_job_id <- function(output) {
  text <- trimws(paste(as.character(output), collapse = "\n"))
  first <- strsplit(text, "[;[:space:]]+", perl = TRUE)[[1]][[1]]
  if (grepl("^[0-9]+$", first)) first else ""
}




sarek_submission_write_values <- function(path, values) {
  values <- as.character(values)
  names(values) <- names(values)
  sarek_write_text_atomic(
    c(
      "field\tvalue",
      paste(names(values), values, sep = "\t")
    ),
    path
  )
}


sarek_validate_slurm_time <- function(value, default = "2-00:00:00") {
  value <- sarek_text(value, default)
  valid <- grepl(
    "^([0-9]+-([01]?[0-9]|2[0-3])|[0-9]+):[0-5][0-9]:[0-5][0-9]$",
    value,
    perl = TRUE
  )
  if (!valid) {
    stop(
      "Sarek controller_time compatibility value must use HH:MM:SS or D-HH:MM:SS format: ",
      value,
      call. = FALSE
    )
  }
  value
}

sarek_start_detached_controller <- function(
  paths,
  mode = c("initial", "resume"),
  starter = NULL
) {
  mode <- match.arg(mode)

  if (is.function(starter)) {
    output <- starter(paths, mode)
    pid <- sarek_parse_slurm_job_id(output)

    if (!nzchar(pid)) {
      stop(
        "Detached-controller test starter did not return a numeric PID.",
        call. = FALSE
      )
    }

    return(list(
      pid = pid,
      host = sarek_text(Sys.info()[["nodename"]]),
      output = as.character(output)
    ))
  }

  bash <- Sys.which("bash")
  nohup <- Sys.which("nohup")
  setsid <- Sys.which("setsid")

  if (!nzchar(bash) || !nzchar(nohup) || !nzchar(setsid)) {
    stop(
      "Detached Sarek execution requires bash, nohup, and setsid on the application host.",
      call. = FALSE
    )
  }

  resume_value <- if (identical(mode, "resume")) "true" else "false"

  command <- paste(
    shQuote(nohup),
    shQuote(setsid),
    "env",
    paste0("CSL_SAREK_RESUME=", resume_value),
    paste0("CSL_SAREK_CONTROLLER_MODE=", mode),
    shQuote(paths$launch_script),
    "</dev/null",
    paste0(">", shQuote(paths$stdout)),
    paste0("2>", shQuote(paths$stderr)),
    "& echo $!"
  )

  output <- tryCatch(
    system2(
      bash,
      c("-lc", shQuote(command)),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(error) structure(conditionMessage(error), status = 1L)
  )

  status <- attr(output, "status")
  if (is.null(status)) status <- 0L

  pid <- sarek_parse_slurm_job_id(output)

  if (!identical(as.integer(status), 0L) || !nzchar(pid)) {
    stop(
      "Could not start the detached Sarek Nextflow controller: ",
      paste(as.character(output), collapse = " "),
      call. = FALSE
    )
  }

  list(
    pid = pid,
    host = sarek_text(Sys.info()[["nodename"]]),
    output = as.character(output)
  )
}

sarek_submit_run <- function(
  manifest,
  nextflow_input,
  launcher,
  config,
  nxf_home,
  singularity_cache,
  queue = "cpuq",
  controller_time = "2-00:00:00",
  nextflow_version = "25.10.2",
  sbatch = "sbatch",
  submitter = NULL
) {
  runtime <- sarek_submission_validate_runtime(
    launcher,
    config,
    nxf_home,
    singularity_cache,
    sbatch
  )

  # Kept for API compatibility. The detached controller itself no longer
  # consumes a Slurm allocation.
  sarek_validate_slurm_time(controller_time)

  bundle <- sarek_build_submission_bundle(
    manifest,
    nextflow_input,
    runtime,
    nextflow_version
  )

  paths <- bundle$paths

  starter <- if (is.function(submitter)) {
    function(paths, mode) submitter("detached", c(mode, paths$launch_script))
  } else {
    NULL
  }

  controller <- sarek_start_detached_controller(
    paths,
    mode = "initial",
    starter = starter
  )

  submitted_at <- format(
    Sys.time(),
    "%Y-%m-%dT%H:%M:%SZ",
    tz = "UTC"
  )

  sarek_write_text_atomic(
    c(
      "field\tvalue",
      "status\tsubmitted",
      "job_id\t",
      paste0("controller_pid\t", controller$pid),
      paste0("controller_host\t", controller$host),
      "controller_mode\tdetached",
      "attempt\t1",
      "mode\tinitial",
      paste0("child_tag\t", paths$child_tag),
      paste0("submitted_at\t", submitted_at),
      paste0("run_dir\t", paths$run_dir),
      paste0("output_dir\t", paths$output_dir),
      paste0("work_dir\t", paths$work_dir),
      paste0("step\t", bundle$params$step),
      paste0("tools\t", bundle$params$tools),
      paste0(
        "requested_tools\t",
        bundle$tool_resolution$requested_tools
      ),
      paste0(
        "skipped_tools\t",
        paste(bundle$tool_resolution$skipped_tools, collapse = ",")
      ),
      paste0(
        "tool_warning\t",
        paste(bundle$tool_resolution$warnings, collapse = " ")
      ),
      "controller_time\tnot_applicable",
      paste0(
        "submission_version\t",
        SAREK_SUBMISSION_VERSION
      )
    ),
    paths$submission_record
  )

  list(
    status = "submitted",
    job_id = "",
    controller_pid = controller$pid,
    controller_host = controller$host,
    controller_mode = "detached",
    child_tag = paths$child_tag,
    attempt = 1L,
    mode = "initial",
    run_dir = paths$run_dir,
    output_dir = paths$output_dir,
    work_dir = paths$work_dir,
    tools = bundle$params$tools,
    requested_tools = bundle$tool_resolution$requested_tools,
    skipped_tools = bundle$tool_resolution$skipped_tools,
    warnings = bundle$tool_resolution$warnings,
    step = bundle$params$step,
    controller_time = "not_applicable"
  )
}

sarek_resume_run <- function(
  run_dir,
  launcher,
  config,
  nxf_home,
  singularity_cache,
  queue = "cpuq",
  controller_time = "2-00:00:00",
  sbatch = "sbatch",
  submitter = NULL,
  status_runner = NULL,
  squeue = "squeue",
  sacct = "sacct"
) {
  run_dir <- sarek_text(run_dir)

  if (!sarek_is_absolute_path(run_dir) || !dir.exists(run_dir)) {
    stop(
      "Cannot resume Sarek run because its run directory is missing: ",
      run_dir,
      call. = FALSE
    )
  }

  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)

  internal_dir <- file.path(run_dir, ".codespring")
  log_dir <- file.path(internal_dir, "logs")

  submission_record <- file.path(internal_dir, "submission.tsv")
  launch_script <- file.path(internal_dir, "launch.sh")
  params_path <- file.path(internal_dir, "params.json")
  run_config <- file.path(internal_dir, "nextflow.config")
  runtime_status <- file.path(internal_dir, "runtime_status.tsv")
  active_children <- file.path(internal_dir, "active_children.tsv")

  required <- c(
    submission_record,
    launch_script,
    params_path,
    run_config
  )

  missing <- required[!file.exists(required)]

  if (length(missing)) {
    stop(
      "Cannot resume Sarek because required run state is missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  launch_lines <- tryCatch(
    readLines(launch_script, warn = FALSE),
    error = function(error) character(0)
  )

  if (!any(grepl("CSL_SAREK_RESUME", launch_lines, fixed = TRUE))) {
    stop(
      "This Sarek run was created before detached resume support was added.",
      call. = FALSE
    )
  }

  values <- sarek_read_key_value_file(submission_record)

  controller_info <- sarek_read_key_value_file(
    file.path(internal_dir, "controller.tsv")
  )

  existing_pid <- sarek_text(
    controller_info["pid"],
    sarek_text(values["controller_pid"])
  )
  existing_host <- sarek_text(
    controller_info["host"],
    sarek_text(values["controller_host"])
  )

  if (isTRUE(sarek_pid_alive(existing_pid, existing_host))) {
    stop(
      "Cannot resume Sarek while its Nextflow controller PID ",
      existing_pid,
      " is still running.",
      call. = FALSE
    )
  }

  runtime <- sarek_submission_validate_runtime(
    launcher,
    config,
    nxf_home,
    singularity_cache,
    sbatch
  )

  sarek_validate_slurm_time(controller_time)

  work_dir <- sarek_text(values["work_dir"])

  if (!sarek_is_absolute_path(work_dir) || !dir.exists(work_dir)) {
    stop(
      "Cannot resume Sarek because its original work directory is missing: ",
      work_dir,
      call. = FALSE
    )
  }

  child_tag <- sarek_text(values["child_tag"])

  if (!nzchar(child_tag)) {
    stop(
      "Cannot resume Sarek because its child-job tag was not recorded.",
      call. = FALSE
    )
  }

  paths <- list(
    manifest_id = basename(run_dir),
    run_dir = run_dir,
    output_dir = sarek_text(
      values["output_dir"],
      file.path(run_dir, "results")
    ),
    internal_dir = internal_dir,
    log_dir = log_dir,
    work_dir = work_dir,
    params_path = params_path,
    run_config = run_config,
    launch_script = launch_script,
    nextflow_log = file.path(log_dir, "nextflow.log"),
    trace_path = file.path(log_dir, "trace.tsv"),
    stdout = file.path(log_dir, "controller.out"),
    stderr = file.path(log_dir, "controller.err"),
    submission_record = submission_record,
    runtime_status = runtime_status,
    active_children = active_children,
    controller_info = file.path(internal_dir, "controller.tsv"),
    child_tag = child_tag
  )

  starter <- if (is.function(submitter)) {
    function(paths, mode) submitter("detached", c(mode, paths$launch_script))
  } else {
    NULL
  }

  controller <- sarek_start_detached_controller(
    paths,
    mode = "resume",
    starter = starter
  )

  previous_attempt <- suppressWarnings(
    as.integer(sarek_text(values["attempt"], "1"))
  )

  if (is.na(previous_attempt) || previous_attempt < 1L) {
    previous_attempt <- 1L
  }

  attempt <- previous_attempt + 1L

  submitted_at <- format(
    Sys.time(),
    "%Y-%m-%dT%H:%M:%SZ",
    tz = "UTC"
  )

  values["status"] <- "submitted"
  values["job_id"] <- ""
  values["controller_pid"] <- controller$pid
  values["controller_host"] <- controller$host
  values["controller_mode"] <- "detached"
  values["attempt"] <- as.character(attempt)
  values["mode"] <- "resume"
  values["submitted_at"] <- submitted_at
  values["controller_time"] <- "not_applicable"
  values["submission_version"] <- SAREK_SUBMISSION_VERSION

  sarek_submission_write_values(
    submission_record,
    values
  )

  list(
    status = "submitted",
    job_id = "",
    controller_pid = controller$pid,
    controller_host = controller$host,
    controller_mode = "detached",
    child_tag = child_tag,
    attempt = attempt,
    mode = "resume",
    run_dir = run_dir,
    output_dir = paths$output_dir,
    work_dir = paths$work_dir,
    controller_time = "not_applicable"
  )
}
