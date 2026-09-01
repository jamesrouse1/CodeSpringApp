args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tests/test_sarek_live_activity.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

source(file.path(repo_root, "R", "sarek_manifest.R"), local = FALSE)
source(file.path(repo_root, "R", "sarek_nextflow_input.R"), local = FALSE)
source(file.path(repo_root, "R", "sarek_submission.R"), local = FALSE)

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

test_root <- tempfile("sarek-live-activity-")
run_dir <- file.path(test_root, "live_run")
log_dir <- file.path(run_dir, ".codespring", "logs")
dir.create(log_dir, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

log_path <- file.path(log_dir, "nextflow.log")
trace_path <- file.path(log_dir, "trace.tsv")
writeLines(
  c(
    "Aug-20 12:00:00 [Task submitter] INFO nextflow.Session - [aa/111111] Submitted process > NFCORE_SAREK:SAREK:ALIGN (sample_a)",
    "Aug-20 12:10:00 [Task submitter] INFO nextflow.Session - [bb/222222] Submitted process > NFCORE_SAREK:SAREK:ALIGN (sample_b)",
    "Aug-20 12:20:00 [Task submitter] INFO nextflow.Session - [cc/333333] Submitted process > NFCORE_SAREK:SAREK:ALIGN (sample_c)"
  ),
  log_path
)
writeLines(
  c(
    "task_id\tname\tstatus\trealtime\tstart",
    "1\tNFCORE_SAREK:SAREK:ALIGN (sample_a)\tCOMPLETED\t10m\t2026-08-20 12:00:00",
    "2\tNFCORE_SAREK:SAREK:ALIGN (sample_b)\tCOMPLETED\t14m\t2026-08-20 12:10:00",
    "3\tNFCORE_SAREK:SAREK:ALIGN (sample_c)\tRUNNING\t\t2026-08-20 12:20:00"
  ),
  trace_path
)

duration <- sarek_duration_seconds(c("420ms", "1m 30s", "1:02:03", "2h"))
assert_true(
  isTRUE(all.equal(unname(duration), c(0.42, 90, 3723, 7200), tolerance = 0.001)),
  "Nextflow duration values were not converted to seconds correctly."
)
assert_true(sarek_format_duration(3723) == "1 hr 2 min", "Readable duration formatting is incorrect.")

trace <- sarek_read_nextflow_trace(trace_path)
assert_true(NROW(trace) == 3L, "The live Nextflow trace was not read.")
estimate <- sarek_task_time_estimate(
  trace,
  "NFCORE_SAREK:SAREK:ALIGN (sample_c)",
  now = as.POSIXct("2026-08-20 12:25:00", tz = "UTC")
)
assert_true(isTRUE(estimate$available), "A supported current-task estimate was not produced.")
assert_true(estimate$expected_seconds == 720, "The current-task estimate did not use the median completed duration.")
assert_true(estimate$remaining_seconds == 420, "The estimated current-task time remaining is incorrect.")

activity <- sarek_run_activity(
  list(run_dir = run_dir, nextflow_log = log_path, trace_path = trace_path),
  state = "RUNNING"
)
assert_true(activity$completed == 2L, "Completed live task count is incorrect.")
assert_true(activity$running == 1L, "Running live task count is incorrect.")
assert_true(activity$failed == 0L, "Failed live task count is incorrect.")
assert_true(identical(activity$current_labels, "ALIGN (sample_c)"), "The active Sarek process was not identified.")
assert_true(length(activity$recent_events) == 3L, "Recent Nextflow submissions were not captured.")

log_only_activity <- sarek_run_activity(
  list(run_dir = run_dir, nextflow_log = log_path, trace_path = file.path(log_dir, "missing-trace.tsv")),
  state = "RUNNING"
)
assert_true(
  identical(log_only_activity$current_labels, "ALIGN (sample_c)"),
  "An older run without a trace did not fall back to its most recently submitted process."
)

completed_activity <- sarek_run_activity(
  list(run_dir = run_dir, nextflow_log = log_path, trace_path = trace_path),
  state = "COMPLETED"
)
assert_true(!length(completed_activity$current_steps), "A completed run still reports an active task.")

manifest <- list(
  manifest_id = "live_run",
  storage = list(results_root = test_root, work_root = file.path(test_root, "work"))
)
paths <- sarek_submission_paths(manifest)
launch <- paste(sarek_submission_launch_script(
  paths,
  runtime = list(
    launcher = "/opt/nextflow",
    config = "/opt/sarek.config",
    nxf_home = "/opt/nxf-home",
    singularity_cache = "/opt/singularity"
  ),
  pipeline = "nf-core/sarek",
  pipeline_version = "3.9.0",
  nextflow_version = "25.10.2"
), collapse = "\n")
assert_true(grepl("-with-trace", launch, fixed = TRUE), "New Sarek launches do not request a live trace.")
assert_true(grepl(paths$trace_path, launch, fixed = TRUE), "The live trace is not written inside the private run log directory.")

ui_lines <- readLines(file.path(repo_root, "R", "sarek_manifest_shiny.R"), warn = FALSE)
assert_true(any(grepl("Live Nextflow activity", ui_lines, fixed = TRUE)), "The run-status UI does not show live Nextflow activity.")
assert_true(any(grepl("Recent Nextflow events", ui_lines, fixed = TRUE)), "The run-status UI does not expose recent Nextflow events.")
assert_true(any(grepl("not a whole-workflow ETA", ui_lines, fixed = TRUE)), "The UI does not qualify the task-duration estimate.")

cat("Sarek live activity tests: PASS\n")
