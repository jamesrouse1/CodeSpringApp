args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[[1]])
} else {
  "tests/test_sarek_detached_lifecycle.R"
}

repo_root <- normalizePath(
  file.path(dirname(script_path), ".."),
  winslash = "/",
  mustWork = TRUE
)

source(file.path(repo_root, "R", "sarek_manifest.R"), local = FALSE)
source(file.path(repo_root, "R", "sarek_nextflow_input.R"), local = FALSE)
source(file.path(repo_root, "R", "sarek_submission.R"), local = FALSE)

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

read_value <- function(path, field) {
  values <- sarek_read_key_value_file(path)
  sarek_text(values[field])
}

test_root <- tempfile("sarek-detached-lifecycle-")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

run_dir <- file.path(test_root, "run")
internal_dir <- file.path(run_dir, ".codespring")
log_dir <- file.path(internal_dir, "logs")
work_dir <- file.path(test_root, "work")
bin_dir <- file.path(test_root, "bin")

dir.create(log_dir, recursive = TRUE)
dir.create(bin_dir, recursive = TRUE)

child_tag <- "codespring_sarek_lifecycle_test"

fake_launcher <- file.path(bin_dir, "nextflow-sarek")
writeLines(
  c(
    "#!/usr/bin/env bash",
    "exit \"${FAKE_NEXTFLOW_EXIT:-0}\""
  ),
  fake_launcher
)
Sys.chmod(fake_launcher, "0700")

fake_squeue <- file.path(bin_dir, "squeue")
writeLines(
  c(
    "#!/usr/bin/env bash",
    "if [[ \"${FAKE_CHILD_MODE:-empty}\" == \"active\" ]]; then",
    paste0(
      "  printf '999001|RUNNING|0:10|8:00:00|",
      child_tag,
      "|nf-TEST\\n'"
    ),
    "fi",
    "exit 0"
  ),
  fake_squeue
)
Sys.chmod(fake_squeue, "0700")

run_config <- file.path(internal_dir, "nextflow.config")
params_path <- file.path(internal_dir, "params.json")
file.create(run_config)
file.create(params_path)

paths <- list(
  manifest_id = "lifecycle_test",
  run_dir = run_dir,
  output_dir = file.path(run_dir, "results"),
  internal_dir = internal_dir,
  log_dir = log_dir,
  work_dir = work_dir,
  params_path = params_path,
  run_config = run_config,
  launch_script = file.path(internal_dir, "launch.sh"),
  nextflow_log = file.path(log_dir, "nextflow.log"),
  trace_path = file.path(log_dir, "trace.tsv"),
  stdout = file.path(log_dir, "controller.out"),
  stderr = file.path(log_dir, "controller.err"),
  submission_record = file.path(internal_dir, "submission.tsv"),
  runtime_status = file.path(internal_dir, "runtime_status.tsv"),
  active_children = file.path(internal_dir, "active_children.tsv"),
  controller_info = file.path(internal_dir, "controller.tsv"),
  child_tag = child_tag
)

runtime <- list(
  launcher = fake_launcher,
  config = run_config,
  nxf_home = file.path(test_root, "nxf"),
  singularity_cache = file.path(test_root, "singularity")
)

writeLines(
  sarek_submission_launch_script(
    paths = paths,
    runtime = runtime,
    pipeline = "nf-core/sarek",
    pipeline_version = "test",
    nextflow_version = "25.10.2"
  ),
  paths$launch_script
)
Sys.chmod(paths$launch_script, "0700")

old_path <- Sys.getenv("PATH")
on.exit(Sys.setenv(PATH = old_path), add = TRUE)

Sys.setenv(
  PATH = paste(bin_dir, old_path, sep = .Platform$path.sep)
)

run_case <- function(exit_code, child_mode, expected_state) {
  Sys.setenv(
    FAKE_NEXTFLOW_EXIT = as.character(exit_code),
    FAKE_CHILD_MODE = child_mode,
    CSL_SAREK_RESUME = "false",
    CSL_SAREK_CONTROLLER_MODE = "initial"
  )

  status <- suppressWarnings(
    system2(
      "bash",
      paths$launch_script,
      stdout = TRUE,
      stderr = TRUE
    )
  )

  observed_exit <- attr(status, "status")
  if (is.null(observed_exit)) observed_exit <- 0L

  state <- read_value(paths$runtime_status, "state")

  assert_true(
    identical(state, expected_state),
    paste0(
      "Expected state ", expected_state,
      " but observed ", state
    )
  )

  if (exit_code == 0L) {
    assert_true(
      observed_exit == 0L,
      "Successful fake Nextflow run returned non-zero."
    )
  } else {
    assert_true(
      observed_exit == exit_code,
      "Failed fake Nextflow exit code was not preserved."
    )
  }
}

# 1. Nextflow succeeds + no active Slurm children => COMPLETED
run_case(
  exit_code = 0L,
  child_mode = "empty",
  expected_state = "COMPLETED"
)

# The actual launch shell must identify itself.
controller_pid <- read_value(paths$controller_info, "pid")
controller_host <- read_value(paths$controller_info, "host")

assert_true(
  grepl("^[0-9]+$", controller_pid),
  "controller.tsv did not contain a numeric controller PID."
)

assert_true(
  nzchar(controller_host),
  "controller.tsv did not contain the controller hostname."
)

# 2. Nextflow exits 0 but tagged Slurm child remains => INCOMPLETE
run_case(
  exit_code = 0L,
  child_mode = "active",
  expected_state = "INCOMPLETE"
)

children <- readLines(paths$active_children, warn = FALSE)

assert_true(
  length(children) == 2L &&
    grepl("999001", children[[2]], fixed = TRUE),
  "Active child snapshot was not recorded."
)

# 3. Nextflow fails => FAILED
run_case(
  exit_code = 7L,
  child_mode = "empty",
  expected_state = "FAILED"
)

assert_true(
  identical(read_value(paths$runtime_status, "exit_code"), "7"),
  "Nextflow failure exit code was not written to runtime status."
)

# 4. R-side Slurm child lookup recognizes the run tag.
children <- sarek_query_active_children(
  child_tag,
  runner = function(command, args) {
    paste0(
      "999002|RUNNING|00:12|08:00:00|",
      child_tag,
      "|nf-FAKE_TASK"
    )
  }
)

assert_true(
  NROW(children) == 1L &&
    identical(children$job_id[[1]], "999002"),
  "Tagged Slurm child lookup failed."
)

# 5. A dead detached controller that still says RUNNING is INCOMPLETE,
#    never falsely COMPLETED.
writeLines(
  c(
    "field\tvalue",
    "state\tRUNNING",
    paste0("child_tag\t", child_tag)
  ),
  paths$runtime_status
)

run_record <- list(
  run_dir = run_dir,
  runtime_status = paths$runtime_status,
  controller_pid = "999999999",
  controller_host = sarek_text(Sys.info()[["nodename"]]),
  child_tag = child_tag,
  status = "submitted"
)

status <- sarek_run_status(
  run_record,
  runner = function(command, args) character(0)
)

assert_true(
  identical(status$state, "INCOMPLETE"),
  paste0(
    "Dead RUNNING controller should be INCOMPLETE, observed ",
    status$state
  )
)

cat("Sarek detached lifecycle: PASS\n")
