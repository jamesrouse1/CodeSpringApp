args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tests/test_sarek_submission.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

if (!requireNamespace("jsonlite", quietly = TRUE)) stop("The jsonlite package is required for this test.")

source(file.path(repo_root, "R", "sarek_manifest.R"), local = FALSE)
source(file.path(repo_root, "R", "sarek_nextflow_input.R"), local = FALSE)
source(file.path(repo_root, "R", "sarek_submission.R"), local = FALSE)

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

expect_error <- function(expression, pattern, message) {
  observed <- tryCatch({ force(expression); "" }, error = function(error) conditionMessage(error))
  if (!nzchar(observed) || !grepl(pattern, observed, ignore.case = TRUE, perl = TRUE)) {
    stop(message, " Observed: ", observed, call. = FALSE)
  }
}

app_lines <- readLines(file.path(repo_root, "app.R"), warn = FALSE)
assert_true(
  any(grepl("CSL_SAREK_CONTROLLER_TIME", app_lines, fixed = TRUE)) &&
    any(grepl("controller_time = SAREK_CONTROLLER_TIME", app_lines, fixed = TRUE)),
  "app.R does not pass its administrator-configured Sarek controller time to submission."
)

make_file <- function(path, executable = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.create(path)) stop("Could not create test file: ", path)
  if (isTRUE(executable)) Sys.chmod(path, mode = "0700")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

test_root <- tempfile("sarek-submission-test-")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

input_root <- file.path(test_root, "inputs")
results_root <- file.path(test_root, "results")
work_root <- file.path(test_root, "work")
fastq_r1 <- make_file(file.path(input_root, "sample_L001_R1.fastq.gz"))
fastq_r2 <- make_file(file.path(input_root, "sample_L001_R2.fastq.gz"))

manifest <- list(
  schema_version = SAREK_MANIFEST_SCHEMA_VERSION,
  manifest_id = "submission_test",
  status = "confirmed",
  created_at = "2026-08-18T00:00:00Z",
  created_by = "test_user",
  reference = list(species = "human", assembly = "GRCh38", sarek_genome = "GATK.GRCh38"),
  assay = list(type = "WGS", intervals = NULL),
  analysis = list(mode = "germline", preset = "core"),
  storage = list(results_root = results_root, work_root = work_root),
  patients = list(list(
    patient_id = "patient1",
    relationships = list(),
    samples = list(list(
      sample_id = "sample1",
      role = "germline",
      input_format = "fastq",
      processing_state = "unmapped",
      files = list(
        list(path = fastq_r1, lane = "L001", read = 1L),
        list(path = fastq_r2, lane = "L001", read = 2L)
      )
    ))
  )),
  provenance = list(source_paths = input_root, inspector_version = SAREK_INSPECTOR_VERSION)
)
nextflow_input <- sarek_build_nextflow_input(manifest)

assert_true(sarek_submission_tools("germline") == "haplotypecaller", "Germline core tools are incorrect.")
assert_true(
  sarek_submission_tools("tumor_only") == "mutect2,controlfreec",
  "Tumor-only core tools must include Mutect2 and Control-FREEC."
)
assert_true(
  sarek_submission_tools("matched_tumor_normal") == "mutect2,manta,ascat",
  "Matched core tools must include Mutect2, Manta, and ASCAT."
)
assert_true(sarek_submission_tools("annotation_only") == "vep", "Annotation-only core tools are incorrect.")

launcher <- make_file(file.path(test_root, "backend", "bin", "nextflow-sarek"), executable = TRUE)
config <- make_file(file.path(test_root, "backend", "pipelines", "sarek", "conf", "cshl_sarek.config"))
sbatch <- make_file(file.path(test_root, "bin", "sbatch"), executable = TRUE)
nxf_home <- file.path(test_root, "backend", "runtime", "nextflow", "sarek")
singularity_cache <- file.path(test_root, "backend", "cache", "singularity")

captured <- new.env(parent = emptyenv())
result <- sarek_submit_run(
  manifest = manifest,
  nextflow_input = nextflow_input,
  launcher = launcher,
  config = config,
  nxf_home = nxf_home,
  singularity_cache = singularity_cache,
  queue = "cpuq",
  nextflow_version = "25.10.2",
  sbatch = sbatch,
  submitter = function(command, args) {
    captured$command <- command
    captured$args <- args
    "424242"
  }
)

assert_true(
  identical(result$status, "submitted"),
  "Submission result did not report submitted status."
)

assert_true(
  identical(result$job_id, ""),
  "Detached Sarek submission should not return a Slurm controller job ID."
)

assert_true(
  identical(result$controller_pid, "424242"),
  "Detached controller PID was not preserved."
)

assert_true(
  identical(result$controller_mode, "detached"),
  "Sarek controller was not recorded as detached."
)

assert_true(
  identical(captured$command, "detached"),
  "Sarek submission did not use the detached-controller starter."
)

assert_true(
  length(captured$args) >= 2L &&
    identical(captured$args[[1]], "initial"),
  "Detached controller was not started in initial mode."
)

assert_true(
  file.exists(captured$args[[2]]),
  "Detached controller launch script was not created."
)

assert_true(
  nzchar(result$child_tag),
  "Sarek submission did not create a child-job tag."
)

assert_true(
  identical(result$controller_time, "not_applicable"),
  "Detached controller should not have a Slurm time limit."
)

run_config <- file.path(
  result$run_dir,
  ".codespring",
  "nextflow.config"
)

assert_true(
  file.exists(run_config),
  "Per-run Nextflow config was not created."
)

run_config_lines <- readLines(run_config, warn = FALSE)

assert_true(
  any(grepl(
    paste0("--comment=", result$child_tag),
    run_config_lines,
    fixed = TRUE
  )),
  "Per-run Nextflow config does not tag Sarek child Slurm jobs."
)

launch_script <- file.path(
  result$run_dir,
  ".codespring",
  "launch.sh"
)

launch_lines <- readLines(launch_script, warn = FALSE)

assert_true(
  any(grepl("CSL_SAREK_RESUME", launch_lines, fixed = TRUE)),
  "Generated Sarek launcher does not support resume mode."
)

assert_true(
  any(grepl("controller_info_file", launch_lines, fixed = TRUE)),
  "Generated Sarek launcher does not record its real controller PID."
)
assert_true(file.exists(file.path(result$run_dir, ".codespring", "manifest.json")), "Internal manifest was not recorded.")
assert_true(file.exists(file.path(result$run_dir, ".codespring", "input", "samplesheet.csv")), "Internal samplesheet was not recorded.")
assert_true(file.exists(file.path(result$run_dir, ".codespring", "params.json")), "Nextflow parameter file was not recorded.")
assert_true(file.exists(file.path(result$run_dir, ".codespring", "submission.tsv")), "Submission record was not written.")

record <- sarek_read_key_value_file(file.path(result$run_dir, ".codespring", "submission.tsv"))

assert_true(
  identical(sarek_text(record[["job_id"]]), ""),
  "Detached submission record should not contain a Slurm controller job ID."
)

assert_true(
  record[["controller_pid"]] == "424242",
  "Submission record did not preserve the detached controller PID."
)

assert_true(
  record[["controller_mode"]] == "detached",
  "Submission record does not identify the controller as detached."
)

assert_true(
  nzchar(record[["child_tag"]]),
  "Submission record does not preserve the Sarek child-job tag."
)

assert_true(
  record[["tools"]] == "haplotypecaller",
  "Submission record does not preserve the selected tools."
)

assert_true(
  record[["controller_time"]] == "not_applicable",
  "Detached controller should not preserve a Slurm controller wall time."
)

catalog <- sarek_submission_catalog(results_root)

assert_true(
  NROW(catalog) == 1L &&
    catalog$run_id[[1]] == "submission_test",
  "Submitted run was not found in run history."
)

assert_true(
  identical(sarek_text(catalog$job_id[[1]]), ""),
  "Run history should not contain a Slurm controller job ID."
)

assert_true(
  catalog$controller_pid[[1]] == "424242",
  "Run history lost the detached controller PID."
)

assert_true(
  catalog$child_tag[[1]] == record[["child_tag"]],
  "Run history lost the Sarek child-job tag."
)

dir.create(file.path(result$output_dir, "variant_calling"), recursive = TRUE)
writeLines("<html><body>MultiQC</body></html>", file.path(result$output_dir, "multiqc_report.html"))
writeLines("##fileformat=VCFv4.2", file.path(result$output_dir, "variant_calling", "sample.vcf.gz"))
result_files <- sarek_result_file_catalog(result$output_dir)
assert_true(NROW(result_files) == 2L, "The result catalog did not find the public run outputs.")
assert_true(all(c("Report", "Variants") %in% result_files$type), "The result catalog did not classify reports and variants.")
assert_true(!any(grepl(".codespring", result_files$path, fixed = TRUE)), "Private run metadata leaked into the result catalog.")
limited_result_files <- sarek_result_file_catalog(result$output_dir, max_files = 1L)
assert_true(NROW(limited_result_files) == 1L && isTRUE(attr(limited_result_files, "truncated")), "The result catalog file limit is not enforced.")

params <- jsonlite::read_json(file.path(result$run_dir, ".codespring", "params.json"), simplifyVector = TRUE)
assert_true(params$step == "mapping", "Submission parameters contain the wrong starting step.")
assert_true(params$tools == "haplotypecaller", "Submission parameters contain the wrong tool preset.")
assert_true(params$genome == "GATK.GRCh38", "Submission parameters contain the wrong reference key.")
assert_true(params$outdir == result$output_dir, "Submission parameters contain the wrong output directory.")

launch_lines <- readLines(file.path(result$run_dir, ".codespring", "launch.sh"), warn = FALSE)
launch_text <- paste(launch_lines, collapse = "\n")
assert_true(grepl("nf-core/sarek", launch_text, fixed = TRUE), "Launch script does not pin the Sarek project.")
assert_true(grepl("3.9.0", launch_text, fixed = TRUE), "Launch script does not pin Sarek 3.9.0.")
assert_true(
  grepl(run_config, launch_text, fixed = TRUE),
  "Launch script does not use the generated per-run Sarek config."
)

run_config_text <- paste(
  readLines(run_config, warn = FALSE),
  collapse = "\n"
)

assert_true(
  grepl(config, run_config_text, fixed = TRUE),
  "Per-run Sarek config does not include the CSHL cluster configuration."
)

assert_true(
  grepl(
    paste0("--comment=", result$child_tag),
    run_config_text,
    fixed = TRUE
  ),
  "Per-run Sarek config does not tag child Slurm jobs."
)
assert_true(grepl("-profile singularity", launch_text, fixed = TRUE), "Launch script does not select Singularity.")
assert_true(grepl("-params-file", launch_text, fixed = TRUE), "Launch script does not use the generated parameter file.")
assert_true(grepl(result$work_dir, launch_text, fixed = TRUE), "Launch script does not use user work storage.")
assert_true(grepl("runtime_status.tsv", launch_text, fixed = TRUE), "Launch script does not record durable runtime status.")
assert_true(grepl("trap finish_controller EXIT", launch_text, fixed = TRUE), "Launch script does not record controller completion or failure.")

runtime_status <- file.path(result$run_dir, ".codespring", "runtime_status.tsv")
writeLines(
  c(
    "field\tvalue",
    "state\tCOMPLETED",
    "started_at\t2026-08-18T12:00:00Z",
    "ended_at\t2026-08-18T13:00:00Z",
    "exit_code\t0"
  ),
  runtime_status
)
catalog <- sarek_submission_catalog(results_root)
status <- sarek_run_status(as.list(catalog[1, , drop = FALSE]))
assert_true(status$state == "COMPLETED" && status$source == "runtime", "Durable completed status was not preferred over scheduler lookup.")
progress <- sarek_run_progress(status$state)
assert_true(progress$percent == 100L && progress$kind == "success", "Completed run progress is incorrect.")

queue_status <- sarek_query_slurm_job(
  "424242",
  runner = function(command, args) {
    if (grepl("squeue$", command)) return("RUNNING|00:12:34")
    character(0)
  }
)
assert_true(queue_status$state == "RUNNING" && queue_status$elapsed == "00:12:34", "Live Slurm status was not parsed.")

assert_true(
  sarek_validate_slurm_time("12:34:56") == "12:34:56" &&
    sarek_validate_slurm_time("1-02:03:04") == "1-02:03:04",
  "Valid Slurm controller time formats were rejected."
)
expect_error(
  sarek_validate_slurm_time("7 days"),
  "HH:MM:SS or D-HH:MM:SS",
  "An invalid Slurm controller time format was accepted."
)

expect_error(
  sarek_submit_run(
    manifest, nextflow_input, launcher, config, nxf_home, singularity_cache,
    queue = "cpuq", sbatch = sbatch, submitter = function(command, args) "424243"
  ),
  "already exists",
  "Duplicate run names should be rejected before Slurm submission."
)

wes_manifest <- manifest
wes_manifest$manifest_id <- "wes_without_intervals"
wes_manifest$assay$type <- "WES"
expect_error(
  sarek_build_submission_bundle(
    wes_manifest,
    nextflow_input,
    sarek_submission_validate_runtime(launcher, config, nxf_home, singularity_cache, sbatch)
  ),
  "require an intervals field",
  "WES submission without target intervals should be blocked."
)

cat("Sarek submission tests: PASS\n")
