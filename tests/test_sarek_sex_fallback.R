args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tests/test_sarek_sex_fallback.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

source(file.path(repo_root, "R", "sarek_manifest.R"), local = FALSE)
source(file.path(repo_root, "R", "sarek_nextflow_input.R"), local = FALSE)
source(file.path(repo_root, "R", "sarek_submission.R"), local = FALSE)

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

test_root <- tempfile("sarek-sex-fallback-")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

make_fastq_pair <- function(sample_id) {
  paths <- file.path(test_root, paste0(sample_id, c("_R1.fastq.gz", "_R2.fastq.gz")))
  file.create(paths)
  list(
    list(path = paths[[1]], lane = "L001", read = 1L),
    list(path = paths[[2]], lane = "L001", read = 2L)
  )
}

make_sample <- function(sample_id, role) {
  list(
    sample_id = sample_id,
    role = role,
    input_format = "fastq",
    processing_state = "unmapped",
    files = make_fastq_pair(sample_id)
  )
}

make_manifest <- function(mode, sex = "NA", matched = FALSE) {
  samples <- list(make_sample("tumor1", "tumor"))
  relationships <- list()
  if (isTRUE(matched)) {
    samples <- c(samples, list(make_sample("normal1", "normal")))
    relationships <- list(list(
      type = "matched_tumor_normal",
      sample_ids = c("tumor1", "normal1")
    ))
  }
  list(
    schema_version = SAREK_MANIFEST_SCHEMA_VERSION,
    manifest_id = paste0("sex_test_", mode),
    status = "confirmed",
    analysis = list(mode = mode, preset = "core"),
    patients = list(list(
      patient_id = "patient1",
      sex = sex,
      samples = samples,
      relationships = relationships
    ))
  )
}

assert_true(
  identical(unname(vapply(c("xx", "XY", "unknown"), sarek_normalize_sex, character(1))), c("XX", "XY", "NA")),
  "Sex chromosome normalization is incorrect."
)

tumor_missing <- make_manifest("tumor_only")
tumor_input <- sarek_build_nextflow_input(tumor_missing)
tumor_resolution <- sarek_submission_tool_resolution("tumor_only", "core", tumor_missing)
assert_true(all(tumor_input$samplesheet$sex == "NA"), "Missing sex was not written as NA in the Sarek samplesheet.")
assert_true(length(tumor_input$warnings) == 1L, "Tumor-only input did not warn about missing sex.")
assert_true(tumor_resolution$tools == "mutect2", "Control-FREEC was not omitted when sex was missing.")
assert_true(identical(tumor_resolution$skipped_tools, "controlfreec"), "Skipped Control-FREEC was not recorded.")
assert_true(grepl("run will continue", tumor_resolution$warnings, ignore.case = TRUE), "Caller warning does not explain that the run continues.")

tumor_with_sex <- make_manifest("tumor_only", sex = "XX")
tumor_with_sex_input <- sarek_build_nextflow_input(tumor_with_sex)
assert_true(all(tumor_with_sex_input$samplesheet$sex == "XX"), "Confirmed XX was not preserved in the Sarek samplesheet.")
assert_true(
  sarek_submission_tools("tumor_only", "core", tumor_with_sex) == "mutect2,controlfreec",
  "Control-FREEC was not retained when sex was provided."
)

matched_missing <- make_manifest("matched_tumor_normal", matched = TRUE)
matched_resolution <- sarek_submission_tool_resolution("matched_tumor_normal", "core", matched_missing)
assert_true(matched_resolution$tools == "mutect2,manta", "ASCAT was not omitted when sex was missing.")
assert_true(identical(matched_resolution$skipped_tools, "ascat"), "Skipped ASCAT was not recorded.")

matched_with_sex <- make_manifest("matched_tumor_normal", sex = "XY", matched = TRUE)
assert_true(
  sarek_submission_tools("matched_tumor_normal", "core", matched_with_sex) == "mutect2,manta,ascat",
  "ASCAT was not retained when sex was provided."
)

ui_lines <- readLines(file.path(repo_root, "R", "sarek_manifest_shiny.R"), warn = FALSE)
assert_true(any(grepl('ns("edit_sex")', ui_lines, fixed = TRUE)), "The manifest editor does not expose sex chromosomes.")
assert_true(any(grepl("Caller adjusted", ui_lines, fixed = TRUE)), "The confirmation UI does not display the caller adjustment warning.")

cat("Sarek sex and caller fallback tests: PASS\n")
