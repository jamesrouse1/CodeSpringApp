args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tests/test_sarek_manifest_shiny.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

if (!requireNamespace("shiny", quietly = TRUE)) stop("The shiny package is required for this test.")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("The jsonlite package is required for this test.")

source(file.path(repo_root, "R", "sarek_manifest.R"), local = FALSE)
source(file.path(repo_root, "R", "sarek_bam_inspector.R"), local = FALSE)
source(file.path(repo_root, "R", "sarek_nextflow_input.R"), local = FALSE)
source(file.path(repo_root, "R", "sarek_submission.R"), local = FALSE)
source(file.path(repo_root, "R", "sarek_manifest_shiny.R"), local = FALSE)

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

app_lines <- readLines(file.path(repo_root, "app.R"), warn = FALSE)
nextflow_source_line <- grep("source(SAREK_NEXTFLOW_INPUT_HELPERS, local = FALSE)", app_lines, fixed = TRUE)
bam_inspector_source_line <- grep("source(SAREK_BAM_INSPECTOR_HELPERS, local = FALSE)", app_lines, fixed = TRUE)
submission_source_line <- grep("source(SAREK_SUBMISSION_HELPERS, local = FALSE)", app_lines, fixed = TRUE)
shiny_source_line <- grep("source(SAREK_MANIFEST_SHINY, local = FALSE)", app_lines, fixed = TRUE)
assert_true(length(nextflow_source_line) == 1L, "app.R does not source the Sarek Nextflow input helpers.")
assert_true(
  length(bam_inspector_source_line) == 1L && length(submission_source_line) == 1L && length(shiny_source_line) == 1L &&
    bam_inspector_source_line[[1]] < nextflow_source_line[[1]] &&
    nextflow_source_line[[1]] < submission_source_line[[1]] && submission_source_line[[1]] < shiny_source_line[[1]],
  "app.R must source the BAM inspector, converter, and submission helpers before the Sarek Shiny module."
)
assert_true(
  any(grepl("browse_handler = function", app_lines, fixed = TRUE)) &&
    any(grepl("updateTextAreaInput(session, path_browser$target", app_lines, fixed = TRUE)),
  "app.R does not connect the Sarek browse buttons to the shared server folder browser."
)

paths <- sarek_parse_path_input(" /tmp/a.fastq.gz\n\n/tmp/b.fastq.gz\n/tmp/a.fastq.gz ")
assert_true(identical(paths, c("/tmp/a.fastq.gz", "/tmp/b.fastq.gz")), "Path input was not trimmed and deduplicated.")
assert_true(isTRUE(sarek_parse_include_value("yes")), "Truthy include value was not recognized.")
assert_true(identical(sarek_parse_include_value("0"), FALSE), "False include value was not recognized.")
assert_true(is.na(sarek_parse_include_value("maybe")), "Ambiguous include value should remain invalid.")

widths <- sarek_manifest_column_widths()
assert_true(all(c("patient_id", "sample_id", "path", "warning") %in% names(widths)), "Manifest column widths are incomplete.")
assert_true(identical(widths[["path"]], "360px"), "The input path column should remain readable.")
assert_true(sarek_manifest_status_kind("Enter input paths.") == "info", "Instruction status should use the information style.")
assert_true(sarek_manifest_status_kind("ACTION REQUIRED: Review fields.") == "review", "Review status should use the attention style.")
assert_true(sarek_manifest_status_kind("ERROR: Missing path.") == "error", "Error status should use the error style.")
assert_true(sarek_manifest_validation_kind("No manifest has been validated yet.") == "info", "Pending validation should use the information style.")
assert_true(sarek_manifest_validation_kind("VALID\n- Ready") == "success", "Valid output should use the success style.")
assert_true(sarek_manifest_validation_kind("VALID\n- Ready\n\nWARNINGS\n- Review") == "warning", "Warnings should override the success style.")
assert_true(sarek_manifest_validation_kind("ERRORS\n- Fix this") == "error", "Validation errors should use the error style.")

test_root <- tempfile("sarek-shiny-test-")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

fastq_r1 <- file.path(test_root, "patient_T_L001_R1.fastq.gz")
fastq_r2 <- file.path(test_root, "patient_T_L001_R2.fastq.gz")
invisible(file.create(fastq_r1, fastq_r2))

table <- sarek_build_discovery_table(test_root, allowed_roots = test_root)
assert_true(NROW(table) == 2L, "Synthetic paired FASTQs were not discovered.")
summary <- sarek_confirmation_summary(table)
assert_true(summary$Value[summary$Measure == "Files"] == "2", "Discovery summary file count is incorrect.")
assert_true(summary$Value[summary$Measure == "Samples"] == "1", "Discovery summary sample count is incorrect.")

sample_review <- sarek_sample_review_table(table)
assert_true(NROW(sample_review) == 1L, "File rows were not condensed into one sample-level review row.")
assert_true(sample_review$FASTQ_pairing[[1]] == "Complete: 1 lane", "Complete paired FASTQs were not recognized.")
sample_display <- sarek_sample_review_display_table(table)
assert_true(
  identical(
    names(sample_display),
    c(
      "Include", "Patient ID", "Sex chromosomes", "Sample ID", "Role", "Matched normal", "Input format",
      "Processing state", "Index", "BAM inspection", "Files", "FASTQ pairing"
    )
  ),
  "The one-sample display table does not have twelve stable frontend columns."
)
assert_true(NROW(sample_display) == 1L && NCOL(sample_display) == 12L, "The one-sample display table collapsed during rendering.")
assert_true(startsWith(sarek_fastq_pairing_status(table[1, , drop = FALSE]), "Needs attention:"), "A missing FASTQ mate was not surfaced before validation.")
assert_true(sarek_recommend_analysis_mode(table) == "tumor_only", "Tumor FASTQs did not recommend tumor-only mode.")

bam <- file.path(test_root, "patient_T.recalibrated.bam")
bai <- paste0(bam, ".bai")
invisible(file.create(bam, bai))
bam_table <- sarek_build_discovery_table(bam, allowed_roots = test_root)
bam_display <- sarek_sample_review_display_table(bam_table)
assert_true(
  identical(bam_display$Index[[1]], paste0("Detected: ", basename(bai))),
  "The sample-facing table does not expose the detected BAM index."
)
assert_true(bam_display$`BAM inspection`[[1]] == "Not run", "An uninspected BAM has the wrong frontend status.")
bam_validation <- sarek_validate_confirmation_table(bam_table)
assert_true(
  !isTRUE(bam_validation$valid) && any(grepl("processing state", bam_validation$errors, ignore.case = TRUE)),
  "An unknown BAM processing state was not blocked during validation."
)

sample_key <- sarek_sample_key(table$patient_id[[1]], table$sample_id[[1]])
sample_updated <- sarek_apply_sample_update(
  table,
  sample_key = sample_key,
  include = TRUE,
  patient_id = "patient",
  sample_id = "patient_tumor",
  role = "tumor",
  matched_normal_id = "",
  processing_state = "unmapped"
)
assert_true(all(sample_updated$sample_id == "patient_tumor"), "A sample-level edit was not applied to every associated file.")

numeric_r1 <- file.path(test_root, "test_1.fastq.gz")
numeric_r2 <- file.path(test_root, "test_2.fastq.gz")
invisible(file.create(numeric_r1, numeric_r2))
numeric_pair <- sarek_build_discovery_table(
  c(numeric_r1, numeric_r2),
  allowed_roots = test_root
)
assert_true(
  NROW(numeric_pair) == 2L && length(unique(numeric_pair$sample_id)) == 1L &&
    unique(numeric_pair$sample_id) == "test",
  "Numeric _1/_2 FASTQs were not grouped into one sample."
)
assert_true(
  identical(sort(numeric_pair$read), c(1L, 2L)),
  "Numeric _1/_2 FASTQs were not detected as read 1 and read 2."
)
assert_true(
  sarek_fastq_pairing_status(numeric_pair) == "Complete: 1 lane",
  "Numeric _1/_2 FASTQs were not recognized as a complete pair."
)

file_updated <- sarek_apply_file_pairing_update(
  table,
  path = table$path[[1]],
  lane = "L999",
  read = 2L
)
assert_true(file_updated$lane[[1]] == "L999" && file_updated$read[[1]] == 2L, "A FASTQ lane/read correction was not applied.")

normal_r1 <- file.path(test_root, "patient_N_L001_R1.fastq.gz")
normal_r2 <- file.path(test_root, "patient_N_L001_R2.fastq.gz")
invisible(file.create(normal_r1, normal_r2))
matched_table <- sarek_build_discovery_table(
  c(fastq_r1, fastq_r2, normal_r1, normal_r2),
  allowed_roots = test_root
)
assert_true(sarek_recommend_analysis_mode(matched_table) == "matched_tumor_normal", "Tumor-normal roles did not recommend matched mode.")

patient_col <- match("patient_id", names(table)) - 1L
edited <- sarek_apply_confirmation_edit(table, list(row = 1L, col = patient_col, value = "patient_renamed"))
assert_true(edited$patient_id[[1]] == "patient_renamed", "Editable confirmation field was not updated.")

path_col <- match("path", names(table)) - 1L
readonly <- sarek_apply_confirmation_edit(table, list(row = 1L, col = path_col, value = "/tmp/replaced.fastq.gz"))
assert_true(identical(readonly$path, table$path), "Read-only file path was unexpectedly editable.")

include_col <- match("include", names(table)) - 1L
invalid_include <- sarek_apply_confirmation_edit(table, list(row = 1L, col = include_col, value = "maybe"))
assert_true(is.na(invalid_include$include[[1]]), "Ambiguous include edit should remain invalid for validation.")
include_validation <- sarek_validate_confirmation_table(invalid_include)
assert_true(!include_validation$valid, "Invalid include edit unexpectedly passed validation.")

ui_text <- paste(as.character(sarek_manifest_ui("sarek_test")), collapse = "\n")
assert_true(grepl("sarek_test-paths", ui_text, fixed = TRUE), "Shiny module inputs were not namespaced.")
assert_true(grepl("review and correct samples", ui_text, ignore.case = TRUE), "UI does not clearly identify the sample-review step.")
assert_true(grepl("sarek-confirmation-table", ui_text, fixed = TRUE), "Manifest table is missing its horizontal-scroll container.")
assert_true(grepl("sarek-status-banner", ui_text, fixed = TRUE), "Manifest status is missing its visible banner style.")
assert_true(grepl("sarek-review-checklist", ui_text, fixed = TRUE), "The pre-validation checklist is missing.")
assert_true(grepl("sarek_test-sample_editor", ui_text, fixed = TRUE), "The sample-level editor output is missing.")
assert_true(grepl("sarek_test-bam_inspection", ui_text, fixed = TRUE), "The BAM inspection panel output is missing.")
assert_true(grepl("sarek_test-file_editor", ui_text, fixed = TRUE), "The FASTQ lane/read correction output is missing.")
assert_true(grepl("Field definitions and rules", ui_text, fixed = TRUE), "The manifest field guide is missing.")
assert_true(grepl("Results root", ui_text, fixed = TRUE) && grepl("Work root", ui_text, fixed = TRUE), "Storage fields are not explained.")
assert_true(grepl("click any sample row once", ui_text, ignore.case = TRUE), "The table does not clearly explain how to edit a sample.")
assert_true(grepl("sarek_test-validation", ui_text, fixed = TRUE), "The visible validation panel output is missing.")
assert_true(grepl("sarek_test-nextflow_input", ui_text, fixed = TRUE), "The Nextflow input panel output is missing.")
assert_true(grepl("sarek_test-run_review", ui_text, fixed = TRUE), "The front-end run review is missing.")
assert_true(grepl("sarek_test-submission_ui", ui_text, fixed = TRUE), "The submission controls are missing.")
assert_true(grepl("Previous and active Sarek runs", ui_text, fixed = TRUE), "The Sarek run-history section is missing.")
assert_true(grepl("sarek_test-selected_run_status", ui_text, fixed = TRUE), "The selected-run progress output is missing.")
assert_true(grepl("sarek_test-sarek_sections", ui_text, fixed = TRUE), "The top-level Sarek tab navigation is missing.")
assert_true(
  all(vapply(c("Prepare job", "Run status", "Results"), grepl, logical(1), x = ui_text, fixed = TRUE)),
  "The Prepare job, Run status, and Results tabs are not all visible."
)
assert_true(grepl("sarek_test-results_run_selector", ui_text, fixed = TRUE), "The results run selector is missing.")
assert_true(grepl("sarek_test-results_file_table", ui_text, fixed = TRUE), "The searchable result file table is missing.")
assert_true(grepl("sarek_test-browse_input_folder", ui_text, fixed = TRUE), "Input folder browsing is missing.")
assert_true(grepl("sarek_test-browse_results_root", ui_text, fixed = TRUE), "Results-root browsing is missing.")
assert_true(grepl("sarek_test-browse_work_root", ui_text, fixed = TRUE), "Work-root browsing is missing.")
assert_true(grepl("sarek-file-details", ui_text, fixed = TRUE), "Selected-file details do not have an isolated layout container.")
assert_true(grepl("Show input and index paths", ui_text, fixed = TRUE), "The selected-file section does not advertise index visibility.")
assert_true(!grepl("JSON", ui_text, fixed = TRUE), "Developer-facing JSON remains visible in the end-user UI.")
assert_true(!grepl("download_manifest", ui_text, fixed = TRUE), "The end-user JSON download remains visible.")
nextflow_position <- regexpr("sarek_test-nextflow_input", ui_text, fixed = TRUE)[[1]]
review_position <- regexpr("sarek_test-run_review", ui_text, fixed = TRUE)[[1]]
submission_position <- regexpr("sarek_test-submission_ui", ui_text, fixed = TRUE)[[1]]
assert_true(
  nextflow_position > 0L && review_position > nextflow_position && submission_position > review_position,
  "Run review and submission controls must appear after the generated input."
)

results_root <- file.path(test_root, "results")
work_root <- file.path(test_root, "work")
older_internal <- file.path(results_root, "older_run", ".codespring")
dir.create(older_internal, recursive = TRUE)
writeLines(
  c(
    "field\tvalue", "status\tsubmitted", "job_id\t111111",
    "submitted_at\t2026-08-17T12:00:00Z",
    paste0("run_dir\t", file.path(results_root, "older_run")),
    paste0("output_dir\t", file.path(results_root, "older_run", "results")),
    paste0("work_dir\t", file.path(work_root, "older_run")),
    "step\tmapping", "tools\tmutect2"
  ),
  file.path(older_internal, "submission.tsv")
)
writeLines(
  c("field\tvalue", "state\tCOMPLETED", "exit_code\t0"),
  file.path(older_internal, "runtime_status.tsv")
)
older_results <- file.path(results_root, "older_run", "results")
dir.create(file.path(older_results, "variant_calling"), recursive = TRUE)
writeLines("<html><body>MultiQC</body></html>", file.path(older_results, "multiqc_report.html"))
writeLines("##fileformat=VCFv4.2", file.path(older_results, "variant_calling", "older_run.vcf"))
submission_capture <- new.env(parent = emptyenv())
submission_capture$count <- 0L
browse_capture <- new.env(parent = emptyenv())
browse_capture$count <- 0L

shiny::testServer(
  sarek_manifest_server,
  args = list(
    default_results_root = results_root,
    default_work_root = work_root,
    created_by = "test_user",
    allowed_input_roots = test_root,
    allowed_results_roots = test_root,
    allowed_work_roots = test_root,
    max_files = 10L,
    run_refresh_ms = 60000L,
    bam_inspector = function(path, index, samtools) {
      sarek_bam_inspection_result(
        path = path,
        index = index,
        status = "passed",
        summary = "Inspection passed; suggested processing state: recalibrated (high confidence).",
        quickcheck_ok = TRUE,
        index_ok = TRUE,
        sample_ids = "TUMOR_01",
        sort_order = "coordinate",
        grch38_compatible = TRUE,
        reference_summary = "chr1 length matches GRCh38 and the contig naming matches GATK.GRCh38.",
        recommendation = "recalibrated",
        confidence = "high",
        evidence = c("samtools quickcheck passed.", "Applied BQSR was detected.")
      )
    },
    browse_handler = function(target, mode, current, input_type, append) {
      browse_capture$count <- browse_capture$count + 1L
      browse_capture$target <- target
      browse_capture$mode <- mode
      browse_capture$input_type <- input_type
      browse_capture$append <- append
      invisible(TRUE)
    },
    submit_handler = function(manifest, nextflow_input) {
      submission_capture$count <- submission_capture$count + 1L
      submission_capture$manifest <- manifest
      submission_capture$nextflow_input <- nextflow_input
      list(
        status = "submitted",
        job_id = "424242",
        run_dir = file.path(results_root, manifest$manifest_id),
        output_dir = file.path(results_root, manifest$manifest_id, "results"),
        work_dir = file.path(work_root, manifest$manifest_id),
        tools = sarek_submission_tools(manifest$analysis$mode, manifest$analysis$preset),
        step = nextflow_input$step
      )
    }
  ),
  {
    session$setInputs(
      paths = paste(fastq_r1, fastq_r2, sep = "\n"),
      recursive = FALSE,
      manifest_id = "tumor_test",
      assay_type = "WGS",
      analysis_mode = "tumor_only",
      preset = "core",
      results_root = results_root,
      work_root = work_root,
      species = "human",
      assembly = "GRCh38",
      sarek_genome = "GATK.GRCh38"
    )
    session$setInputs(discover = 1L)
    session$flushReact()
    assert_true(NROW(confirmation_state()) == 2L, "Module discovery did not populate the confirmation table.")
    assert_true(is.null(confirmed_manifest()), "Discovery should not automatically confirm a manifest.")
    editor_html <- paste(as.character(output$sample_editor), collapse = "\n")
    assert_true(grepl("edit_include", editor_html, fixed = TRUE), "The selected-sample editor does not expose an include checkbox.")
    sample_table_html <- paste(as.character(output$sample_review_table), collapse = "\n")
    assert_true(!grepl("names.*attribute", sample_table_html, ignore.case = TRUE), "The one-sample frontend still throws the column-name rendering error.")

    history_html <- paste(as.character(output$run_history_selector), collapse = "\n")
    assert_true(grepl("older_run", history_html, fixed = TRUE), "An older submitted Sarek run is not available in the run selector.")
    session$setInputs(selected_run = "older_run")
    session$flushReact()
    progress_html <- paste(as.character(output$selected_run_status), collapse = "\n")
    assert_true(grepl("Completed", progress_html, fixed = TRUE), "Completed run status is not shown in the history panel.")
    assert_true(grepl("sarek-progress-bar", progress_html, fixed = TRUE), "The Sarek progress bar is missing.")
    results_selector_html <- paste(as.character(output$results_run_selector), collapse = "\n")
    assert_true(grepl("older_run", results_selector_html, fixed = TRUE), "An older run is not available in the Results tab.")
    session$setInputs(results_run = "older_run")
    session$flushReact()
    results_summary_html <- paste(as.character(output$selected_results_summary), collapse = "\n")
    assert_true(grepl("2 readable result files", results_summary_html, fixed = TRUE), "The Results tab does not summarize readable outputs.")
    assert_true(grepl(older_results, results_summary_html, fixed = TRUE), "The selected run results location is not visible.")
    assert_true(NROW(results_file_catalog()) == 2L, "The Results tab did not catalog the selected run outputs.")

    session$setInputs(browse_input_folder = 1L)
    session$flushReact()
    assert_true(browse_capture$count == 1L, "Input folder browse action did not call the shared server browser.")
    assert_true(identical(browse_capture$mode, "dir") && identical(browse_capture$input_type, "textarea") && isTRUE(browse_capture$append), "Input browse action did not request an appended server folder.")

    session$setInputs(sample_review_table_rows_selected = 1L)
    session$flushReact()
    expected_key <- sarek_sample_key(confirmation_state()$patient_id[[1]], confirmation_state()$sample_id[[1]])
    assert_true(selected_sample_state() == expected_key, "A single table-row selection did not select the sample for editing.")

    reviewed_sample <- confirmation_state()[1, , drop = FALSE]
    session$setInputs(
      edit_include = TRUE,
      edit_patient_id = reviewed_sample$patient_id[[1]],
      edit_sample_id = reviewed_sample$sample_id[[1]],
      edit_sex = "XX",
      edit_role = "tumor",
      edit_matched_normal_id = "",
      edit_processing_state = "unmapped"
    )
    session$setInputs(apply_sample = 1L)
    session$flushReact()
    assert_true(
      all(confirmation_state()$sex == "XX"),
      "The corroborated patient sex chromosomes were not applied to the FASTQ sample."
    )

    session$setInputs(confirm = 1L)
    session$flushReact()
    manifest <- confirmed_manifest()
    assert_true(!is.null(manifest), "Valid reviewed inputs did not produce a confirmed manifest.")
    validation_html <- paste(as.character(output$validation), collapse = "\n")
    assert_true(grepl("sarek-validation-panel-success", validation_html, fixed = TRUE), "Successful validation did not produce a green success panel.")
    nextflow_state <- nextflow_input_state()
    assert_true(isTRUE(nextflow_state$valid), "A valid confirmed manifest did not generate a Nextflow input.")
    assert_true(isTRUE(submission_readiness()$valid), "A valid confirmed WGS run was not marked ready for submission.")
    assert_true(nextflow_state$input$step == "mapping", "FASTQ input did not select the Sarek mapping step.")
    assert_true(NROW(nextflow_state$input$samplesheet) == 1L, "One FASTQ pair did not produce one Sarek samplesheet row.")
    nextflow_html <- paste(as.character(output$nextflow_input), collapse = "\n")
    assert_true(grepl("Nextflow input ready", nextflow_html, fixed = TRUE), "The Nextflow readiness panel was not rendered.")
    assert_true(grepl("Starting step: mapping", nextflow_html, fixed = TRUE), "The resolved Sarek starting step is not visible.")
    review_html <- paste(as.character(output$run_review), collapse = "\n")
    assert_true(grepl("tumor_test", review_html, fixed = TRUE), "The confirmed run name is not visible in the front-end review.")
    assert_true(grepl("mutect2", review_html, fixed = TRUE), "The mode-aware tool preset is not visible in the run review.")
    assert_true(grepl("controlfreec", review_html, fixed = TRUE), "The tumor-only CNV caller is not visible in the run review.")
    samples_html <- paste(as.character(output$confirmed_samples), collapse = "\n")
    assert_true(grepl("patient_T", samples_html, fixed = TRUE), "The confirmed sample table is not visible.")
    submission_html <- paste(as.character(output$submission_ui), collapse = "\n")
    assert_true(grepl("validation_reviewed", submission_html, fixed = TRUE), "Review acknowledgement is missing before submission.")
    assert_true(!grepl("request_submit", submission_html, fixed = TRUE), "Submit was enabled before review acknowledgement.")
    assert_true(!grepl("download_nextflow_input", submission_html, fixed = TRUE), "Advanced samplesheet download was enabled before review acknowledgement.")
    session$setInputs(validation_reviewed = TRUE)
    session$flushReact()
    submission_html <- paste(as.character(output$submission_ui), collapse = "\n")
    assert_true(grepl("request_submit", submission_html, fixed = TRUE), "Submit was not enabled after review acknowledgement.")
    assert_true(grepl("download_nextflow_input", submission_html, fixed = TRUE), "Advanced samplesheet download was not enabled after review acknowledgement.")
    assert_true(!grepl("download_manifest", submission_html, fixed = TRUE), "Developer-facing JSON download is still present.")
    assert_true(manifest$analysis$mode == "tumor_only", "Confirmed manifest has the wrong analysis mode.")
    assert_true(manifest$manifest_id == "tumor_test", "Confirmed manifest has the wrong ID.")
    assert_true(length(manifest$patients[[1]]$samples[[1]]$files) == 2L, "Confirmed manifest lost a FASTQ mate.")

    session$setInputs(request_submit = 1L)
    session$flushReact()
    session$setInputs(confirm_submit = 1L)
    session$flushReact()
    assert_true(submission_capture$count == 1L, "The confirmed submission did not call the backend exactly once.")
    assert_true(submission_capture$nextflow_input$step == "mapping", "The backend received the wrong starting step.")
    assert_true(identical(submission_state()$result$job_id, "424242"), "The submitted Slurm job ID was not recorded.")
    submission_status_html <- paste(as.character(output$submission_status), collapse = "\n")
    assert_true(grepl("sarek-validation-panel-success", submission_status_html, fixed = TRUE), "Successful submission did not render a green status panel.")
    assert_true(grepl("424242", submission_status_html, fixed = TRUE), "The submitted Slurm job ID is not visible.")

    session$setInputs(preset = "changed_after_confirmation")
    session$flushReact()
    assert_true(is.null(confirmed_manifest()), "Changing settings should invalidate the prior confirmation.")

    session$setInputs(
      paths = bam,
      analysis_mode = "tumor_only",
      preset = "core"
    )
    session$setInputs(discover = 2L)
    session$flushReact()
    inspected_bam <- confirmation_state()
    assert_true(NROW(inspected_bam) == 1L, "BAM discovery did not produce one input row.")
    assert_true(inspected_bam$inspection_status[[1]] == "passed", "Automatic BAM inspection did not run during discovery.")
    bam_table_html <- paste(as.character(output$sample_review_table), collapse = "\n")
    assert_true(grepl(basename(bai), bam_table_html, fixed = TRUE), "The main sample table does not show the detected BAM index.")
    bam_panel_html <- paste(as.character(output$bam_inspection), collapse = "\n")
    assert_true(grepl("TUMOR_01", bam_panel_html, fixed = TRUE), "The BAM inspection panel does not show the header sample ID.")
    assert_true(grepl("recalibrated", bam_panel_html, fixed = TRUE), "The BAM inspection recommendation is not visible.")
    bam_editor_html <- paste(as.character(output$sample_editor), collapse = "\n")
    assert_true(
      grepl("not confirmed until you select Apply sample changes", bam_editor_html, fixed = TRUE),
      "The UI does not explain that the BAM recommendation requires corroboration."
    )
    session$setInputs(
      edit_include = TRUE,
      edit_patient_id = inspected_bam$patient_id[[1]],
      edit_sample_id = inspected_bam$sample_id[[1]],
      edit_role = "tumor",
      edit_matched_normal_id = "",
      edit_processing_state = "recalibrated"
    )
    session$setInputs(apply_sample = 1L)
    session$flushReact()
    assert_true(
      confirmation_state()$processing_state[[1]] == "recalibrated",
      "The corroborated BAM processing state was not applied."
    )
  }
)

cat("Sarek Shiny manifest tests: PASS\n")
