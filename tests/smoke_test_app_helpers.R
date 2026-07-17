args <- commandArgs(trailingOnly = TRUE)
repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), ".."), mustWork = TRUE)
lab_root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(file.path(repo_root, "..", "CodeSpringLab-fix"), mustWork = TRUE)
Sys.setenv(CSL_CODESPRINGLAB_ROOT = lab_root)

app_env <- new.env(parent = globalenv())
sys.source(file.path(repo_root, "app.R"), envir = app_env)

assert <- function(value, message) if (!isTRUE(value)) stop("ASSERTION FAILED: ", message, call. = FALSE)
assert(app_env$is_codespring_process_command("Rscript -e shiny::runApp('/home/user/CodeSpringWeb', port=8601)"), "CodeSpringApp process command recognized")
assert(!app_env$is_codespring_process_command("Rscript unrelated_analysis.R"), "unrelated Rscript process is not treated as CodeSpringApp")
assert(!app_env$is_codespring_process_command("Rscript -e shiny::runApp('/home/user/another_app')"), "unrelated Shiny app is not treated as CodeSpringApp")
root <- tempfile("codespring-app-smoke-")
dir.create(root, recursive = TRUE)
on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

design_path <- file.path(root, "design_matrix.txt")
design <- data.frame(
  sample = c("A1", "I1", "A2", "I2", "B1", "I3", "B2", "I4"),
  treatment = rep(c("A", "A", "B", "B"), each = 2),
  reference = rep(c("chip", "input"), 4),
  condition = rep(c("A", "A", "B", "B"), each = 2),
  replicate = rep(c(1, 1, 2, 2), each = 2),
  control_sample = c("I1", "", "I2", "", "I3", "", "I4", ""),
  filename = paste0(c("A1", "I1", "A2", "I2", "B1", "I3", "B2", "I4"), ".fastq.gz"),
  stringsAsFactors = FALSE
)
write.table(design, design_path, sep = "\t", row.names = FALSE, quote = FALSE)

chip_project <- list(
  id = "fake-chip", name = "fake-chip", analysis_key = "chip", analysis = "ChIP-seq",
  design_matrix_path = design_path, data_dir = root, results_root = dirname(root),
  fastq_dir = root, fastq_dirs = root, paired_end = FALSE, genome = "mouse"
)
assert(identical(app_env$chip_control_sample_for(chip_project, "A1"), "I1"), "explicit ChIP control resolution")
assert(nrow(app_env$chip_target_design(chip_project)) == 4L, "input rows excluded from ChIP targets")
mouse_chip_ref <- app_env$chip_reference_resources(chip_project)
human_chip_project <- chip_project
human_chip_project$genome <- "human"
human_chip_ref <- app_env$chip_reference_resources(human_chip_project)
assert(identical(mouse_chip_ref$genome_version, "mouse_gencodeM39") && grepl("mouse_gencodeM39", mouse_chip_ref$bowtie2_index), "ChIP mouse reference uses GRCm39/GENCODE M39")
assert(identical(human_chip_ref$genome_version, "human_gencode50") && grepl("human_gencode50", human_chip_ref$bowtie2_index), "ChIP human reference uses GRCh38/GENCODE v50")
assert(length(app_env$genome_reference_choices("mouse", "ChIP-seq")) == 1L, "ChIP setup offers only the current mouse reference")
assert(length(app_env$genome_reference_choices("human", "ChIP-seq")) == 1L, "ChIP setup offers only the current human reference")

duplicate_design <- data.frame(
  include = TRUE, sample = c("sample-A", "sample A"), cell_type = "", condition = c("A", "B"),
  replicate = 1:2, filename = c("a.fastq.gz", "b.fastq.gz"), status = "", stringsAsFactors = FALSE
)
atac_design_project <- chip_project
atac_design_project$analysis_key <- "atac"
atac_design_project$analysis <- "ATAC-seq"
duplicate_error <- tryCatch({
  app_env$write_design_matrix(atac_design_project, duplicate_design, c("condition", "replicate"))
  ""
}, error = conditionMessage)
assert(grepl("remain unique", duplicate_error), "filesystem-safe sample collisions rejected")

blank_design <- duplicate_design[1, , drop = FALSE]
blank_design$sample <- "sample1"
blank_design$filename <- ""
blank_error <- tryCatch({ app_env$write_design_matrix(atac_design_project, blank_design, c("condition", "replicate")); "" }, error = conditionMessage)
assert(grepl("FASTQ filename", blank_error), "blank included FASTQ filenames rejected")

valid_design <- duplicate_design[1, , drop = FALSE]
valid_design$sample <- "sample1"
saved_design <- app_env$write_design_matrix(atac_design_project, valid_design, c("condition", "replicate"))
assert(file.exists(saved_design) && file.info(saved_design)$size > 0, "design matrix saved atomically")
saved_table <- app_env$safe_read_table(saved_design)
assert(all(c("cell_type", "condition", "replicate") %in% names(saved_table)), "required ATAC metadata columns preserved")

unsafe_design <- valid_design
unsafe_design$condition <- "A\tB"
unsafe_error <- tryCatch({ app_env$write_design_matrix(atac_design_project, unsafe_design, c("condition", "replicate")); "" }, error = conditionMessage)
assert(grepl("tabs or line breaks", unsafe_error), "tab characters rejected before TSV save")

for (key in c("rna", "atac", "chip")) {
  example <- app_env$example_dataset_paths(key)
  example_design <- file.path(example$design_dir, "design_matrix.txt")
  assert(dir.exists(example$fastq_dir) && file.exists(example_design), paste(key, "bundled example paths exist"))
  table <- app_env$safe_read_table(example_design)
  assert(NROW(table) > 0 && !anyDuplicated(table$sample), paste(key, "bundled example design has unique samples"))
  reads <- trimws(unlist(strsplit(as.character(table$filename), "[;,]")))
  read_paths <- file.path(example$fastq_dir, reads[nzchar(reads)])
  assert(all(file.exists(read_paths)), paste(key, "bundled example FASTQs match the design"))
  readable <- vapply(read_paths, function(path) {
    connection <- gzfile(path, open = "rt")
    on.exit(close(connection), add = TRUE)
    length(readLines(connection, n = 4L, warn = FALSE)) == 4L
  }, logical(1))
  assert(all(readable), paste(key, "bundled example FASTQs are readable gzip data"))
}

atac_project <- chip_project
atac_project$analysis_key <- "atac"
atac_project$analysis <- "ATAC-seq"
sample_dir <- file.path(root, "macs2", "A1")
dir.create(sample_dir, recursive = TRUE)
legacy_peak <- file.path(sample_dir, "A1_peaks.narrowPeak")
run_log <- file.path(sample_dir, "A1_macs2.log")
marker <- file.path(sample_dir, "A1_macs2_complete.txt")
writeLines("chr1\t1\t2", legacy_peak)
assert(identical(app_env$atac_macs2_completion_target(atac_project, "A1"), legacy_peak), "legacy ATAC peaks remain recognized")
writeLines("Traceback (most recent call last):\nOSError: No space left on device", run_log)
assert(identical(app_env$atac_macs2_completion_target(atac_project, "A1"), marker), "new ATAC runs require a completion marker")
assert(app_env$cutrun_macs_fatal_error_signal(atac_project, data.frame(), "MACS2 Peaks", "A1"), "ATAC internal MACS2 exception detection")
writeLines("status\tcomplete", marker)
assert(identical(app_env$atac_macs2_completion_target(atac_project, "A1"), marker), "completed ATAC marker selected")
assert(identical(app_env$atac_macs2_peak_file(atac_project, "A1"), legacy_peak), "validated ATAC peak selected for DiffBind")

unlink(marker)
assert(identical(app_env$chip_macs2_peak_file(chip_project, "A1"), ""), "partial ChIP MACS2 output rejected")
writeLines("status\tcomplete", marker)
assert(identical(app_env$chip_macs2_peak_file(chip_project, "A1"), legacy_peak), "completed ChIP MACS2 peak accepted")
chip_peaks <- app_env$chip_peak_summary_table(chip_project)
assert(NROW(chip_peaks) == 4L && chip_peaks$status[chip_peaks$sample == "A1"] == "Completed", "ChIP matched-input peak summary reports completion")
alignment_dir <- file.path(root, "bowtie2", "A1")
dir.create(alignment_dir, recursive = TRUE)
writeLines(c("sample\tA1", "mapped_reads\t100", "deduplicated_reads\t80", "bigwig_normalization\tCPM"), file.path(alignment_dir, "A1_alignment_summary.txt"))
signal_file <- file.path(alignment_dir, "A1Aligned.sortedByCoord_removeDup.out.bw")
writeLines("fake bigWig", signal_file)
chip_alignment <- app_env$chip_alignment_summary_table(chip_project)
assert(NROW(chip_alignment) == 1L && all(c("role", "condition", "matched_input") %in% names(chip_alignment)) && chip_alignment$matched_input[[1]] == "I1", "ChIP alignment summary includes experimental roles")
chip_signal <- app_env$peak_signal_track_table(chip_project)
assert(NROW(chip_signal) == 1L && chip_signal$role[[1]] == "chip" && chip_signal$normalization[[1]] == "CPM", "ChIP signal table reports role and saved normalization")
assert(identical(app_env$validated_project_result_path(chip_project, signal_file), normalizePath(signal_file)), "current-project result path accepted")
outside_file <- tempfile("outside-result-")
writeLines("outside", outside_file)
assert(identical(app_env$validated_project_result_path(chip_project, outside_file), ""), "result path outside the current project is rejected")
assert(inherits(app_env$atac_summary_cards_ui(chip_project), "shiny.tag"), "ChIP summary cards render with fake results")
assert(inherits(app_env$chip_results_explorer_ui(), "shiny.tag"), "ChIP Results Explorer UI renders locally")

atac_ui_text <- as.character(app_env$atac_results_explorer_ui())
chip_ui_text <- as.character(app_env$chip_results_explorer_ui())
cutrun_ui_text <- as.character(app_env$cutrun_results_explorer_ui())
for (ui_check in list(
  ATAC = atac_ui_text,
  ChIP = chip_ui_text,
  CUTRUN = cutrun_ui_text
)) {
  assert(grepl("Developed by CSHL's Bioinformatics Shared Resource", ui_check, fixed = TRUE), paste(names(ui_check), "Results Explorer uses the shared branded header"))
  assert(grepl("Overview", ui_check, fixed = TRUE) && grepl("QC", ui_check, fixed = TRUE) && grepl("Files", ui_check, fixed = TRUE), "custom Results Explorer exposes the standard navigation")
  assert(grepl("Signal &amp; Peaks", ui_check, fixed = TRUE) || grepl("Signal & Peaks", ui_check, fixed = TRUE), "custom Results Explorer exposes standardized signal and peak navigation")
}
assert(grepl("Initial QC", chip_ui_text, fixed = TRUE) && grepl("Fragment Size", chip_ui_text, fixed = TRUE), "ChIP Results Explorer includes RNA-style QC navigation")
assert(grepl("Signal Tracks", atac_ui_text, fixed = TRUE) && grepl("Signal Tracks", chip_ui_text, fixed = TRUE), "ATAC and ChIP Results Explorers expose signal-track navigation")

rna_project <- chip_project
rna_project$id <- "fake-rna"
rna_project$name <- "fake-rna"
rna_project$analysis_key <- "rna"
rna_project$analysis <- "RNA-seq"
rna_project$genome <- "mouse"
rna_project$genome_version <- "mouse_gencodeM39"
old_app_home <- app_env$APP_HOME
app_env$APP_HOME <- file.path(root, "fake-app-home")
rna_config <- app_env$write_native_shiny_config(rna_project)
rna_config_text <- paste(readLines(rna_config, warn = FALSE), collapse = "\n")
assert(grepl('genome_species <- "mouse"', rna_config_text, fixed = TRUE), "RNA Results Explorer config records the analysis species")
assert(grepl('genome_version <- "mouse_gencodeM39"', rna_config_text, fixed = TRUE), "RNA Results Explorer config records the analysis reference")
assert(grepl("gencode.vM39.primary_assembly.annotation.gtf", rna_config_text, fixed = TRUE), "RNA Results Explorer config receives the analysis GTF")
rna_viewer <- app_env$load_native_rnaseq_viewer(rna_project)
app_env$APP_HOME <- old_app_home
assert(inherits(rna_viewer$ui, "shiny.tag") && is.function(rna_viewer$server), "RNA Results Explorer loads against a synthetic project")
assert(grepl("RNA-Seq Results Explorer", as.character(rna_viewer$ui), fixed = TRUE), "RNA Results Explorer branded UI is present")

fake_jobs <- data.frame(
  step = c("Bowtie2", "Bowtie2", "Bowtie2", "Bowtie2", "FastQC"),
  sample = c("A1", "A2", "A3", "A4", "A1"),
  job_id = c("101", "102", "103", "104", "105"),
  slurm_state = c("RUNNING", "PENDING", "COMPLETED", "CANCELLED", "RUNNING"),
  stringsAsFactors = FALSE
)
active_bowtie <- app_env$active_step_jobs_from_jobs(fake_jobs, "Bowtie2")
assert(identical(sort(active_bowtie$job_id), c("101", "102")), "active-job filtering excludes completed, cancelled, and other-step jobs")
assert(identical(unname(app_env$active_step_sample_choices(fake_jobs, "Bowtie2")), c("A1", "A2")), "cancellation choices contain only samples with active jobs")
assert(identical(app_env$filter_active_jobs_by_samples(active_bowtie, "A2")$job_id, "102"), "selected-sample cancellation resolves only the requested active job")
assert(inherits(app_env$active_jobs_modal_table(active_bowtie), "shiny.tag"), "active-job cancellation summary renders locally")

assay_jobs <- data.frame(
  step = c("STAR", "Bowtie2", "SEACR", "MACS2 Peaks"),
  sample = c("rna_sample", "atac_sample", "cutrun_sample", "chip_sample"),
  job_id = c("201", "202", "203", "204"),
  slurm_state = rep("RUNNING", 4),
  stringsAsFactors = FALSE
)
for (step in assay_jobs$step) {
  expected <- assay_jobs$sample[assay_jobs$step == step]
  assert(identical(unname(app_env$active_step_sample_choices(assay_jobs, step)), expected), paste(step, "supports active sample cancellation"))
}

assert(system2("bash", c("-n", shQuote(file.path(repo_root, "run_codespringweb.sh")))) == 0L, "CodeSpringApp launcher shell syntax is valid")

bad_q <- app_env$submit_atac_macs2_jobs(atac_project, "not-a-number", "A1")
assert(grepl("q-value must be", bad_q), "invalid ATAC MACS2 q-value rejected before submission")
assert(grepl("two different", app_env$submit_atac_diffbind_job(atac_project, "condition", "A", "A")), "identical ATAC DiffBind conditions rejected")
assert(grepl("two different", app_env$submit_chip_diffbind_job(chip_project, "condition", "A", "A")), "identical ChIP DiffBind conditions rejected")

comparison_dir <- file.path(root, "diffbind", "B_vs_A")
dir.create(comparison_dir, recursive = TRUE)
legacy_result <- file.path(comparison_dir, "DifferentialPeaks_B_vs_A_ref.txt")
writeLines("Fold\tFDR\n1\t0.01", legacy_result)
assert(identical(app_env$peak_diffbind_status(atac_project), "Complete"), "legacy DiffBind comparison remains recognized")
writeLines(character(0), legacy_result)
assert(!app_env$diffbind_comparison_complete(comparison_dir), "empty legacy result is not accepted")
writeLines("Fold\tFDR\n1\t0.01", legacy_result)
writeLines("status\trunning", file.path(comparison_dir, "_RUN_STARTED"))
assert(identical(app_env$peak_diffbind_status(atac_project), "Likely failed"), "partial DiffBind output is not accepted")
assert(!app_env$diffbind_comparison_complete(comparison_dir), "started DiffBind comparison hidden from Results Explorer")
active_jobs <- data.frame(
  step = "Differential Peaks", slurm_state = "RUNNING", sample = basename(comparison_dir),
  target = file.path(comparison_dir, "_COMPLETE"), stringsAsFactors = FALSE
)
assert(app_env$diffbind_comparison_active(atac_project, comparison_dir, jobs = active_jobs), "active DiffBind comparison recognized")
unlink(file.path(comparison_dir, "_RUN_STARTED"))
writeLines("status\tcomplete", file.path(comparison_dir, "_COMPLETE"))
assert(app_env$diffbind_comparison_complete(comparison_dir), "final DiffBind marker accepted")

cat("CodeSpringApp fake-data helper smoke tests passed.\n")
