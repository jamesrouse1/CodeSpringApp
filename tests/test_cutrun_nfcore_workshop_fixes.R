repo_root <- normalizePath(
  ".",
  winslash = "/",
  mustWork = TRUE
)

app_env <- new.env(parent = globalenv())

sys.source(
  file.path(repo_root, "app.R"),
  envir = app_env
)

# app.R sources runtime helpers globally. Re-source them into the isolated
# test environment so dependencies (and the mocks below) resolve locally.
for (helper in c("cutrun_nfcore_submission.R", "cutrun_nfcore_results.R")) {
  sys.source(
    file.path(repo_root, "R", helper),
    envir = app_env
  )
}

assert <- function(value, message) {
  if (!isTRUE(value)) {
    stop("ASSERTION FAILED: ", message, call. = FALSE)
  }
}

expect_error_matching <- function(expr, pattern, message) {
  result <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  assert(
    !is.null(result) && grepl(pattern, result),
    paste0(message, " (got: ", result %||% "no error", ")")
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

root <- tempfile("cutrun_nfcore_workshop_")
dir.create(root, recursive = TRUE)

on.exit(
  unlink(root, recursive = TRUE, force = TRUE),
  add = TRUE
)

# -------------------------------------------------------------------------
# Samplesheet: designs without IgG controls
# -------------------------------------------------------------------------

fastq_dir <- file.path(root, "fastq")
dir.create(fastq_dir)

fastqs <- c(
  "H3K4me3_rep1_R1.fastq.gz", "H3K4me3_rep1_R2.fastq.gz",
  "H3K4me3_rep2_R1.fastq.gz", "H3K4me3_rep2_R2.fastq.gz",
  "IgG_rep1_R1.fastq.gz", "IgG_rep1_R2.fastq.gz"
)
invisible(file.create(file.path(fastq_dir, fastqs)))

sheet_project <- list(
  id = "cutrun/no_igg",
  name = "no_igg",
  analysis = "CUT&RUN",
  analysis_key = "cutrun",
  cutrun_backend = "nfcore",
  paired_end = TRUE,
  fastq_dir = fastq_dir,
  fastq_dirs = fastq_dir,
  data_dir = root,
  design_matrix_path = file.path(fastq_dir, "design_matrix.txt")
)

no_igg_design <- data.frame(
  sample = c("H3K4me3_rep1", "H3K4me3_rep2"),
  cell_type = "Epithelial",
  mark = "H3K4me3",
  target_class = "histone_narrow",
  condition = "Vehicle",
  replicate = c(1, 2),
  control_sample = c("", ""),
  filename = c(
    "H3K4me3_rep1_R1.fastq.gz,H3K4me3_rep1_R2.fastq.gz",
    "H3K4me3_rep2_R1.fastq.gz,H3K4me3_rep2_R2.fastq.gz"
  ),
  stringsAsFactors = FALSE
)

no_igg <- app_env$cutrun_nfcore_samplesheet(sheet_project, no_igg_design)

assert(
  NROW(no_igg$samplesheet) == 2L &&
    all(!nzchar(no_igg$samplesheet$control)),
  "targets without IgG should produce an empty nf-core control column"
)

mixed_design <- rbind(
  no_igg_design,
  data.frame(
    sample = "IgG_rep1",
    cell_type = "Epithelial",
    mark = "IgG",
    target_class = "control",
    condition = "Vehicle",
    replicate = 1,
    control_sample = "",
    filename = "IgG_rep1_R1.fastq.gz,IgG_rep1_R2.fastq.gz",
    stringsAsFactors = FALSE
  )
)
mixed_design$control_sample[[1]] <- "IgG_rep1"

expect_error_matching(
  app_env$cutrun_nfcore_samplesheet(sheet_project, mixed_design),
  "every target, or for none.*H3K4me3_rep2",
  "a mix of controlled and uncontrolled targets should be rejected clearly"
)

# -------------------------------------------------------------------------
# Bundle: a rejected design must not leave a run directory behind
# -------------------------------------------------------------------------

bundle_root <- file.path(root, "bundle")
dir.create(bundle_root)

launcher <- file.path(bundle_root, "nextflow-sarek")
writeLines(c("#!/usr/bin/env bash", "exit 0"), launcher)
Sys.chmod(launcher, mode = "0700")

cluster_config <- file.path(bundle_root, "cshl_slurm.config")
writeLines("process { executor = 'slurm' }", cluster_config)

runtime <- list(
  launcher = launcher,
  config = cluster_config,
  nxf_home = file.path(bundle_root, "nextflow_home"),
  singularity_cache = file.path(bundle_root, "singularity_cache")
)

bundle_project <- list(
  id = "cutrun/bundle",
  name = "bundle",
  analysis = "CUT&RUN",
  analysis_key = "cutrun",
  cutrun_backend = "nfcore",
  data_dir = bundle_root
)

real_samplesheet <- app_env$cutrun_nfcore_samplesheet

app_env$cutrun_nfcore_samplesheet <- function(project, design = NULL) {
  stop("design rejected for test")
}

for (attempt in 1:2) {
  expect_error_matching(
    app_env$cutrun_nfcore_build_bundle(
      bundle_project,
      run_id = "rejected",
      runtime = runtime
    ),
    "design rejected for test",
    paste("attempt", attempt, "should report the design error, not 'run already exists'")
  )
}

assert(
  !dir.exists(file.path(bundle_root, "rejected")),
  "a rejected design left a run directory behind"
)

# A failure after directories exist must also clean up.
app_env$cutrun_nfcore_samplesheet <- function(project, design = NULL) {
  list(
    samplesheet = data.frame(group = "target__x", replicate = 1L, fastq_1 = "a", fastq_2 = "b", control = ""),
    mapping = data.frame(codespring_sample = "x", nfcore_group = "target__x", replicate = 1L)
  )
}
app_env$cutrun_nfcore_params <- function(...) list(input = "x")
app_env$cutrun_nfcore_launch_script <- function(...) stop("launch script failed for test")

expect_error_matching(
  app_env$cutrun_nfcore_build_bundle(
    bundle_project,
    run_id = "late_failure",
    runtime = runtime
  ),
  "launch script failed for test",
  "late bundle failure should surface its own error"
)

assert(
  !dir.exists(file.path(bundle_root, "late_failure")),
  "a late bundle failure left a partial run directory behind"
)

app_env$cutrun_nfcore_samplesheet <- real_samplesheet

# -------------------------------------------------------------------------
# Results: multi-lane samples and MACS2 peak selection
# -------------------------------------------------------------------------

results_root <- file.path(root, "results_project")
dir.create(results_root)

results_project <- list(
  id = "cutrun/results",
  name = "results",
  analysis = "CUT&RUN",
  analysis_key = "cutrun",
  cutrun_backend = "nfcore",
  data_dir = results_root
)

result_paths <- app_env$cutrun_nfcore_results_paths(results_project)
dir.create(dirname(result_paths$mapping_path), recursive = TRUE)

group <- "target__epithelial_h3k4me3_vehicle"

utils::write.table(
  data.frame(
    codespring_sample = c("H3K4me3_rep1", "H3K4me3_rep1", "H3K4me3_rep2"),
    nfcore_group = group,
    replicate = c(1L, 1L, 2L),
    lane = c(1L, 2L, 1L),
    role = "target",
    stringsAsFactors = FALSE
  ),
  result_paths$mapping_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

meta <- app_env$cutrun_nfcore_sample_meta(results_project, "H3K4me3_rep1")

assert(
  !is.null(meta) && identical(as.integer(meta$replicate), 1L),
  "a sample sequenced on two lanes should still resolve to one mapping row"
)

assert(
  NROW(app_env$cutrun_nfcore_result_mapping(results_project)) == 2L,
  "the results mapping should hold one row per sample, not per lane"
)

macs2_dir <- file.path(result_paths$called_peaks, "macs2")
seacr_dir <- file.path(result_paths$called_peaks, "seacr")
dir.create(macs2_dir, recursive = TRUE)
dir.create(seacr_dir, recursive = TRUE)

for (rep in c("R1", "R2")) {
  prefix <- paste0(group, "_", rep)
  invisible(file.create(file.path(
    macs2_dir,
    paste0(prefix, c(".macs2_peaks.narrowPeak", ".macs2_summits.bed", ".macs2.peaks.cut.bed", ".macs2_peaks.xls"))
  )))
  invisible(file.create(file.path(seacr_dir, paste0(prefix, ".seacr.peaks.stringent.bed"))))
}

macs2_peak <- app_env$cutrun_nfcore_peak_path(results_project, "H3K4me3_rep1", "MACS2")

assert(
  identical(basename(macs2_peak), paste0(group, "_R1.macs2_peaks.narrowPeak")),
  paste("MACS2 peaks should resolve to the narrowPeak file, got:", macs2_peak)
)

seacr_peak <- app_env$cutrun_nfcore_peak_path(results_project, "H3K4me3_rep2", "SEACR")

assert(
  identical(basename(seacr_peak), paste0(group, "_R2.seacr.peaks.stringent.bed")),
  paste("SEACR peaks should still resolve, got:", seacr_peak)
)

# -------------------------------------------------------------------------
# Runtime: private NXF_HOME when the shared one is not writable
# -------------------------------------------------------------------------

shared_backend <- file.path(root, "backend")
shared_home <- file.path(shared_backend, "runtime", "nextflow", "cutandrun")
dir.create(shared_home, recursive = TRUE)

assert(
  identical(app_env$cutrun_nfcore_default_nxf_home(shared_backend), shared_home),
  "a writable shared NXF_HOME should be used as-is"
)

assert(
  identical(
    app_env$cutrun_nfcore_default_nxf_home(file.path(root, "missing_backend")),
    file.path(path.expand("~"), ".codespring", "nextflow", "cutandrun")
  ),
  "an unusable shared NXF_HOME should fall back to a private one"
)

# -------------------------------------------------------------------------
# Progress: tolerated or retried failures and per-replicate task attribution
# -------------------------------------------------------------------------

trace <- data.frame(
  label = c(
    "NFCORE_CUTANDRUN:CUTANDRUN:PRESEQ_LCEXTRAP (target__K562_H3K27me3_Untreated_R2)",
    "NFCORE_CUTANDRUN:CUTANDRUN:BOWTIE2_ALIGN (target__K562_H3K27me3_Untreated_R1)",
    "NFCORE_CUTANDRUN:CUTANDRUN:BOWTIE2_ALIGN (target__K562_H3K27me3_Untreated_R1)"
  ),
  stage = c("Preprocessing", "Alignment", "Alignment"),
  task_status = c("FAILED", "FAILED", "COMPLETED"),
  status_class = c("failed", "failed", "complete"),
  stringsAsFactors = FALSE
)

running <- app_env$cutrun_nfcore_resolve_failed_attempts(trace, "RUNNING")
assert(
  identical(app_env$cutrun_nfcore_stage_state(running[running$stage == "Alignment", ]), "Complete"),
  "a failed attempt that a retry completed should not mark its stage Failed"
)
assert(
  identical(app_env$cutrun_nfcore_stage_state(running[running$stage == "Preprocessing", ]), "Failed"),
  "an unresolved failure in a running workflow should still show as Failed"
)

completed <- app_env$cutrun_nfcore_resolve_failed_attempts(trace, "COMPLETED")
assert(
  identical(app_env$cutrun_nfcore_stage_state(completed[completed$stage == "Preprocessing", ]), "Complete"),
  "a failure nf-core tolerated in a completed workflow should not mark its stage Failed"
)

aliases <- list(
  K562_H3K4me3_rep1 = c("K562_H3K4me3_rep1", "target__K562_H3K4me3_Untreated_R1"),
  K562_H3K4me3_rep2 = c("K562_H3K4me3_rep2", "target__K562_H3K4me3_Untreated_R2"),
  K562_IgG_rep1 = c("K562_IgG_rep1", "control__K562_IgG_Untreated_R1"),
  Sample_rep10 = c("Sample_rep10", "target__K562_H3K4me3_Untreated_R10")
)
task_sample <- function(label) app_env$cutrun_nfcore_task_sample(label, aliases)

assert(
  identical(task_sample("NFCORE_CUTANDRUN:CUTANDRUN:FASTQC (target__K562_H3K4me3_Untreated_R2)"), "K562_H3K4me3_rep2"),
  "replicate 2 tasks should be attributed to replicate 2"
)
assert(
  identical(task_sample("NFCORE_CUTANDRUN:CUTANDRUN:FASTQC (target__K562_H3K4me3_Untreated_R10)"), "Sample_rep10"),
  "_R1 must not claim a task labelled _R10"
)
assert(
  identical(task_sample("NFCORE_CUTANDRUN:CUTANDRUN:CONSENSUS_PEAKS:SORT (target__K562_H3K4me3_Untreated)"), ""),
  "group-level tasks should not be attributed to a single replicate"
)

cat("CUT&RUN nf-core workshop fixes: PASS\n")
