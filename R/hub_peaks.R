# CodeSpringHub bridge for bulk ATAC-seq and ChIP-seq (sourced after hub_bridge.R).
#
# Same rule as RNA-seq: the Hub runs exactly what CodeSpringApp runs. Every step calls the app's own
# submit function (submit_cutadapt_jobs, submit_atac_bowtie2_jobs, submit_chip_macs2_jobs, ...), with
# the app's defaults (Nextera adapters for ATAC; MACS2 q 0.05 for ATAC, 0.01 narrow for ChIP), and the
# app's own completion checks. Nothing here computes anything new.
#
#   ATAC-seq  Cutadapt -> FastQC -> Bowtie2 -> MACS2 Peaks -> Differential Peaks (DiffBind)
#   ChIP-seq  Cutadapt -> FastQC -> Bowtie2 -> MACS2 Peaks (targets vs matched input) ->
#             Differential Peaks (DiffBind, targets) -> Peak Annotation (HOMER)

HUB_PEAK_KINDS <- list(
  atac = list(analysis = "ATAC-seq", label = "Bulk ATAC-seq", adapter = "Nextera Transposase ATAC", qvalue = "0.05",
              needs_paired = TRUE, annotate = FALSE),
  chip = list(analysis = "ChIP-seq", label = "Bulk ChIP-seq", adapter = "Illumina Universal TruSeq RNA", qvalue = "0.01",
              needs_paired = FALSE, annotate = TRUE)
)
HUB_PEAK_TYPES <- c("narrow", "broad")

hub_peak_catalog <- function(kind) {
  k <- HUB_PEAK_KINDS[[kind]]
  example <- file.path(CSL_ROOT, "scripts_DoNotTouch", "test")
  list(
    species = c("human", "mouse"),
    references = list(human = "GRCh38 / GENCODE v50 (Bowtie2)", mouse = "GRCm39 / GENCODE M39 (Bowtie2)"),
    adapters_r1 = as.list(adapter_choices_r1()[adapter_choices_r1() != "__custom__"]),
    adapters_r2 = as.list(adapter_choices_r2()[adapter_choices_r2() != "__custom__"]),
    default_adapter = k$adapter,
    macs2 = c(list(default_qvalue = k$qvalue, qvalue_range = "greater than 0 and at most 1"),
              if (identical(kind, "chip")) list(peak_types = HUB_PEAK_TYPES, default_peak_type = "narrow")),
    design = if (identical(kind, "chip")) {
      paste("Columns: sample, grouping columns, reference (chip for targets, input for controls), control_sample",
            "(for each target, the sample name of its matched input row), filename (last).")
    } else "Columns: sample, grouping columns, filename (last). Paired-end reads only (as in CodeSpringLab).",
    differential = "DiffBind on one grouping column, comparison vs reference, at least 2 replicates in each group",
    run_modes = HUB_RUN_MODES,
    notify = list(email_to = hub_email_address(), every_hours = "0 to 24", default = "off"),
    limits = list(max_samples = HUB_MAX_SAMPLES, max_fastq_dirs = HUB_MAX_FASTQ_DIRS, min_replicates_per_group = 2L),
    example = list(fastq_dir = file.path(example, paste0("fastq_", kind)),
                   design_matrix = file.path(example, paste0("manifest_", kind), "design_matrix.txt"),
                   species = "mouse", paired_end = identical(kind, "atac"),
                   diffbind = if (identical(kind, "atac")) list(column = "condition", reference = "day1", comparison = "day3")
                              else list(column = "condition", reference = "control", comparison = "siGATA"))
  )
}

hub_qvalue <- function(x, default, add = function(...) NULL) {
  q <- suppressWarnings(as.numeric(x %||% default))
  if (length(q) != 1 || !is.finite(q) || q <= 0 || q > 1) { add("macs2.qvalue must be a number greater than 0 and at most 1."); return(default) }
  format(q, scientific = FALSE, trim = TRUE)
}

hub_check_peaks <- function(spec, kind) {
  k <- HUB_PEAK_KINDS[[kind]]
  errors <- character(0)
  add <- function(...) errors <<- c(errors, paste0(...))
  b <- hub_check_basics(spec, add, default_adapter = k$adapter)
  if (isTRUE(k$needs_paired) && !isTRUE(b$paired)) add("ATAC-seq needs paired-end reads (CodeSpringLab's ATAC workflow is paired-end only).")
  if (length(b$species) && !b$species %in% c("human", "mouse")) add("ATAC-seq and ChIP-seq run on human or mouse.")
  design <- b$design
  macs2 <- spec$macs2 %||% list()
  qvalue <- hub_qvalue(macs2$qvalue, k$qvalue, add)
  peak_type <- if (identical(kind, "chip")) tolower(trimws(as.character(macs2$peak_type %||% "narrow"))) else "narrow"
  if (!peak_type %in% HUB_PEAK_TYPES) add("macs2.peak_type must be narrow or broad.")

  targets <- design
  if (!is.null(design) && identical(kind, "chip")) {
    for (col in c("reference", "control_sample")) if (!col %in% names(design)) add("A ChIP-seq design needs a '", col, "' column.")
    if (all(c("reference", "control_sample") %in% names(design))) {
      role <- tolower(trimws(design$reference))
      if (any(!role %in% c("chip", "input"))) add("reference must be 'chip' (target) or 'input' (control) in every row.")
      inputs <- design$sample[role == "input"]
      targets <- design[role == "chip", , drop = FALSE]
      if (!NROW(targets)) add("The design has no ChIP target rows (reference = chip).")
      bad_ctrl <- targets$sample[!trimws(targets$control_sample) %in% inputs]
      if (length(bad_ctrl)) add("Each target's control_sample must name an input row; not for: ", paste(head(bad_ctrl, 5), collapse = ", "), ".")
    }
  }
  diff <- spec$diffbind
  if (!is.null(diff) && !is.null(design)) {
    groups <- setdiff(names(design)[-c(1, NCOL(design))], c("reference", "control_sample"))
    col <- as.character(diff$column %||% "")
    if (!col %in% groups) add("diffbind.column must be one of the design's grouping columns: ", paste(groups, collapse = ", "), ".")
    else {
      ref <- as.character(diff$reference %||% ""); cmp <- as.character(diff$comparison %||% "")
      vals <- trimws(as.character(targets[[col]]))
      if (!ref %in% vals || !cmp %in% vals) add("diffbind.reference and diffbind.comparison must be values of '", col, "'",
                                               if (identical(kind, "chip")) " among the ChIP targets" else "", ": ", paste(unique(vals), collapse = ", "), ".")
      else if (identical(ref, cmp)) add("diffbind.reference and diffbind.comparison must differ.")
      else if (sum(vals == ref) < 2 || sum(vals == cmp) < 2) add("DiffBind needs at least 2 replicates in each group (", ref, ": ", sum(vals == ref), ", ", cmp, ": ", sum(vals == cmp), ").")
    }
  }
  annotate <- tryCatch(hub_flag(spec$annotate, k$annotate), error = function(e) { add("annotate: ", conditionMessage(e)); k$annotate })
  if (length(errors)) return(list(ok = FALSE, errors = errors))

  project <- hub_new_project(b, k$analysis)
  if (dir_has_contents(project_result_dir(project))) {
    return(list(ok = FALSE, errors = paste0("A results folder for '", b$name, "' already exists (", project_result_dir(project),
      "). Choose another project name; the Hub never clears existing results.")))
  }
  reads <- hub_check_reads(project, design, b$paired)
  if (is.character(reads) && length(reads) == 1 && startsWith(reads, "ERROR:")) return(list(ok = FALSE, errors = sub("^ERROR: ", "", reads)))
  list(ok = TRUE, pipeline = kind, project = project, design = design, paired = b$paired, species = b$species,
       version = b$version, trim = b$trim, reads = reads, run_mode = b$run_mode, notify = b$notify,
       macs2 = list(qvalue = qvalue, peak_type = peak_type), diffbind = diff, annotate = isTRUE(annotate),
       targets = if (identical(kind, "chip")) targets$sample else design$sample)
}

hub_peak_scripts <- function(kind, paired) {
  pe <- if (isTRUE(paired)) "PE" else "SE"
  list(
    Cutadapt = file.path(SCRIPTS_DIR, if (paired) "cutadapt_PE/qsub_cutadapt_PE.sh" else "cutadapt_SE/qsub_cutadapt_SE.sh"),
    FastQC = file.path(SCRIPTS_DIR, "FastQC", "qsub_fastqc.sh"),
    Bowtie2 = file.path(SCRIPTS_DIR, "bowtie2", if (identical(kind, "atac")) "qsub_bowtie2_PE.sh" else paste0("qsub_bowtie2_chip_", pe, ".sh")),
    `MACS2 Peaks` = file.path(SCRIPTS_DIR, "MACS2", if (identical(kind, "atac")) "qsub_macs2_PE.sh" else "qsub_macs2_chip_SE.sh"),
    `Differential Peaks` = file.path(SCRIPTS_DIR, "DiffBind", if (identical(kind, "atac")) "qsub_diffbind.sh" else "qsub_diffbind_chip.sh"),
    `Peak Annotation` = file.path(SCRIPTS_DIR, "Homer", "qsub_annotate_peak_results.sh")
  )
}

hub_plan_peaks <- function(checked) {
  kind <- checked$pipeline
  p <- checked$project
  scripts <- hub_peak_scripts(kind, checked$paired)
  n <- NROW(checked$design)
  ref <- if (identical(kind, "atac")) atac_reference_resources(p) else chip_reference_resources(p)
  step <- function(name, jobs, what) {
    script <- scripts[[name]]
    list(step = name, jobs = jobs, what = what, script = basename(script),
         resources_per_job = if (file.exists(script)) hub_sbatch_resources(script) else list())
  }
  steps <- list()
  if (isTRUE(checked$trim$enabled)) steps[[length(steps) + 1]] <- step("Cutadapt", n, paste0("Trim adapters (", checked$trim$adapter_r1, "), minimum length ", checked$trim$min_length))
  steps[[length(steps) + 1]] <- step("FastQC", n, paste("Read quality reports on", if (isTRUE(checked$trim$enabled)) "trimmed" else "raw", "reads"))
  steps[[length(steps) + 1]] <- step("Bowtie2", n, paste("Align to", ref$label, "and remove duplicates"))
  steps[[length(steps) + 1]] <- if (identical(kind, "atac")) {
    step("MACS2 Peaks", n, paste0("Call peaks with MACS2 (ATAC shift -100, extension 200, q ", checked$macs2$qvalue, ")"))
  } else {
    step("MACS2 Peaks", length(checked$targets), paste0("Call ", checked$macs2$peak_type, " peaks for each target against its matched input (q ", checked$macs2$qvalue, ")"))
  }
  d <- checked$diffbind
  if (!is.null(d)) steps[[length(steps) + 1]] <- step("Differential Peaks", 1L, paste0("DiffBind: ", d$comparison, " vs ", d$reference, " (reference) on '", d$column, "'"))
  if (isTRUE(checked$annotate)) steps[[length(steps) + 1]] <- step("Peak Annotation", 1L, "Annotate peaks to the nearest genes with HOMER")
  roles <- if (identical(kind, "chip")) as.list(table(tolower(checked$design$reference))) else NULL
  list(
    pipeline = kind,
    project_id = p$id, project_name = p$name, results_dir = project_result_dir(p),
    reference = ref$label, paired_end = checked$paired, samples = checked$design$sample,
    roles = roles,
    groups = lapply(checked$design[, setdiff(names(checked$design)[-c(1, NCOL(checked$design))], c("reference", "control_sample")), drop = FALSE],
                    function(x) as.list(table(x))),
    steps = steps,
    run_mode = checked$run_mode, run_mode_what = HUB_RUN_MODE_WHAT[[checked$run_mode]],
    notify = if (isTRUE(checked$notify$email)) list(email_to = hub_email_address(), every_hours = checked$notify$every_hours) else "off",
    batching = list(max_jobs_at_once = hub_submit_limit(),
      note = paste0("Per-sample jobs are submitted in waves of at most ", max(1L, hub_submit_limit() - HUB_SLOT_RESERVE), " so you stay under your Slurm limit.")),
    controller = list(what = "A small background process on the login node that submits each step when the previous one has finished (it uses no Slurm resources itself)"),
    runs_as = CURRENT_USER, codespringapp = APP_ROOT, codespringlab = CSL_ROOT
  )
}

hub_recheck_peaks <- function(spec, add) {
  # The controller re-checks what it will pass to the app (the stored file could have been edited).
  kind <- spec$pipeline
  k <- HUB_PEAK_KINDS[[kind]]
  if (is.null(k)) { add("unknown pipeline."); return(spec) }
  spec$macs2$qvalue <- hub_qvalue(spec$macs2$qvalue, k$qvalue, add)
  if (!(spec$macs2$peak_type %||% "narrow") %in% HUB_PEAK_TYPES) add("peak type outside narrow/broad.")
  if (!is.null(spec$trim$adapter_r1) && !spec$trim$adapter_r1 %in% names(adapter_choices_r1())) add("adapter outside the app's list.")
  spec
}

hub_peak_slug <- function(kind, d) {
  # The folder name the app gives a DiffBind comparison (submit_atac_diffbind_job / submit_chip_diffbind_job).
  if (identical(kind, "atac")) clean_name(paste0(d$comparison, "_vs_", d$reference), "comparison")
  else clean_name(paste(c(clean_name(d$comparison, "comparison"), "vs", clean_name(d$reference, "reference")), collapse = "_"), "comparison")
}

hub_peak_entries <- function(p, spec) {
  kind <- spec$pipeline
  k <- HUB_PEAK_KINDS[[kind]]
  trim <- isTRUE(hub_flag(spec$trim$enabled, TRUE))
  a1 <- adapter_choices_r1()[[spec$trim$adapter_r1 %||% k$adapter]]
  a2 <- adapter_choices_r2()[[spec$trim$adapter_r2 %||% spec$trim$adapter_r1 %||% k$adapter]]
  min_len <- as.integer(spec$trim$min_length %||% 20L)
  q <- spec$macs2$qvalue %||% k$qvalue
  peak_type <- spec$macs2$peak_type %||% "narrow"
  d <- spec$diffbind
  all_samples <- function() project_samples(p)
  targets <- function() if (identical(kind, "chip")) as.character(chip_target_design(p)$sample) else all_samples()
  every <- function(step, samples) { s <- samples(); length(s) > 0 && setequal(completed_samples_for_step(p, step, s), s) }
  list(
    Cutadapt = list(run = function(s) submit_cutadapt_jobs(p, a1, a2, min_len, s), done = function() every("Cutadapt", all_samples)),
    FastQC = list(run = function(s) submit_fastqc_jobs(p, trimmed = trim, samples = s), done = function() every("FastQC", all_samples)),
    Bowtie2 = list(run = function(s) if (identical(kind, "atac")) submit_atac_bowtie2_jobs(p, trim, s) else submit_chip_bowtie2_jobs(p, trim, s),
                   done = function() every("Bowtie2", all_samples)),
    `MACS2 Peaks` = list(samples = targets,
                         run = function(s) if (identical(kind, "atac")) submit_atac_macs2_jobs(p, q, s) else submit_chip_macs2_jobs(p, q, peak_type, s),
                         done = function() every("MACS2 Peaks", targets)),
    `Differential Peaks` = list(run = function(s) if (identical(kind, "atac")) submit_atac_diffbind_job(p, d$column, d$reference, d$comparison)
                                                  else submit_chip_diffbind_job(p, d$column, d$reference, d$comparison),
                                explain = function() {
                                  # The app refuses DiffBind when a sample has no peaks; say which, in plain words.
                                  counts <- vapply(targets(), function(x) { m <- hub_macs2_stats(p, x); if (is.null(m$peaks)) NA_real_ else m$peaks }, numeric(1))
                                  empty <- names(counts)[!is.na(counts) & counts == 0]
                                  if (length(empty)) paste0("MACS2 found 0 peaks for: ", paste(empty, collapse = ", "),
                                    ". DiffBind needs peaks in every sample; this usually means too few reads or a weak signal (ask BSR).") else ""
                                },
                                done = function() diffbind_comparison_complete(file.path(p$data_dir, "diffbind", hub_peak_slug(kind, d)))),
    `Peak Annotation` = list(run = function(s) submit_peak_annotation_job(p), done = function() peak_annotation_is_current(p))
  )
}

# ---- QC from each tool's own output ----------------------------------------------------------------

hub_bowtie2_stats <- function(p, sample) {
  x <- metric_file_to_named_list(file.path(p$data_dir, "bowtie2", sample, paste0(sample, "_alignment_summary.txt")))
  if (!length(x)) return(NULL)
  mapped <- suppressWarnings(as.numeric(gsub(",", "", x$mapped_reads %||% NA)))
  dedup <- suppressWarnings(as.numeric(gsub(",", "", x$deduplicated_reads %||% NA)))
  list(mapped_reads = if (is.finite(mapped)) mapped else NULL, deduplicated_reads = if (is.finite(dedup)) dedup else NULL,
       kept_after_dedup_pct = if (is.finite(mapped) && mapped > 0 && is.finite(dedup)) round(100 * dedup / mapped, 1) else NULL)
}

hub_macs2_stats <- function(p, sample) {
  dir <- file.path(p$data_dir, "macs2", sample)
  x <- metric_file_to_named_list(file.path(dir, paste0(sample, "_macs2_summary.txt")))
  count <- suppressWarnings(as.numeric(gsub(",", "", x$peak_count %||% NA)))
  if (!is.finite(count)) {
    peak <- list.files(dir, pattern = "_peaks\\.(narrowPeak|broadPeak)$", full.names = TRUE)
    if (length(peak) == 1) count <- length(readLines(peak, warn = FALSE))
  }
  frip <- suppressWarnings(as.numeric(sub("%$", "", x$frip %||% x$FRiP %||% NA)))
  out <- list(peaks = if (is.finite(count)) count else NULL, frip_pct = if (is.finite(frip)) (if (frip <= 1) round(100 * frip, 1) else frip) else NULL)
  if (length(Filter(Negate(is.null), out))) out else NULL
}
