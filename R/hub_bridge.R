# CodeSpringHub bridge: run CodeSpringApp pipelines without the browser.
#
# CodeSpringHub (the BSR MCP server) lets researchers start analyses from a Claude chat. It must run
# them exactly the way CodeSpringApp does: the same CodeSpringLab scripts, references, resource
# requests, folders and job records. So instead of re-implementing anything, this file drives the
# app's own functions (new_project_from_inputs, submit_*_jobs, job_history, ...). It is sourced into
# the same environment as app.R by bin/codespring-run.
#
# Guardrails live here, next to the code they protect:
#   * only whitelisted pipelines, references, adapters, gene sets and options are accepted;
#   * results always go to the user's own CodeSpringApp results root, never an arbitrary path;
#   * an existing results folder is never cleared or overwritten (the app's UI can; this can't);
#   * the CodeSpringLab scripts and runners are chosen by the app, never by the caller.
#
# Bulk RNA-seq is the first pipeline wired up; others are listed in the catalog as "app only".
#
# A run can go straight through ("auto") or pause after every step ("step") so the person can look
# at that step's QC and decide: continue, add an optional step (Kallisto, RSEM: the app's own), skip
# one, or cancel. Per-sample jobs are submitted in waves that stay under the person's Slurm job
# limit. Optional emails go only to the person's own address.

HUB_SPEC_VERSION <- 2L
HUB_MAX_SAMPLES <- local({ n <- suppressWarnings(as.integer(Sys.getenv("CSL_HUB_MAX_SAMPLES", unset = "96"))); if (is.na(n) || n < 2) 96L else n })
HUB_MAX_FASTQ_DIRS <- 10L
HUB_POLL_SECONDS <- as.numeric(Sys.getenv("CSL_HUB_POLL_SECONDS", unset = "30"))
HUB_NAME_RE <- "^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$"
HUB_COLUMN_RE <- "^[A-Za-z][A-Za-z0-9_]{0,40}$"
HUB_FAILED_STATES <- c("FAILED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY", "CANCELLED", "CANCELLED+", "CA", "PREEMPTED", "BOOT_FAIL")

# Optional steps are CodeSpringApp's own optional RNA-seq quantifiers, run exactly as the app runs them.
HUB_OPTIONAL <- list(
  kallisto = list(step = "Kallisto", after = "FastQC",
                  what = "Transcript abundance (TPM) with Kallisto pseudo-alignment of the same reads"),
  rsem = list(step = "RSEM", after = "STAR",
              what = "Gene and transcript quantification with RSEM from STAR's transcriptome alignments")
)
HUB_STEP_ORDER <- c("Cutadapt", "FastQC", "STAR", "featureCounts", "Kallisto", "RSEM", "DESeq2", "GSEA")
HUB_PER_SAMPLE_STEPS <- c("Cutadapt", "FastQC", "STAR", "featureCounts", "Kallisto", "RSEM")
HUB_SKIPPABLE <- c("Kallisto", "RSEM", "GSEA")
HUB_MODES <- c(auto = "Run every step one after another", step = "Pause after each step: show its QC and wait for the go-ahead")
HUB_EMAIL_DOMAIN <- Sys.getenv("CSL_HUB_EMAIL_DOMAIN", unset = "cshl.edu")
HUB_SLOT_RESERVE <- 2L  # Slurm submit slots kept free (matrix-build jobs, the person's other work)
HUB_DONE_STATES <- c("COMPLETED", "COMPLETED+", "CD")

hub_app_step <- function(step) switch(step, Kallisto = "Kallisto (optional)", RSEM = "RSEM (optional)", step)
hub_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
hub_email_address <- function() paste0(CURRENT_USER, "@", HUB_EMAIL_DOMAIN)

hub_error <- function(...) structure(class = c("hub_error", "error", "condition"), list(message = paste0(...), call = NULL))

# ---- what can be run --------------------------------------------------------------------------

hub_sbatch_resources <- function(script) {
  # The #SBATCH header of a CodeSpringLab wrapper: what each job of that step asks Slurm for.
  if (!file.exists(script)) return(list())
  lines <- grep("^#SBATCH", readLines(script, warn = FALSE), value = TRUE)
  out <- list()
  for (line in lines) {
    opt <- trimws(sub("^#SBATCH", "", line))
    kv <- regmatches(opt, regexec("^--?([A-Za-z-]+)[= ]?(.*)$", opt))[[1]]
    if (length(kv) == 3 && kv[2] %in% c("cpus-per-task", "mem", "mem-per-cpu", "time", "partition", "c", "t", "p")) {
      key <- switch(kv[2], c = "cpus-per-task", t = "time", p = "partition", kv[2])
      out[[key]] <- trimws(kv[3])
    }
  }
  out
}

hub_rna_step_scripts <- function(paired = TRUE) {
  pe <- if (isTRUE(paired)) "PE" else "SE"
  list(
    Cutadapt = file.path(SCRIPTS_DIR, if (paired) "cutadapt_PE/qsub_cutadapt_PE.sh" else "cutadapt_SE/qsub_cutadapt_SE.sh"),
    FastQC = file.path(SCRIPTS_DIR, "FastQC", "qsub_fastqc.sh"),
    STAR = file.path(SCRIPTS_DIR, "STAR", paste0("qsub_star_", pe, ".sh")),
    featureCounts = file.path(SCRIPTS_DIR, "featureCounts", paste0("qsub_featurecounts_", pe, ".sh")),
    Kallisto = file.path(SCRIPTS_DIR, "Kallisto", paste0("qsub_kallisto_", pe, ".sh")),
    RSEM = file.path(SCRIPTS_DIR, "RSEM", paste0("qsub_RSEM_", pe, ".sh")),
    DESeq2 = file.path(SCRIPTS_DIR, "DESeq2", "qsub_deseq2.sh"),
    GSEA = ""
  )
}

hub_catalog <- function() {
  rna_refs <- genome_reference_catalog()
  refs <- lapply(names(rna_refs), function(species) {
    lapply(names(rna_refs[[species]]), function(key) list(species = species, genome_version = key, label = rna_refs[[species]][[key]]$label))
  })
  list(
    spec_version = HUB_SPEC_VERSION,
    codespringapp = APP_ROOT,
    codespringlab = CSL_ROOT,
    results_root = DEFAULT_RESULTS_ROOT,
    pipelines = list(
      list(id = "rnaseq", name = "Bulk RNA-seq (CodeSpringLab: cutadapt, FastQC, STAR, featureCounts, DESeq2, GSEApy)",
           runnable_from_hub = TRUE, results_explorer = "RNA-seq Results Explorer"),
      list(id = "atac", name = "Bulk ATAC-seq", runnable_from_hub = FALSE, results_explorer = "ATAC Results Explorer"),
      list(id = "chip", name = "Bulk ChIP-seq", runnable_from_hub = FALSE, results_explorer = "ChIP Results Explorer"),
      list(id = "cutrun", name = "CUT&RUN (nf-core/cutandrun 3.2.2 via CodeSpringApp, or native)", runnable_from_hub = FALSE, results_explorer = "CUT&RUN Results Explorer"),
      list(id = "scrna", name = "Single-cell RNA-seq (Seurat or Scanpy)", runnable_from_hub = FALSE, results_explorer = "scRNA Results Explorer"),
      list(id = "sarek", name = "Whole-genome germline/somatic variants (nf-core/sarek 3.9.0 via CodeSpringApp)", runnable_from_hub = FALSE, results_explorer = "Sarek Explorer"),
      list(id = "fetchngs", name = "Public data download (nf-core/fetchngs 1.12.0 via CodeSpringApp)", runnable_from_hub = FALSE, results_explorer = "FetchNGS Outputs")
    ),
    rnaseq = list(
      references = unlist(refs, recursive = FALSE),
      default_reference = list(human = "human_gencode50", mouse = "mouse_gencodeM39"),
      adapters_r1 = as.list(adapter_choices_r1()[adapter_choices_r1() != "__custom__"]),
      adapters_r2 = as.list(adapter_choices_r2()[adapter_choices_r2() != "__custom__"]),
      features = c("gene_name", "gene_id"),
      gsea_genesets = GSEAPY_GENESET_OPTIONS,
      optional_steps = lapply(names(HUB_OPTIONAL), function(id) list(id = id, step = HUB_OPTIONAL[[id]]$step,
        what = HUB_OPTIONAL[[id]]$what, needs = HUB_OPTIONAL[[id]]$after, can_add = "when planning, at any pause, or after the run finishes")),
      modes = as.list(HUB_MODES),
      notify = list(email_to = hub_email_address(), every_hours = "0 (only when a step finishes, fails or pauses) to 24", default = "off"),
      batching = list(max_jobs_at_once = hub_submit_limit(), note = "Per-sample jobs go out in waves that stay under your Slurm job limit."),
      limits = list(max_samples = HUB_MAX_SAMPLES, max_fastq_dirs = HUB_MAX_FASTQ_DIRS, min_replicates_per_group = 2L),
      example = list(fastq_dir = RNA_EXAMPLE_FASTQ_DIR, design_matrix = file.path(RNA_EXAMPLE_DESIGN_DIR, "design_matrix.txt"),
                     species = "mouse", paired_end = TRUE, deseq2 = list(column = "treatment", reference = "control", comparison = "treated"))
    )
  )
}

# ---- checking a request -----------------------------------------------------------------------

hub_flag <- function(x, default) {
  if (is.null(x)) return(default)
  if (is.logical(x) && length(x) == 1 && !is.na(x)) return(x)
  hub_error_value <- tolower(trimws(as.character(x)))
  if (hub_error_value %in% c("true", "yes", "y", "1", "paired", "pe")) return(TRUE)
  if (hub_error_value %in% c("false", "no", "n", "0", "single", "se")) return(FALSE)
  stop(hub_error("Expected yes/no, got: ", x))
}

hub_read_design <- function(spec) {
  if (!is.null(spec[["design"]][["rows"]])) {
    cols <- as.character(unlist(spec[["design"]][["columns"]]))
    rows <- spec[["design"]][["rows"]]
    df <- as.data.frame(do.call(rbind, lapply(rows, function(r) as.character(unlist(r)))), stringsAsFactors = FALSE)
    if (NCOL(df) != length(cols)) stop(hub_error("Every design row needs one value per column (", paste(cols, collapse = ", "), ")."))
    names(df) <- cols
    return(df)
  }
  path <- trimws(as.character(spec$design_matrix_path %||% ""))
  if (!nzchar(path)) stop(hub_error("Give the sample design: either design.columns/design.rows or design_matrix_path."))
  path <- normalizePath(path.expand(path), winslash = "/", mustWork = FALSE)
  if (dir.exists(path)) path <- file.path(path, "design_matrix.txt")
  if (!file.exists(path) || file.access(path, 4) != 0) stop(hub_error("Can't read the design matrix: ", path))
  utils::read.delim(path, colClasses = "character", check.names = FALSE)
}

hub_check_rnaseq <- function(spec) {
  errors <- character(0)
  add <- function(...) errors <<- c(errors, paste0(...))
  cat_ <- hub_catalog()$rnaseq

  name <- trimws(as.character(spec$project_name %||% ""))
  if (!grepl(HUB_NAME_RE, name)) add("project_name must be 1-64 letters, digits, '.', '_' or '-' and start with a letter or digit.")
  if (!is.null(spec$results_root)) add("results_root can't be chosen: results always go to your CodeSpringApp results folder (", DEFAULT_RESULTS_ROOT, ").")

  fastq_dirs <- unique(normalizePath(path.expand(as.character(unlist(spec$fastq_dirs %||% character(0)))), winslash = "/", mustWork = FALSE))
  if (!length(fastq_dirs)) add("Give at least one FASTQ folder (fastq_dirs).")
  if (length(fastq_dirs) > HUB_MAX_FASTQ_DIRS) add("At most ", HUB_MAX_FASTQ_DIRS, " FASTQ folders.")
  unreadable <- fastq_dirs[!dir.exists(fastq_dirs) | file.access(fastq_dirs, 5) != 0]
  if (length(unreadable)) add("These FASTQ folders don't exist or you can't read them: ", paste(unreadable, collapse = ", "))

  paired <- tryCatch(hub_flag(spec$paired_end, TRUE), error = function(e) { add("paired_end: ", conditionMessage(e)); TRUE })
  species <- tolower(trimws(as.character(spec$species %||% "")))
  version <- trimws(as.character(spec$genome_version %||% ""))
  refs <- genome_reference_catalog()
  if (!species %in% names(refs)) add("species must be one of: ", paste(names(refs), collapse = ", "), ".")
  else {
    if (!nzchar(version)) version <- cat_$default_reference[[species]] %||% names(refs[[species]])[[1]]
    if (!version %in% names(refs[[species]])) add("genome_version for ", species, " must be one of: ", paste(names(refs[[species]]), collapse = ", "), ".")
  }

  design <- tryCatch(hub_read_design(spec), error = function(e) { add(conditionMessage(e)); NULL })
  if (!is.null(design)) {
    cols <- names(design)
    if (length(cols) < 3 || cols[[1]] != "sample" || cols[[length(cols)]] != "filename") {
      add("The design needs 'sample' as the first column, 'filename' as the last, and at least one grouping column in between.")
    }
    bad_cols <- cols[!grepl(HUB_COLUMN_RE, cols)]
    if (length(bad_cols)) add("Column names must be simple words (letters, digits, _): ", paste(bad_cols, collapse = ", "))
    if (NROW(design) < 2) add("The design needs at least two samples.")
    if (NROW(design) > HUB_MAX_SAMPLES) add("At most ", HUB_MAX_SAMPLES, " samples per run from Claude; ask BSR for larger studies.")
    if ("sample" %in% cols) {
      bad <- design$sample[!grepl(HUB_NAME_RE, design$sample)]
      if (length(bad)) add("Sample names must be letters, digits, '.', '_' or '-': ", paste(head(bad, 5), collapse = ", "))
      dup <- unique(design$sample[duplicated(design$sample)])
      if (length(dup)) add("Duplicate sample names: ", paste(dup, collapse = ", "))
    }
    if ("filename" %in% cols && any(grepl("/", design$filename))) add("filename must name files inside the FASTQ folders, not paths.")
  }

  trim <- spec$trim %||% list()
  trim_on <- tryCatch(hub_flag(trim$enabled, TRUE), error = function(e) { add("trim.enabled: ", conditionMessage(e)); TRUE })
  a1_name <- trim$adapter_r1 %||% "Illumina Universal TruSeq RNA"
  a2_name <- trim$adapter_r2 %||% a1_name
  if (trim_on) {
    if (!a1_name %in% names(cat_$adapters_r1)) add("trim.adapter_r1 must be one of: ", paste(names(cat_$adapters_r1), collapse = "; "), ".")
    if (paired && !a2_name %in% names(cat_$adapters_r2)) add("trim.adapter_r2 must be one of: ", paste(names(cat_$adapters_r2), collapse = "; "), ".")
  }
  min_len <- suppressWarnings(as.integer(trim$min_length %||% 20L))
  if (is.na(min_len) || min_len < 10 || min_len > 100) add("trim.min_length must be a whole number from 10 to 100.")

  feature <- spec$feature %||% "gene_name"
  if (!feature %in% cat_$features) add("feature must be gene_name or gene_id.")

  de <- spec$deseq2
  if (!is.null(de) && !is.null(design)) {
    col <- as.character(de$column %||% "")
    if (!col %in% names(design)[-c(1, NCOL(design))]) add("deseq2.column must be one of the design's grouping columns: ", paste(names(design)[-c(1, NCOL(design))], collapse = ", "), ".")
    else {
      ref <- as.character(de$reference %||% ""); cmp <- as.character(de$comparison %||% "")
      vals <- design[[col]]
      if (!ref %in% vals || !cmp %in% vals) add("deseq2.reference and deseq2.comparison must be values of '", col, "': ", paste(unique(vals), collapse = ", "), ".")
      else if (identical(ref, cmp)) add("deseq2.reference and deseq2.comparison must differ.")
      else if (sum(vals == ref) < 2 || sum(vals == cmp) < 2) add("DESeq2 needs at least 2 replicates in each group (", ref, ": ", sum(vals == ref), ", ", cmp, ": ", sum(vals == cmp), ").")
    }
  }
  gsea <- spec$gsea
  if (!is.null(gsea)) {
    if (is.null(de)) add("GSEA needs a DESeq2 comparison.")
    if (!(gsea$geneset %||% "") %in% GSEAPY_GENESET_OPTIONS) add("gsea.geneset must be one of: ", paste(GSEAPY_GENESET_OPTIONS, collapse = ", "), ".")
  }
  mode <- tolower(trimws(as.character(spec$mode %||% "auto")))
  if (!mode %in% names(HUB_MODES)) add("mode must be 'auto' (run straight through) or 'step' (pause after each step).")
  optional <- hub_optional_ids(spec$optional, add)
  notify <- hub_check_notify(spec$notify, add)
  if (length(errors)) return(list(ok = FALSE, errors = errors))

  project <- new_project_from_inputs(list(
    new_project_analysis = "RNA-seq", new_project_name = name, new_project_mode = "new",
    new_results_root = DEFAULT_RESULTS_ROOT, new_fastq_location_mode = "multiple",
    new_fastq_dirs = paste(fastq_dirs, collapse = "\n"), new_paired_end = if (paired) "paired" else "single",
    new_species = species, new_genome_version = version
  ))
  if (length(optional) && !rna_optional_quantifiers_available(project)) {
    return(list(ok = FALSE, errors = "Kallisto and RSEM aren't available for this reference."))
  }
  if (dir_has_contents(project_result_dir(project))) {
    return(list(ok = FALSE, errors = paste0("A results folder for '", name, "' already exists (", project_result_dir(project),
      "). Choose another project name; the Hub never clears existing results.")))
  }

  # Check every read file resolves, using the app's own pairing logic on a scratch copy of the design.
  scratch <- tempfile("hub_design_"); dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  probe <- project; probe$design_matrix_path <- file.path(scratch, "design_matrix.txt")
  utils::write.table(design, probe$design_matrix_path, sep = "\t", quote = FALSE, row.names = FALSE)
  pairs <- tryCatch(sample_fastq_pairs(probe, FALSE), error = function(e) e)
  if (inherits(pairs, "error")) return(list(ok = FALSE, errors = conditionMessage(pairs)))
  reads <- unique(unlist(lapply(c(pairs$r1, if (paired) pairs$r2), split_fastq_path_list)))
  missing <- reads[!file.exists(reads)]
  if (length(missing)) return(list(ok = FALSE, errors = paste0("These read files from the design aren't in the FASTQ folders: ", paste(head(basename(missing), 10), collapse = ", "))))

  list(ok = TRUE, project = project, design = design, paired = paired, species = species, version = version,
       trim = list(enabled = trim_on, adapter_r1 = a1_name, adapter_r2 = if (paired) a2_name else "", min_length = min_len),
       feature = feature, deseq2 = de, gsea = gsea, reads = reads, mode = mode, optional = optional, notify = notify)
}

hub_optional_ids <- function(x, add = function(...) NULL) {
  if (is.null(x) || !length(x)) return(character(0))
  ids <- if (is.list(x) && !is.null(names(x))) {
    names(x)[vapply(x, function(v) isTRUE(tryCatch(hub_flag(v, FALSE), error = function(e) FALSE)), logical(1))]
  } else tolower(trimws(as.character(unlist(x))))
  bad <- setdiff(ids, names(HUB_OPTIONAL))
  if (length(bad)) add("optional steps must be from: ", paste(names(HUB_OPTIONAL), collapse = ", "), " (got: ", paste(bad, collapse = ", "), ").")
  intersect(names(HUB_OPTIONAL), ids)
}

hub_check_notify <- function(x, add = function(...) NULL) {
  off <- list(email = FALSE, every_hours = 0)
  if (is.null(x)) return(off)
  if (!is.list(x)) { add("notify must look like {\"email\": true, \"every_hours\": 6}."); return(off) }
  extra <- setdiff(names(x), c("email", "every_hours"))
  if (length(extra)) add("notify only takes 'email' and 'every_hours'; emails always go to your own address (", hub_email_address(), ").")
  email <- tryCatch(hub_flag(x$email, FALSE), error = function(e) { add("notify.email: ", conditionMessage(e)); FALSE })
  every <- suppressWarnings(as.numeric(x$every_hours %||% 0))
  if (length(every) != 1 || is.na(every) || (every != 0 && (every < 1 || every > 24))) {
    add("notify.every_hours must be 0 (only when a step finishes, fails or pauses) or 1 to 24."); every <- 0
  }
  list(email = isTRUE(email), every_hours = every)
}

hub_plan <- function(checked) {
  p <- checked$project
  scripts <- hub_rna_step_scripts(checked$paired)
  n <- NROW(checked$design)
  step <- function(name, jobs, what, script = scripts[[name]]) {
    list(step = name, jobs = jobs, what = what, script = if (nzchar(script)) basename(script) else "",
         resources_per_job = if (nzchar(script)) hub_sbatch_resources(script) else list(note = "small helper job"))
  }
  steps <- list()
  if (isTRUE(checked$trim$enabled)) steps[[length(steps) + 1]] <- step("Cutadapt", n, paste0("Trim adapters (", checked$trim$adapter_r1, "), minimum length ", checked$trim$min_length))
  steps[[length(steps) + 1]] <- step("FastQC", n, paste("Read quality reports on", if (isTRUE(checked$trim$enabled)) "trimmed" else "raw", "reads"))
  steps[[length(steps) + 1]] <- step("STAR", n, paste("Align to", genome_resources(p)$label))
  steps[[length(steps) + 1]] <- step("featureCounts", n + 1L, paste0("Count reads per gene (", checked$feature, ") and build the count matrix"))
  for (id in checked$optional) steps[[length(steps) + 1]] <- step(HUB_OPTIONAL[[id]]$step, n + 1L, paste("Optional:", HUB_OPTIONAL[[id]]$what))
  if (!is.null(checked$deseq2)) steps[[length(steps) + 1]] <- step("DESeq2", 1L, paste0("Differential expression: ", checked$deseq2$comparison, " vs ", checked$deseq2$reference, " (reference) on '", checked$deseq2$column, "'"))
  if (!is.null(checked$gsea)) steps[[length(steps) + 1]] <- step("GSEA", 1L, paste("Pathway enrichment with", checked$gsea$geneset))
  list(
    pipeline = "rnaseq",
    project_id = p$id,
    project_name = p$name,
    results_dir = project_result_dir(p),
    reference = genome_resources(p)$label,
    paired_end = checked$paired,
    samples = checked$design$sample,
    groups = lapply(checked$design[, -c(1, NCOL(checked$design)), drop = FALSE], function(x) as.list(table(x))),
    steps = steps,
    mode = list(id = checked$mode, what = HUB_MODES[[checked$mode]]),
    optional_not_chosen = lapply(setdiff(names(HUB_OPTIONAL), checked$optional), function(id)
      list(id = id, what = HUB_OPTIONAL[[id]]$what, note = "Can be added at any pause or after the run finishes")),
    notify = if (isTRUE(checked$notify$email)) list(email_to = hub_email_address(), every_hours = checked$notify$every_hours,
      when = paste0("when each step finishes, fails or pauses", if (checked$notify$every_hours > 0) paste0(", and a summary every ", checked$notify$every_hours, " h") else "")) else "off",
    batching = list(max_jobs_at_once = hub_submit_limit(),
      note = paste0("Per-sample jobs are submitted in waves of at most ", max(1L, hub_submit_limit() - HUB_SLOT_RESERVE), " so you stay under your Slurm limit.")),
    controller = list(what = "A small background process on the login node that submits each step when the previous one has finished (it uses no Slurm resources itself)"),
    runs_as = CURRENT_USER,
    codespringapp = APP_ROOT,
    codespringlab = CSL_ROOT
  )
}

# ---- starting a run ---------------------------------------------------------------------------

hub_status_path <- function(project) file.path(project_result_dir(project), "log", "codespring_hub_run.json")
hub_spec_path <- function(project) file.path(project_result_dir(project), "log", "codespring_hub_spec.json")
hub_cancel_flag <- function(project) file.path(project_result_dir(project), "log", "codespring_hub_cancel")
hub_cancel_requested <- function(project) file.exists(hub_cancel_flag(project))

hub_write_json <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null"), tmp)
  invisible(file.rename(tmp, path))
}

hub_read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)

hub_update_status <- function(project, ...) {
  path <- hub_status_path(project)
  st <- if (file.exists(path)) hub_read_json(path) else list()
  changes <- list(...)
  for (k in names(changes)) st[k] <- list(changes[[k]])  # NULL clears a field
  st$updated <- hub_now()
  hub_write_json(st, path)
  invisible(st)
}

hub_set_step <- function(project, name, state, message = NULL) {
  path <- hub_status_path(project)
  st <- hub_read_json(path)
  for (i in seq_along(st$steps)) if (identical(st$steps[[i]]$step, name)) {
    st$steps[[i]]$state <- state
    if (identical(state, "running")) { st$steps[[i]]$started <- hub_now(); st$steps[[i]]$message <- NULL }
    if (state %in% c("complete", "failed", "cancelled")) st$steps[[i]]$finished <- hub_now()
    if (!is.null(message)) st$steps[[i]]$message <- message
  }
  st$updated <- hub_now()
  hub_write_json(st, path)
}

hub_start <- function(spec, starter = NULL) {
  checked <- hub_check_rnaseq(spec)
  if (!isTRUE(checked$ok)) return(checked)
  plan <- hub_plan(checked)
  p <- checked$project
  log_dir <- file.path(project_result_dir(p), "log")
  dir.create(file.path(p$data_dir, "manifest"), recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  p$design_matrix_path <- file.path(p$data_dir, "manifest", "design_matrix.txt")
  utils::write.table(checked$design, p$design_matrix_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cfg <- write_project_config(p)
  # Store the checked choices, so the controller never re-reads anything the caller could have changed.
  spec$mode <- checked$mode
  spec$optional <- as.list(checked$optional)
  spec$notify <- checked$notify
  spec$checked_at <- hub_now()
  hub_write_json(spec, hub_spec_path(p))
  hub_write_json(list(
    pipeline = "rnaseq", project_id = p$id, project_name = p$name, results_dir = project_result_dir(p),
    state = "queued", started_by = "CodeSpringHub", user = CURRENT_USER, mode = checked$mode,
    steps = lapply(plan$steps, function(s) list(step = s$step, state = "waiting", what = s$what)),
    created = hub_now()
  ), hub_status_path(p))

  launched <- hub_launch_controller(p, starter)
  if (!isTRUE(launched$ok)) {
    # Nothing was started: remove only what this call just created (the folder was empty before).
    unlink(project_result_dir(p), recursive = TRUE)
    unlink(cfg)
    return(launched)
  }
  list(ok = TRUE, project_id = p$id, results_dir = project_result_dir(p), controller_pid = launched$pid,
       controller_host = launched$host, status_file = hub_status_path(p), project_config = cfg, plan = plan)
}

hub_launch_controller <- function(p, starter = NULL) {
  log_dir <- file.path(project_result_dir(p), "log")
  out <- file.path(log_dir, "output_hub_controller.txt")
  err <- file.path(log_dir, "error_hub_controller.txt")
  runner <- file.path(APP_ROOT, "bin", "codespring-run")
  # Like CodeSpringApp's Sarek and nf-core CUT&RUN controllers, the controller is a small detached
  # process on this login node: it mostly sleeps, and it needs the app's R packages, which the
  # compute nodes don't have. setsid + nohup keep it running after Claude Desktop disconnects;
  # setsid also makes it a process group, so cancelling can stop it and anything it started.
  command <- paste(shQuote(Sys.which("nohup")), shQuote(Sys.which("setsid")), shQuote(Sys.which("Rscript")), shQuote(runner), "drive",
                   shQuote(hub_spec_path(p)), "</dev/null", paste0(">>", shQuote(out)), paste0("2>>", shQuote(err)), "& echo $!")
  res <- if (is.function(starter)) starter(command) else
    suppressWarnings(tryCatch(system2("bash", c("-c", shQuote(command)), stdout = TRUE, stderr = TRUE), error = function(e) conditionMessage(e)))
  pid <- trimws(res[length(res)] %||% "")
  if (!grepl("^[0-9]+$", pid)) return(list(ok = FALSE, errors = paste("Couldn't start the run controller:", paste(res, collapse = " "))))
  host <- Sys.info()[["nodename"]]
  hub_update_status(p, controller_pid = pid, controller_host = host, state = "queued", error = NULL)
  save_job(p, "Hub controller", c("codespring-run", "drive", hub_spec_path(p)),
           paste(c(paste("controller_pid:", pid), paste("host:", host), paste("stdout:", out), paste("stderr:", err)), collapse = "\n"))
  list(ok = TRUE, pid = pid, host = host)
}

hub_controller_alive <- function(st, spec_path) {
  # TRUE only if the recorded controller is still this run's controller: right host, a numeric pid,
  # and that process's command line names this run's spec file. A reused pid is never mistaken
  # for it (and nothing from the status file reaches a shell).
  pid <- as.character(st$controller_pid %||% "")
  if (!grepl("^[0-9]+$", pid) || !identical(st$controller_host, Sys.info()[["nodename"]])) return(FALSE)
  cmd <- suppressWarnings(tryCatch(readBin(file.path("/proc", pid, "cmdline"), "raw", 65536L), error = function(e) raw(0)))
  if (!length(cmd)) return(FALSE)
  cmd[cmd == as.raw(0)] <- as.raw(32)
  grepl(spec_path, rawToChar(cmd), fixed = TRUE)
}

hub_project <- function(project_id) {
  p <- discover_projects()[[project_id]]  # only the person's own CodeSpringApp projects
  if (is.null(p) || !file.exists(hub_status_path(p))) stop(hub_error("No run started from the Hub called ", project_id, " in your CodeSpringApp projects."))
  p
}

hub_relaunch <- function(p, raw, starter = NULL) {
  if (hub_controller_alive(raw, hub_spec_path(p))) return(list(ok = FALSE, errors = "This run's controller is still running; wait for it to stop."))
  raw$state <- "queued"
  raw$error <- NULL
  raw$paused_after <- NULL
  raw$next_step <- NULL
  hub_write_json(raw, hub_status_path(p))
  unlink(hub_cancel_flag(p))
  launched <- hub_launch_controller(p, starter)
  if (!isTRUE(launched$ok)) return(launched)
  list(ok = TRUE, project_id = p$id, controller_pid = launched$pid, controller_host = launched$host)
}

hub_resume <- function(project_id, starter = NULL) {
  st <- hub_status(project_id)
  if (!isTRUE(st$ok)) return(st)
  if (!st$state %in% c("failed", "stopped", "cancelled")) {
    return(list(ok = FALSE, errors = paste0("This run is ", st$state, "; only a failed, stopped or cancelled run can be resumed",
      if (identical(st$state, "paused")) " (a paused run is continued instead)" else "", ".")))
  }
  p <- hub_project(project_id)
  raw <- hub_read_json(hub_status_path(p))
  for (i in seq_along(raw$steps)) if (!raw$steps[[i]]$state %in% c("complete", "skipped")) raw$steps[[i]]$state <- "waiting"
  res <- hub_relaunch(p, raw, starter)
  if (isTRUE(res$ok)) res$resumed <- TRUE
  res
}

# ---- continuing a paused run, adding or skipping steps ------------------------------------------

hub_continue <- function(project_id, add = character(0), skip = character(0), mode = NULL, starter = NULL) {
  st <- hub_status(project_id)
  if (!isTRUE(st$ok)) return(st)
  errors <- character(0)
  adder <- function(...) errors <<- c(errors, paste0(...))
  add <- hub_optional_ids(add, adder)
  skip <- unique(as.character(unlist(skip %||% character(0))))
  if (!st$state %in% c("paused", "complete")) {
    adder("This run is ", st$state, ". Only a paused run can be continued, or a finished one given an optional step",
          if (st$state %in% c("failed", "stopped", "cancelled")) "; resume it instead" else "", ".")
  }
  if (identical(st$state, "complete") && !length(add)) adder("This run has finished; say which optional step to add (", paste(names(HUB_OPTIONAL), collapse = ", "), ").")
  if (!is.null(mode) && !mode %in% names(HUB_MODES)) adder("mode must be 'auto' or 'step'.")
  if (length(errors)) return(list(ok = FALSE, errors = errors))

  p <- hub_project(project_id)
  if (length(add) && !rna_optional_quantifiers_available(p)) return(list(ok = FALSE, errors = "Kallisto and RSEM aren't available for this reference."))
  raw <- hub_read_json(hub_status_path(p))
  spec <- hub_read_json(hub_spec_path(p))
  have <- vapply(raw$steps, `[[`, "", "step")
  waiting <- have[vapply(raw$steps, function(s) identical(s$state, "waiting"), logical(1))]
  bad_skip <- skip[!skip %in% HUB_SKIPPABLE | !skip %in% waiting]
  if (length(bad_skip)) {
    return(list(ok = FALSE, errors = paste0("Only steps that haven't run yet and aren't needed later can be skipped (",
      paste(intersect(HUB_SKIPPABLE, waiting), collapse = ", "), "); not: ", paste(bad_skip, collapse = ", "), ".")))
  }
  for (id in add) {
    step <- HUB_OPTIONAL[[id]]$step
    if (step %in% have) next
    raw$steps[[length(raw$steps) + 1]] <- list(step = step, state = "waiting", what = paste("Optional:", HUB_OPTIONAL[[id]]$what), added = hub_now())
  }
  raw$steps <- raw$steps[order(match(vapply(raw$steps, `[[`, "", "step"), HUB_STEP_ORDER))]
  for (i in seq_along(raw$steps)) if (raw$steps[[i]]$step %in% skip) raw$steps[[i]]$state <- "skipped"
  spec$optional <- as.list(union(as.character(unlist(spec$optional)), add))
  if (!is.null(mode)) { spec$mode <- mode; raw$mode <- mode }
  hub_write_json(spec, hub_spec_path(p))
  if (!any(vapply(raw$steps, function(s) identical(s$state, "waiting"), logical(1)))) {
    raw$state <- "complete"
    hub_write_json(raw, hub_status_path(p))
    return(list(ok = TRUE, project_id = p$id, note = "Nothing left to run."))
  }
  res <- hub_relaunch(p, raw, starter)
  if (isTRUE(res$ok)) res$steps <- lapply(raw$steps, function(s) list(step = s$step, state = s$state))
  res
}

# ---- cancelling ----------------------------------------------------------------------------------

hub_cancel_jobs <- function(p, steps) {
  # Cancel the queued or running jobs recorded in THIS project's job history, step by step. (No
  # job-name sweep: names of projects like 'a.b' and 'a_b' share a prefix, so a sweep could hit
  # another run.) scancel can only ever touch the person's own jobs.
  msgs <- character(0)
  for (step in steps) {
    m <- tryCatch(cancel_active_step_jobs(p, hub_app_step(step)), error = function(e) paste("ERROR:", conditionMessage(e)))
    if (grepl("Requested cancellation|ERROR", m)) msgs <- c(msgs, m)
  }
  msgs
}

hub_mark_cancelled <- function(p) {
  raw <- hub_read_json(hub_status_path(p))
  raw$controller_pid <- NULL
  for (i in seq_along(raw$steps)) if (identical(raw$steps[[i]]$state, "running")) {
    raw$steps[[i]]$state <- "cancelled"
    raw$steps[[i]]$finished <- hub_now()
  }
  raw$state <- "cancelled"
  raw$cancelled <- hub_now()
  raw$paused_after <- NULL
  raw$next_step <- NULL
  raw$updated <- hub_now()
  hub_write_json(raw, hub_status_path(p))
}

hub_cancel <- function(project_id) {
  st <- hub_status(project_id)
  if (!isTRUE(st$ok)) return(st)
  if (st$state %in% c("complete", "cancelled")) return(list(ok = FALSE, errors = paste0("This run is already ", st$state, "; there's nothing to cancel.")))
  if (identical(st$state, "cancelling")) return(list(ok = FALSE, errors = "This run is already being cancelled."))
  p <- hub_project(project_id)
  file.create(hub_cancel_flag(p))  # a controller on another login node stops at its next check
  raw <- hub_read_json(hub_status_path(p))
  stopped <- "not running"
  steps <- vapply(st$steps, `[[`, "", "step")
  if (st$state %in% c("queued", "running") && !identical(raw$controller_host, Sys.info()[["nodename"]]) && !is.null(raw$controller_pid)) {
    # Can't signal another node: its controller cancels its jobs and marks the run cancelled itself.
    msgs <- hub_cancel_jobs(p, steps)
    hub_update_status(p, state = "cancelling")
    return(list(ok = TRUE, project_id = p$id, state = "cancelling",
                controller = paste("asked to stop (it runs on", raw$controller_host, "and checks every", HUB_POLL_SECONDS, "seconds)"),
                jobs = if (length(msgs)) msgs else "No queued or running jobs were found.",
                note = "Finished steps and their results are kept. The run can be resumed once it shows as cancelled."))
  }
  if (hub_controller_alive(raw, hub_spec_path(p))) {
    # The controller was started with setsid, so it leads its own process group: stop the group.
    tools::pskill(-as.integer(raw$controller_pid), tools::SIGTERM)
    stopped <- "stopped"
  }
  msgs <- hub_cancel_jobs(p, steps)
  hub_mark_cancelled(p)
  hub_notify(p, hub_read_json(hub_spec_path(p)), "cancelled", hub_status_text(p))
  list(ok = TRUE, project_id = p$id, state = "cancelled", controller = stopped,
       jobs = if (length(msgs)) msgs else "No queued or running jobs were found.",
       note = "Finished steps and their results are kept. The run can be resumed later.")
}

# ---- driving a run (the controller) --------------------------------------------------------------

hub_submit_limit <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    env <- suppressWarnings(as.integer(Sys.getenv("CSL_HUB_MAX_SUBMIT", unset = "")))
    if (!is.na(env) && env > 0) return(cached <<- env)
    limit <- NA_integer_
    if (nzchar(Sys.which("sacctmgr"))) {
      assoc <- suppressWarnings(tryCatch(system2("sacctmgr", c("-nP", "show", "assoc", paste0("user=", CURRENT_USER), "format=DefaultQOS,MaxSubmit"),
                                                 stdout = TRUE, stderr = FALSE), error = function(e) character(0)))
      parts <- strsplit(assoc[nzchar(assoc)][1] %||% "", "|", fixed = TRUE)[[1]]
      own <- suppressWarnings(as.integer(parts[2] %||% NA))
      if (!is.na(own)) limit <- own
      qos <- trimws(parts[1] %||% "")
      if (nzchar(qos) && grepl("^[A-Za-z0-9_.-]+$", qos)) {
        q <- suppressWarnings(tryCatch(system2("sacctmgr", c("-nP", "show", "qos", qos, "format=MaxSubmitPU"), stdout = TRUE, stderr = FALSE), error = function(e) character(0)))
        qlim <- suppressWarnings(as.integer(trimws(q[nzchar(q)][1] %||% NA)))
        if (!is.na(qlim)) limit <- if (is.na(limit)) qlim else min(limit, qlim)
      }
    }
    cached <<- if (is.na(limit) || limit < 1) 16L else limit
  }
})

HUB_SLOT_WAIT_HOURS <- as.numeric(Sys.getenv("CSL_HUB_SLOT_WAIT_HOURS", unset = "12"))

hub_free_slots <- function() {
  # How many more jobs may be submitted right now without hitting the person's Slurm limit.
  # NA when squeue can't be asked: then nothing is submitted (fail closed).
  if (!nzchar(Sys.which("squeue"))) return(NA_integer_)
  queued <- suppressWarnings(tryCatch(system2("squeue", c("-h", "-u", CURRENT_USER, "-o", "%i"), stdout = TRUE, stderr = FALSE),
                                      error = function(e) NULL))
  if (is.null(queued) || !is.null(attr(queued, "status"))) return(NA_integer_)
  limit <- hub_submit_limit()
  reserve <- min(HUB_SLOT_RESERVE, max(0L, limit - 1L))  # a tiny limit still leaves room for one job
  limit - reserve - length(queued[nzchar(queued)])
}

hub_step_entries <- function(p, spec) {
  trim <- isTRUE(hub_flag(spec$trim$enabled, TRUE))
  a1 <- adapter_choices_r1()[[spec$trim$adapter_r1 %||% "Illumina Universal TruSeq RNA"]]
  a2 <- adapter_choices_r2()[[spec$trim$adapter_r2 %||% spec$trim$adapter_r1 %||% "Illumina Universal TruSeq RNA"]]
  min_len <- as.integer(spec$trim$min_length %||% 20L)
  feature <- spec$feature %||% "gene_name"
  de <- spec$deseq2
  gsea <- spec$gsea
  deg_file <- if (!is.null(de)) file.path(p$data_dir, "deseq2", sprintf("DEG_%s_vs_%s(ref).txt", de$comparison, de$reference)) else ""
  matrix_ready <- function(name) { f <- file.path(p$data_dir, "counts", name); file.exists(f) && file.size(f) > 0 }
  all_samples <- function() project_samples(p)
  list(
    Cutadapt = list(run = function(s) submit_cutadapt_jobs(p, a1, a2, min_len, s),
                    done = function() setequal(completed_samples_for_step(p, "Cutadapt", all_samples()), all_samples())),
    FastQC = list(run = function(s) submit_fastqc_jobs(p, trimmed = trim, samples = s),
                  done = function() setequal(completed_samples_for_step(p, "FastQC", all_samples()), all_samples())),
    STAR = list(run = function(s) submit_star_jobs(p, trimmed = trim, samples = s),
                done = function() { f <- expected_star_bam_files(p); length(f) && all(file.exists(f) & file.size(f) > 0) }),
    featureCounts = list(run = function(s) submit_featurecounts_jobs(p, feature, s),
                         done = function() featurecounts_outputs_ready(p) && matrix_ready("count_matrix.txt")),
    Kallisto = list(run = function(s) submit_kallisto_jobs(p, trimmed = trim, samples = s),
                    done = function() matrix_ready("kallisto_tpm_matrix.txt")),
    RSEM = list(run = function(s) submit_rsem_jobs(p, "gene_id", s),
                done = function() matrix_ready("rsem_tpm_matrix.txt")),
    DESeq2 = list(run = function(s) submit_deseq2_job(p, de$column, de$reference, de$comparison),
                  done = function() file.exists(deg_file) && file.size(deg_file) > 0),
    GSEA = list(run = function(s) submit_gseapy_job(p, de$column, de$reference, de$comparison, gsea$geneset),
                done = function() length(list.files(file.path(p$data_dir, "gseapy"), pattern = "^report\\.gseapy\\..*\\.csv$", recursive = TRUE)) > 0)
  )
}

hub_submit <- function(run, samples) {
  msg <- tryCatch(paste(run(samples), collapse = "\n"), error = function(e) paste("ERROR:", conditionMessage(e)))
  list(ok = !grepl("ERROR", msg, fixed = TRUE), message = msg)
}

hub_wait_idle <- function(p, step, spec) {
  # Wait until this step has no queued or running Slurm jobs. If the job list can't be read, keep
  # waiting (never assume idle), but give up after a while rather than loop forever.
  Sys.sleep(min(20, HUB_POLL_SECONDS))
  errors <- 0L
  repeat {
    if (hub_cancel_requested(p)) return(list(cancelled = TRUE))
    active <- tryCatch(active_jobs_for_step(p, hub_app_step(step)), error = function(e) NULL)
    if (is.null(active)) {
      errors <- errors + 1L
      if (errors * HUB_POLL_SECONDS > 3600) return(list(ok = FALSE, message = paste("Couldn't read the Slurm job list for", step, "for an hour; stopped. Resume when the cluster is back.")))
    } else {
      errors <- 0L
      if (!NROW(active)) return(list(ok = TRUE))
    }
    hub_tick(p, spec)
    Sys.sleep(HUB_POLL_SECONDS)
  }
}

hub_submit_waves <- function(p, name, entry, spec, pending) {
  # Waves: submit as many samples as the person's Slurm limit allows, then more as slots free up.
  total <- length(pending)
  waited_since <- Sys.time()
  while (length(pending)) {
    if (hub_cancel_requested(p)) return(list(cancelled = TRUE))
    free <- hub_free_slots()
    if (is.na(free) || free < 1) {
      if (difftime(Sys.time(), waited_since, units = "hours") > HUB_SLOT_WAIT_HOURS) {
        return(list(ok = FALSE, message = paste0(name, ": waited ", HUB_SLOT_WAIT_HOURS, " h for free Slurm job slots (",
          if (is.na(free)) "squeue didn't answer" else "your other jobs are using them", "). ", total - length(pending), " of ", total,
          " samples were submitted. Resume when slots are free.")))
      }
      hub_set_step(p, name, "running", paste0("Waiting for free Slurm job slots (", total - length(pending), " of ", total, " samples submitted)."))
      hub_tick(p, spec)
      Sys.sleep(HUB_POLL_SECONDS)
      next
    }
    batch <- head(pending, free)
    sub <- hub_submit(entry$run, batch)
    if (!sub$ok) return(list(ok = FALSE, message = paste(name, "was not submitted:", sub$message)))
    pending <- setdiff(pending, batch)
    waited_since <- Sys.time()
    hub_set_step(p, name, "running", if (length(pending)) paste0("Submitted ", total - length(pending), " of ", total, " samples; the rest follow as Slurm slots free up.") else NULL)
    if (length(pending)) Sys.sleep(min(10, HUB_POLL_SECONDS))
  }
  list(ok = TRUE)
}

hub_run_step <- function(p, name, entry, spec) {
  per_sample <- name %in% HUB_PER_SAMPLE_STEPS
  first <- if (per_sample) hub_submit_waves(p, name, entry, spec, project_samples(p)) else {
    sub <- hub_submit(entry$run, NULL)
    if (sub$ok) list(ok = TRUE) else list(ok = FALSE, message = paste(name, "was not submitted:", sub$message))
  }
  if (!isTRUE(first$ok)) return(first)
  waited <- hub_wait_idle(p, name, spec)
  if (!isTRUE(waited$ok)) return(waited)
  invisible(job_history(p, force_refresh = TRUE))
  if (!isTRUE(entry$done()) && hub_step_counts(p, name)$failed == 0) {
    # Every sample finished but a matrix build was held back while earlier waves were still running:
    # one more pass submits just that (the app skips samples whose outputs are complete). It goes
    # through the wave loop too, so it can never exceed the Slurm limit.
    again <- if (per_sample) hub_submit_waves(p, name, entry, spec, project_samples(p)) else {
      sub <- hub_submit(entry$run, NULL)
      if (sub$ok) list(ok = TRUE) else list(ok = FALSE, message = paste(name, "was not finished:", sub$message))
    }
    if (!isTRUE(again$ok)) return(again)
    waited <- hub_wait_idle(p, name, spec)
    if (!isTRUE(waited$ok)) return(waited)
    invisible(job_history(p, force_refresh = TRUE))
  }
  if (isTRUE(entry$done())) return(list(ok = TRUE))
  log_dir <- file.path(project_result_dir(p), "log")
  errs <- list.files(log_dir, pattern = paste0("^error_", tolower(clean_name(name, "step"))), full.names = TRUE, ignore.case = TRUE)
  failed <- hub_step_counts(p, name)
  list(ok = FALSE, message = paste0(name, " finished without all of its expected outputs",
    if (failed$failed > 0) paste0(" (", failed$failed, " job(s) failed: ", paste(head(failed$failed_samples, 5), collapse = ", "), ")") else "",
    ". Check: ", paste(head(errs, 3), collapse = ", ")))
}

hub_drive <- function(spec_path) {
  spec <- hub_read_json(spec_path)
  status_path <- file.path(dirname(spec_path), "codespring_hub_run.json")
  st <- hub_read_json(status_path)
  p <- hub_project(st$project_id)
  if (hub_cancel_requested(p)) { hub_mark_cancelled(p); return(invisible(FALSE)) }
  # Re-check the stored choices (the file sits in the person's results folder and could be edited).
  bad <- character(0)
  spec$notify <- hub_check_notify(spec$notify, function(...) bad <<- c(bad, paste0(...)))
  if (!identical(spec$mode, "step")) spec$mode <- "auto"
  spec$optional <- as.list(hub_optional_ids(spec$optional, function(...) bad <<- c(bad, paste0(...))))
  if (length(bad)) {
    hub_update_status(p, state = "failed", error = paste("The stored request was changed:", paste(bad, collapse = " ")), controller_pid = NULL)
    return(invisible(FALSE))
  }
  hub_update_status(p, state = "running", controller_host = Sys.info()[["nodename"]], paused_after = NULL, next_step = NULL)
  entries <- hub_step_entries(p, spec)
  names_all <- vapply(st$steps, `[[`, "", "step")
  for (i in seq_along(st$steps)) {
    name <- names_all[[i]]
    if (st$steps[[i]]$state %in% c("complete", "skipped")) next  # resuming: keep finished steps
    hub_set_step(p, name, "running")
    res <- hub_run_step(p, name, entries[[name]], spec)
    if (isTRUE(res$cancelled)) {
      hub_cancel_jobs(p, name)
      hub_mark_cancelled(p)
      return(invisible(FALSE))
    }
    if (!isTRUE(res$ok)) {
      hub_set_step(p, name, "failed", res$message)
      hub_update_status(p, state = "failed", error = res$message, controller_pid = NULL)
      hub_notify(p, spec, paste(name, "failed"), hub_status_text(p))
      return(invisible(FALSE))
    }
    hub_set_step(p, name, "complete")
    later <- names_all[seq_along(names_all) > i]
    later <- later[vapply(st$steps[seq_along(names_all) > i], function(s) identical(s$state, "waiting"), logical(1))]
    if (identical(spec$mode, "step") && length(later)) {
      hub_update_status(p, state = "paused", paused_after = name, next_step = later[[1]], controller_pid = NULL)
      hub_notify(p, spec, paste(name, "finished; waiting for you"), hub_status_text(p))
      return(invisible(TRUE))
    }
    if (length(later)) hub_notify(p, spec, paste(name, "finished"), hub_status_text(p))
  }
  # The post-run packet: one HTML page with the app's plots and small result files.
  report <- tryCatch(hub_report(p$id)$report, error = function(e) NULL)
  hub_update_status(p, state = "complete", finished = hub_now(), controller_pid = NULL)
  hub_notify(p, spec, "finished", c(hub_status_text(p), if (!is.null(report)) c("", paste("Report (open in a browser):", report))))
  invisible(TRUE)
}

# ---- what the person sees: status, QC, emails ---------------------------------------------------

hub_step_counts <- function(p, step, jobs = NULL) {
  # Per sample (or per job for single-job steps), the latest job's state for this step.
  jobs <- jobs %||% tryCatch(job_history(p), error = function(e) data.frame())
  out <- list(queued = 0L, running = 0L, done = 0L, failed = 0L, failed_samples = character(0))
  if (!NROW(jobs) || !all(c("step", "slurm_state") %in% names(jobs))) return(out)
  hit <- jobs[canonical_job_step(jobs$step) == canonical_job_step(hub_app_step(step)), , drop = FALSE]
  if (!NROW(hit)) return(out)
  key <- if ("sample" %in% names(hit)) ifelse(nzchar(hit$sample %||% ""), hit$sample, hit$job_id %||% seq_len(NROW(hit))) else seq_len(NROW(hit))
  latest <- hit[!duplicated(key, fromLast = TRUE), , drop = FALSE]
  state <- toupper(trimws(as.character(latest$slurm_state)))
  out$queued <- sum(state %in% c("PENDING", "PD", "SUBMITTED", "QUEUED"))
  out$running <- max(0L, sum(state %in% toupper(active_slurm_states())) - out$queued)
  out$done <- sum(state %in% HUB_DONE_STATES)
  bad <- state %in% HUB_FAILED_STATES & !grepl("cancelled_by_codespringweb", as.character(latest$output %||% ""), fixed = TRUE)
  out$failed <- sum(bad)
  if ("sample" %in% names(latest)) out$failed_samples <- as.character(latest$sample[bad & nzchar(latest$sample)])
  out
}

hub_star_stats <- function(p, sample) {
  lines <- read_metric_lines(file.path(p$data_dir, "star", sample, paste0(sample, "Log.final.out")))
  if (!length(lines)) return(NULL)
  value <- function(label) {
    hit <- grep(label, lines, value = TRUE, fixed = TRUE)
    if (!length(hit)) return(NA_real_)
    clean_metric_number(tail(strsplit(hit[[1]], "|", fixed = TRUE)[[1]], 1))
  }
  list(input_reads = value("Number of input reads"), uniquely_mapped_pct = value("Uniquely mapped reads %"),
       multi_mapped_pct = value("% of reads mapped to multiple loci"), unmapped_too_short_pct = value("% of reads unmapped: too short"),
       unmapped_other_pct = value("% of reads unmapped: other"))
}

hub_fastqc_summary <- function(p, sample, trimmed) {
  # FastQC's own PASS/WARN/FAIL per module, read from summary.txt inside each report zip.
  dir <- file.path(p$data_dir, if (trimmed) "fastqc_cutadapt" else "fastqc")
  pairs <- tryCatch(sample_fastq_pairs(p, trimmed), error = function(e) NULL)
  if (is.null(pairs) || !sample %in% pairs$sample) return(NULL)
  row <- pairs[pairs$sample == sample, , drop = FALSE][1, ]
  reads <- unique(unlist(lapply(c(row$r1, if (isTRUE(p$paired_end)) row$r2), split_fastq_path_list)))
  zips <- file.path(dir, sub(fastq_suffix_regex, "_fastqc.zip", basename(reads), ignore.case = TRUE))
  zips <- zips[file.exists(zips)]
  if (!length(zips)) return(NULL)
  rows <- unlist(lapply(zips, function(z) {
    inner <- paste0(sub("\\.zip$", "", basename(z)), "/summary.txt")
    tryCatch({ con <- unz(z, inner); on.exit(close(con)); readLines(con, warn = FALSE) }, error = function(e) character(0))
  }))
  parts <- strsplit(rows[nzchar(rows)], "\t", fixed = TRUE)
  flag <- vapply(parts, `[`, "", 1)
  module <- vapply(parts, `[`, "", 2)
  list(files = length(zips), pass = sum(flag == "PASS"), warn = sum(flag == "WARN"), fail = sum(flag == "FAIL"),
       failed_modules = I(unique(module[flag == "FAIL"])), warned_modules = I(unique(module[flag == "WARN"])))
}

hub_kallisto_stats <- function(p, sample) {
  f <- file.path(p$data_dir, "kallisto", sample, "run_info.json")
  if (!file.exists(f)) return(NULL)
  x <- tryCatch(hub_read_json(f), error = function(e) NULL)
  if (is.null(x)) return(NULL)
  list(reads = x$n_processed, pseudoaligned_pct = x$p_pseudoaligned)
}

hub_qc <- function(project_id, step = NULL) {
  p <- hub_project(project_id)
  spec <- hub_read_json(hub_spec_path(p))
  trim <- isTRUE(hub_flag(spec$trim$enabled, TRUE))
  jobs <- tryCatch(job_history(p), error = function(e) data.frame())
  want <- function(s) is.null(step) || identical(step, s)
  rows <- lapply(project_samples(p), function(sample) {
    r <- list(sample = sample)
    if (want("Cutadapt") && trim) { m <- tryCatch(extract_cutadapt_metrics(p, sample, jobs), error = function(e) NULL); if (length(m) && any(nzchar(m))) r$cutadapt <- as.list(m) }
    if (want("FastQC")) r$fastqc <- hub_fastqc_summary(p, sample, trim)
    if (want("STAR")) r$star <- hub_star_stats(p, sample)
    if (want("featureCounts")) { m <- tryCatch(extract_featurecounts_metrics(p, sample), error = function(e) NULL); if (length(m) && any(nzchar(m))) r$featurecounts <- as.list(m) }
    if (want("Kallisto")) r$kallisto <- hub_kallisto_stats(p, sample)
    r[!vapply(r, is.null, logical(1))]
  })
  list(ok = TRUE, project_id = p$id, samples = rows, summary = hub_qc_summary(rows),
       note = paste("Numbers come straight from each tool's own report. FastQC's PASS/WARN/FAIL are FastQC's judgements;",
                    "whether a sample is good enough is the researcher's call (ask BSR if unsure)."))
}

hub_qc_summary <- function(rows) {
  # A few lines per tool for a chat or an email: ranges, and which samples FastQC failed.
  pick <- function(tool, field) suppressWarnings(as.numeric(unlist(lapply(rows, function(r) r[[tool]][[field]] %||% NA))))
  rng <- function(x, unit = "") { x <- x[is.finite(x)]; if (!length(x)) return(NULL); paste0(format(round(min(x), 1), big.mark = ","), unit, " to ", format(round(max(x), 1), big.mark = ","), unit) }
  out <- list()
  fq <- Filter(Negate(is.null), lapply(rows, function(r) if (!is.null(r$fastqc)) list(sample = r$sample, fail = r$fastqc$fail, modules = r$fastqc$failed_modules)))
  if (length(fq)) out$fastqc <- list(samples_reported = length(fq),
    samples_with_a_fail = vapply(Filter(function(x) x$fail > 0, fq), function(x) paste0(x$sample, " (", paste(x$modules, collapse = ", "), ")"), ""))
  u <- rng(pick("star", "uniquely_mapped_pct"), "%")
  if (!is.null(u)) out$star <- list(uniquely_mapped = u, input_reads = rng(pick("star", "input_reads")), multi_mapped = rng(pick("star", "multi_mapped_pct"), "%"))
  k <- rng(pick("kallisto", "pseudoaligned_pct"), "%")
  if (!is.null(k)) out$kallisto <- list(pseudoaligned = k)
  out
}

hub_status_text <- function(p) {
  st <- tryCatch(hub_status(p$id), error = function(e) list(state = "unknown", steps = list()))
  lines <- c(paste0("Run ", p$name, ": ", st$state), "")
  for (s in st$steps) {
    j <- s$jobs
    extra <- if (!is.null(j) && (j$queued + j$running + j$done + j$failed) > 0) paste0("  (", j$done, " done, ", j$running, " running, ", j$queued, " queued", if (j$failed) paste0(", ", j$failed, " failed") else "", ")") else ""
    lines <- c(lines, paste0("- ", s$step, ": ", s$state, extra))
  }
  if (identical(st$state, "paused")) lines <- c(lines, "", paste0("Paused after ", st$paused_after, ". Next: ", st$next_step, ". Tell Claude to continue, add an optional step, or cancel."))
  if (!is.null(st$error)) lines <- c(lines, "", paste("Problem:", st$error))
  qc <- tryCatch(hub_qc(p$id)$summary, error = function(e) NULL)
  if (length(qc$star)) lines <- c(lines, "", paste("STAR uniquely mapped:", qc$star$uniquely_mapped))
  if (length(qc$fastqc$samples_with_a_fail)) lines <- c(lines, paste("FastQC FAIL in:", paste(head(qc$fastqc$samples_with_a_fail, 10), collapse = "; ")))
  c(lines, "", paste("Results:", project_result_dir(p)), "Ask Claude about this run, or open it in CodeSpringApp.")
}

hub_notify <- function(p, spec, subject, body) {
  # Only if the person asked for email when planning, and only to their own address.
  if (!isTRUE(spec$notify$email)) return(invisible(FALSE))
  mailer <- Sys.getenv("CSL_HUB_MAILER", unset = Sys.which("mail"))
  if (!nzchar(mailer)) return(invisible(FALSE))
  subj <- paste0("[CodeSpring] ", p$name, ": ", subject)
  ok <- tryCatch(system2(mailer, c("-s", shQuote(subj), hub_email_address()), input = body, stdout = FALSE, stderr = FALSE), error = function(e) 1L)
  hub_update_status(p, last_email = hub_now())
  invisible(identical(ok, 0L))
}

hub_tick <- function(p, spec) {
  # Called while waiting: sends the periodic summary when one is due.
  every <- suppressWarnings(as.numeric(spec$notify$every_hours %||% 0))
  if (!isTRUE(spec$notify$email) || !is.finite(every) || every <= 0) return(invisible(NULL))
  st <- hub_read_json(hub_status_path(p))
  last <- suppressWarnings(as.POSIXct(st$last_summary %||% st$created %||% hub_now(), format = "%Y-%m-%dT%H:%M:%S%z"))
  if (is.na(last) || difftime(Sys.time(), last, units = "hours") >= every) {
    hub_update_status(p, last_summary = hub_now())
    hub_notify(p, spec, "progress", hub_status_text(p))
  }
}

hub_status <- function(project_id) {
  p <- discover_projects()[[project_id]]
  if (is.null(p)) return(list(ok = FALSE, errors = paste("No CodeSpringApp project called", project_id)))
  path <- hub_status_path(p)
  st <- if (file.exists(path)) hub_read_json(path) else list(state = "not started from the Hub")
  if (identical(st$state, "cancelling")) {
    # A controller on another node turns this into "cancelled" itself; if it has gone quiet, do it here.
    last <- suppressWarnings(as.POSIXct(st$updated %||% "", format = "%Y-%m-%dT%H:%M:%S%z"))
    if (is.na(last) || difftime(Sys.time(), last, units = "secs") > 5 * HUB_POLL_SECONDS + 120) {
      hub_mark_cancelled(p)
      st <- hub_read_json(path)
    }
  }
  if (!is.null(st$controller_pid) && st$state %in% c("queued", "running")) {
    here <- identical(st$controller_host, Sys.info()[["nodename"]])
    alive <- here && hub_controller_alive(st, hub_spec_path(p))
    st$controller_state <- if (!here) paste("running on", st$controller_host, "(can't check from here)") else if (alive) "running" else "ended"
    if (identical(st$controller_state, "ended")) {
      st$state <- "stopped"
      st$error <- paste0("The run controller stopped before the run finished (for example the login node restarted). ",
                         "Finished steps are kept; start it again with `codespring-run resume ", p$id, "`. Details: ",
                         file.path(project_result_dir(p), "log", "error_hub_controller.txt"))
    }
  }
  jobs <- tryCatch(job_history(p), error = function(e) data.frame())
  for (i in seq_along(st$steps)) if (!st$steps[[i]]$state %in% c("waiting", "skipped")) {
    st$steps[[i]]$jobs <- hub_step_counts(p, st$steps[[i]]$step, jobs)
  }
  st$next_actions <- switch(st$state,
    paused = c(paste("continue: run", st$next_step), "add an optional step", "skip an optional step or GSEA", "switch to running straight through", "cancel"),
    running = , queued = c("wait (ask again any time)", "look at QC so far", "cancel"),
    complete = c("open the results", "add an optional step"),
    failed = , stopped = , cancelled = c("resume (finished steps are kept)"),
    cancelling = c("wait: the run is stopping"),
    character(0))
  report <- if (exists("hub_report_path")) hub_report_path(p) else ""
  if (nzchar(report) && file.exists(report)) st$report <- report
  if (identical(st$state, "complete")) st$next_actions <- c(st$next_actions, "get the report packet")
  st$ok <- TRUE
  st$project_id <- p$id
  st$results_dir <- project_result_dir(p)
  st
}
