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
# Spec version 2 adds, all optional and backward compatible with version 1:
#   gsea.genesets  several gene-set databases in one run (gsea.geneset, one database, still works);
#   quantifiers    extra quantification after featureCounts: "rsem" and/or "kallisto";
#   run_mode       "all" (default) runs every step in order; "step" pauses after each step until
#                  `codespring-run continue PROJECT_ID`.

HUB_SPEC_VERSION <- 2L
HUB_MAX_GENESETS <- 6L
HUB_QUANTIFIERS <- c(rsem = "RSEM", kallisto = "Kallisto")
HUB_RUN_MODES <- c("all", "step")
HUB_MAX_SAMPLES <- 96L
HUB_MAX_FASTQ_DIRS <- 10L
HUB_POLL_SECONDS <- as.numeric(Sys.getenv("CSL_HUB_POLL_SECONDS", unset = "30"))
HUB_NAME_RE <- "^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$"
HUB_COLUMN_RE <- "^[A-Za-z][A-Za-z0-9_]{0,40}$"
HUB_FAILED_STATES <- c("FAILED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY", "CANCELLED", "CANCELLED+", "CA", "PREEMPTED", "BOOT_FAIL")

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
    DESeq2 = file.path(SCRIPTS_DIR, "DESeq2", "qsub_deseq2.sh"),
    GSEA = "",
    RSEM = file.path(SCRIPTS_DIR, "RSEM", paste0("qsub_RSEM_", pe, ".sh")),
    Kallisto = file.path(SCRIPTS_DIR, "Kallisto", paste0("qsub_kallisto_", pe, ".sh"))
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
      list(id = "rnaseq", name = "Bulk RNA-seq (CodeSpringLab: cutadapt, FastQC, STAR, featureCounts, DESeq2, GSEApy; optional RSEM, Kallisto)",
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
      quantifiers = as.list(HUB_QUANTIFIERS),
      run_modes = HUB_RUN_MODES,
      limits = list(max_samples = HUB_MAX_SAMPLES, max_fastq_dirs = HUB_MAX_FASTQ_DIRS, min_replicates_per_group = 2L,
                    max_genesets = HUB_MAX_GENESETS),
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

hub_gsea_sets <- function(gsea) {
  # gsea.genesets (several) or gsea.geneset (one, spec version 1), in order, without repeats.
  if (is.null(gsea)) return(character(0))
  sets <- as.character(unlist(gsea$genesets %||% gsea$geneset %||% character(0)))
  unique(trimws(sets[!is.na(sets)]))
}

hub_quantifiers <- function(spec) {
  q <- tolower(trimws(as.character(unlist(spec$quantifiers %||% character(0)))))
  unique(q[nzchar(q) & !q %in% c("none", "skip")])
}

hub_gsea_report <- function(project, de, geneset) {
  slug <- clean_name(geneset, "geneset")
  file.path(project$data_dir, "gseapy", paste0(de$comparison, "_vs_", de$reference), slug, paste0("report.gseapy.", slug, ".csv"))
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
  sets <- hub_gsea_sets(gsea)
  if (!is.null(gsea)) {
    if (is.null(de)) add("GSEA needs a DESeq2 comparison.")
    if (!length(sets)) add("gsea.genesets needs at least one gene set.")
    bad_sets <- sets[!sets %in% GSEAPY_GENESET_OPTIONS]
    if (length(bad_sets)) add("gsea.genesets must come from: ", paste(GSEAPY_GENESET_OPTIONS, collapse = ", "), " (not: ", paste(bad_sets, collapse = ", "), ").")
    if (length(sets) > HUB_MAX_GENESETS) add("At most ", HUB_MAX_GENESETS, " gene sets in one run.")
    gsea <- list(genesets = as.list(sets))
  }
  quants <- hub_quantifiers(spec)
  bad_q <- quants[!quants %in% names(HUB_QUANTIFIERS)]
  if (length(bad_q)) add("quantifiers must be rsem and/or kallisto (or none), not: ", paste(bad_q, collapse = ", "), ".")
  run_mode <- tolower(trimws(as.character(spec$run_mode %||% "all")))
  if (!run_mode %in% HUB_RUN_MODES) add("run_mode must be 'all' or 'step'.")
  if (length(errors)) return(list(ok = FALSE, errors = errors))

  project <- new_project_from_inputs(list(
    new_project_analysis = "RNA-seq", new_project_name = name, new_project_mode = "new",
    new_results_root = DEFAULT_RESULTS_ROOT, new_fastq_location_mode = "multiple",
    new_fastq_dirs = paste(fastq_dirs, collapse = "\n"), new_paired_end = if (paired) "paired" else "single",
    new_species = species, new_genome_version = version
  ))
  if (length(quants) && !rna_optional_quantifiers_available(project)) {
    return(list(ok = FALSE, errors = "RSEM and Kallisto aren't available for this reference (maize). Choose no extra quantification."))
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
       feature = feature, deseq2 = de, gsea = gsea, quantifiers = quants, run_mode = run_mode, reads = reads)
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
  if (!is.null(checked$deseq2)) steps[[length(steps) + 1]] <- step("DESeq2", 1L, paste0("Differential expression: ", checked$deseq2$comparison, " vs ", checked$deseq2$reference, " (reference) on '", checked$deseq2$column, "'"))
  if (!is.null(checked$gsea)) {
    sets <- hub_gsea_sets(checked$gsea)
    steps[[length(steps) + 1]] <- step("GSEA", length(sets), paste("Pathway enrichment with", paste(sets, collapse = ", ")))
  }
  if ("rsem" %in% checked$quantifiers) steps[[length(steps) + 1]] <- step("RSEM", n + 1L, "Gene and transcript quantification from the STAR alignments (gene_id), and TPM matrices")
  if ("kallisto" %in% checked$quantifiers) steps[[length(steps) + 1]] <- step("Kallisto", n + 1L, paste("Transcript quantification from the", if (isTRUE(checked$trim$enabled)) "trimmed" else "raw", "reads, and TPM matrices"))
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
    run_mode = checked$run_mode,
    controller = list(what = "A small background process on the login node that submits each step when the previous one has finished (it uses no Slurm resources itself)"),
    runs_as = CURRENT_USER,
    codespringapp = APP_ROOT,
    codespringlab = CSL_ROOT
  )
}

# ---- starting and driving a run ---------------------------------------------------------------

hub_status_path <- function(project) file.path(project_result_dir(project), "log", "codespring_hub_run.json")
hub_spec_path <- function(project) file.path(project_result_dir(project), "log", "codespring_hub_spec.json")

hub_write_json <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null"), tmp)
  file.rename(tmp, path)
}

hub_update_status <- function(project, ...) {
  path <- hub_status_path(project)
  st <- if (file.exists(path)) jsonlite::fromJSON(path, simplifyVector = FALSE) else list()
  changes <- list(...)
  for (k in names(changes)) st[[k]] <- changes[[k]]
  st$updated <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  hub_write_json(st, path)
  invisible(st)
}

hub_set_step <- function(project, name, state, message = NULL) {
  path <- hub_status_path(project)
  st <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  for (i in seq_along(st$steps)) if (identical(st$steps[[i]]$step, name)) {
    st$steps[[i]]$state <- state
    now <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    if (identical(state, "running")) st$steps[[i]]$started <- now
    if (state %in% c("complete", "failed")) st$steps[[i]]$finished <- now
    if (!is.null(message)) st$steps[[i]]$message <- message
  }
  st$updated <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
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
  spec$checked_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  hub_write_json(spec, hub_spec_path(p))
  hub_write_json(list(
    pipeline = "rnaseq", project_id = p$id, project_name = p$name, results_dir = project_result_dir(p),
    state = "queued", started_by = "CodeSpringHub", user = CURRENT_USER,
    steps = lapply(plan$steps, function(s) list(step = s$step, state = "waiting", what = s$what)),
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  ), hub_status_path(p))

  launched <- hub_launch_controller(p, starter)
  if (!isTRUE(launched$ok)) {
    # Nothing was started: remove only what this call just created (the folder was empty before).
    unlink(project_result_dir(p), recursive = TRUE)
    unlink(cfg)
    return(launched)
  }
  pid <- launched$pid; host <- launched$host
  list(ok = TRUE, project_id = p$id, results_dir = project_result_dir(p), controller_pid = pid, controller_host = host,
       status_file = hub_status_path(p), project_config = cfg, plan = plan)
}

hub_launch_controller <- function(p, starter = NULL) {
  log_dir <- file.path(project_result_dir(p), "log")
  out <- file.path(log_dir, "output_hub_controller.txt")
  err <- file.path(log_dir, "error_hub_controller.txt")
  runner <- file.path(APP_ROOT, "bin", "codespring-run")
  # Like CodeSpringApp's Sarek and nf-core CUT&RUN controllers, the controller is a small detached
  # process on this login node: it mostly sleeps, and it needs the app's R packages, which the
  # compute nodes don't have. setsid + nohup keep it running after Claude Desktop disconnects.
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

hub_resume <- function(project_id, starter = NULL) {
  st <- hub_status(project_id)
  if (!isTRUE(st$ok)) return(st)
  if (!st$state %in% c("failed", "stopped")) return(list(ok = FALSE, errors = paste0("This run is ", st$state, "; only a failed or stopped run can be resumed.")))
  p <- discover_projects()[[project_id]]
  path <- hub_status_path(p)
  raw <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  for (i in seq_along(raw$steps)) if (!identical(raw$steps[[i]]$state, "complete")) raw$steps[[i]]$state <- "waiting"
  hub_write_json(raw, path)
  launched <- hub_launch_controller(p, starter)
  if (!isTRUE(launched$ok)) return(launched)
  list(ok = TRUE, project_id = p$id, resumed = TRUE, controller_pid = launched$pid, controller_host = launched$host)
}

hub_continue <- function(project_id, starter = NULL) {
  # A run in one-step-at-a-time mode waits after each step until the person says to go on.
  st <- hub_status(project_id)
  if (!isTRUE(st$ok)) return(st)
  if (!identical(st$state, "paused")) return(list(ok = FALSE, errors = paste0("This run is ", st$state, "; only a paused run can be continued.")))
  p <- discover_projects()[[project_id]]
  launched <- hub_launch_controller(p, starter)
  if (!isTRUE(launched$ok)) return(launched)
  nxt <- Filter(function(s) !identical(s$state, "complete"), st$steps)
  list(ok = TRUE, project_id = p$id, continued = TRUE, next_step = if (length(nxt)) nxt[[1]]$step else NULL,
       controller_pid = launched$pid, controller_host = launched$host)
}

hub_step_progress <- function(project, step) {
  # How many samples a running per-sample step has finished, where the app can tell cheaply.
  samples <- tryCatch(project_samples(project), error = function(e) character(0))
  if (!length(samples)) return(NULL)
  done <- tryCatch(switch(step,
    Cutadapt = , FastQC = length(completed_samples_for_step(project, step, samples)),
    STAR = { f <- file.path(project$data_dir, "star", samples, paste0(samples, "Aligned.sortedByCoord.out.bam")); sum(file.exists(f) & file.size(f) > 0) },
    NULL), error = function(e) NULL)
  if (is.null(done)) return(NULL)
  list(samples_done = as.integer(done), samples_total = length(samples))
}

hub_project_from_status <- function(status_path) {
  st <- jsonlite::fromJSON(status_path, simplifyVector = FALSE)
  projects <- discover_projects()
  p <- projects[[st$project_id]]
  if (is.null(p)) stop(hub_error("Project ", st$project_id, " isn't in your CodeSpringApp projects."))
  p
}

hub_wait_step <- function(project, step, done) {
  # Wait until this step has no queued or running Slurm jobs, then judge it by its outputs.
  Sys.sleep(min(20, HUB_POLL_SECONDS))
  repeat {
    active <- tryCatch(active_jobs_for_step(project, step), error = function(e) data.frame())
    if (!NROW(active)) break
    Sys.sleep(HUB_POLL_SECONDS)
  }
  invisible(job_history(project, force_refresh = TRUE))
  if (isTRUE(done())) return(list(ok = TRUE))
  log_dir <- file.path(project_result_dir(project), "log")
  errs <- list.files(log_dir, pattern = paste0("^error_", tolower(clean_name(step, "step"))), full.names = TRUE, ignore.case = TRUE)
  list(ok = FALSE, message = paste0(step, " finished without all of its expected outputs. Check: ", paste(head(errs, 3), collapse = ", ")))
}

hub_drive <- function(spec_path) {
  spec <- jsonlite::fromJSON(spec_path, simplifyVector = FALSE)
  status_path <- file.path(dirname(spec_path), "codespring_hub_run.json")
  p <- hub_project_from_status(status_path)
  st <- jsonlite::fromJSON(status_path, simplifyVector = FALSE)
  hub_update_status(p, state = "running", controller_host = Sys.info()[["nodename"]])
  samples <- project_samples(p)
  trim <- isTRUE(hub_flag(spec$trim$enabled, TRUE))
  a1 <- adapter_choices_r1()[[spec$trim$adapter_r1 %||% "Illumina Universal TruSeq RNA"]]
  a2 <- adapter_choices_r2()[[spec$trim$adapter_r2 %||% spec$trim$adapter_r1 %||% "Illumina Universal TruSeq RNA"]]
  min_len <- as.integer(spec$trim$min_length %||% 20L)
  feature <- spec$feature %||% "gene_name"
  de <- spec$deseq2
  gsea_sets <- hub_gsea_sets(spec$gsea)
  quants <- hub_quantifiers(spec)
  step_mode <- identical(tolower(as.character(spec$run_mode %||% "all")), "step")

  deg_file <- if (!is.null(de)) file.path(p$data_dir, "deseq2", sprintf("DEG_%s_vs_%s(ref).txt", de$comparison, de$reference)) else ""
  steps <- list(
    Cutadapt = list(run = function() submit_cutadapt_jobs(p, a1, a2, min_len, samples),
                    done = function() setequal(completed_samples_for_step(p, "Cutadapt", samples), samples)),
    FastQC = list(run = function() submit_fastqc_jobs(p, trimmed = trim, samples = samples),
                  done = function() setequal(completed_samples_for_step(p, "FastQC", samples), samples)),
    STAR = list(run = function() submit_star_jobs(p, trimmed = trim, samples = samples),
                done = function() { f <- expected_star_bam_files(p); length(f) && all(file.exists(f) & file.size(f) > 0) }),
    featureCounts = list(run = function() submit_featurecounts_jobs(p, feature, samples),
                         done = function() featurecounts_outputs_ready(p) && file.exists(file.path(p$data_dir, "counts", "count_matrix.txt"))),
    DESeq2 = list(run = function() submit_deseq2_job(p, de$column, de$reference, de$comparison),
                  done = function() file.exists(deg_file) && file.size(deg_file) > 0),
    GSEA = list(run = function() {
                  # One job per gene-set database; each writes to its own gseapy/<comparison>/<set>/ folder.
                  todo <- gsea_sets[!file.exists(vapply(gsea_sets, function(g) hub_gsea_report(p, de, g), ""))]
                  paste(vapply(todo, function(g) paste(submit_gseapy_job(p, de$column, de$reference, de$comparison, g), collapse = "\n"), ""), collapse = "\n")
                },
                done = function() length(gsea_sets) > 0 && all(file.exists(vapply(gsea_sets, function(g) hub_gsea_report(p, de, g), "")))),
    RSEM = list(run = function() submit_rsem_jobs(p, "gene_id", samples),
                done = function() file_size_for(file.path(p$data_dir, "counts", "rsem_tpm_matrix.txt")) > 0),
    Kallisto = list(run = function() submit_kallisto_jobs(p, trimmed = trim, samples = samples),
                    done = function() file_size_for(file.path(p$data_dir, "counts", "kallisto_tpm_matrix.txt")) > 0)
  )
  for (i in seq_along(st$steps)) {
    s <- st$steps[[i]]
    name <- s$step
    if (identical(s$state, "complete")) next  # resuming or continuing: keep finished steps
    hub_set_step(p, name, "running")
    msg <- tryCatch(steps[[name]]$run(), error = function(e) paste("ERROR:", conditionMessage(e)))
    msg <- paste(msg, collapse = "\n")
    if (grepl("ERROR", msg, fixed = TRUE)) {
      hub_set_step(p, name, "failed", msg)
      hub_update_status(p, state = "failed", error = paste(name, "was not submitted:", msg))
      return(invisible(FALSE))
    }
    res <- hub_wait_step(p, name, steps[[name]]$done)
    if (!isTRUE(res$ok)) {
      hub_set_step(p, name, "failed", res$message)
      hub_update_status(p, state = "failed", error = res$message)
      return(invisible(FALSE))
    }
    hub_set_step(p, name, "complete")
    if (step_mode && i < length(st$steps)) {
      # One step at a time: stop here; `codespring-run continue` starts a new controller for the next step.
      hub_update_status(p, state = "paused", paused_after = name, controller_pid = NULL)
      return(invisible(TRUE))
    }
  }
  hub_update_status(p, state = "complete", paused_after = NULL, finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  invisible(TRUE)
}

hub_status <- function(project_id) {
  p <- discover_projects()[[project_id]]
  if (is.null(p)) return(list(ok = FALSE, errors = paste("No CodeSpringApp project called", project_id)))
  path <- hub_status_path(p)
  st <- if (file.exists(path)) jsonlite::fromJSON(path, simplifyVector = FALSE) else list(state = "not started from the Hub")
  if (!is.null(st$controller_pid) && st$state %in% c("queued", "running")) {
    here <- identical(st$controller_host, Sys.info()[["nodename"]])
    alive <- here && identical(suppressWarnings(system2("kill", c("-0", st$controller_pid), stdout = FALSE, stderr = FALSE)), 0L)
    st$controller_state <- if (!here) paste("running on", st$controller_host, "(can't check from here)") else if (alive) "running" else "ended"
    if (identical(st$controller_state, "ended")) {
      st$state <- "stopped"
      st$error <- paste0("The run controller stopped before the run finished (for example the login node restarted). ",
                         "Finished steps are kept; start it again with `codespring-run resume ", p$id, "`. Details: ",
                         file.path(project_result_dir(p), "log", "error_hub_controller.txt"))
    }
  }
  for (i in seq_along(st$steps)) if (identical(st$steps[[i]]$state, "running")) {
    prog <- hub_step_progress(p, st$steps[[i]]$step)
    if (!is.null(prog)) st$steps[[i]] <- c(st$steps[[i]], prog)
  }
  st$ok <- TRUE
  st$project_id <- p$id
  st$results_dir <- project_result_dir(p)
  st
}
