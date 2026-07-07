library(shiny)

DT_AVAILABLE <- requireNamespace("DT", quietly = TRUE)
BASE64_AVAILABLE <- requireNamespace("base64enc", quietly = TRUE)

table_output <- function(output_id) {
  if (DT_AVAILABLE) DT::dataTableOutput(output_id) else tableOutput(output_id)
}

render_csl_table <- function(expr, page_length = 25, editable = FALSE) {
  if (DT_AVAILABLE) {
    DT::renderDataTable({
      df <- expr
      if (!NROW(df)) df <- data.frame()
      DT::datatable(df, editable = editable, rownames = FALSE, options = list(scrollX = TRUE, pageLength = page_length))
    })
  } else {
    renderTable({
      df <- expr
      if (!NROW(df)) return(data.frame())
      utils::head(df, 500)
    }, striped = TRUE, bordered = TRUE, spacing = "s")
  }
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x)) || !nzchar(as.character(x)[1])) y else x
}

clean_name <- function(x, fallback = "sample") {
  x <- gsub("[^A-Za-z0-9_]+", "_", trimws(as.character(x)))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, fallback)
}

find_codespringlab_root <- function() {
  env_root <- Sys.getenv("CSL_CODESPRINGLAB_ROOT", unset = "")
  candidates <- unique(c(
    env_root,
    getwd(),
    dirname(getwd()),
    path.expand("~/CodeSpringLab"),
    path.expand("~/CSH/CodeSpringLab"),
    "/grid/bsr/home/rouse/CodeSpringLab",
    "/Users/rouse/CSH/CodeSpringLab"
  ))
  for (candidate in candidates[nzchar(candidates)]) {
    if (dir.exists(file.path(candidate, "scripts_DoNotTouch"))) {
      return(normalizePath(candidate, winslash = "/", mustWork = FALSE))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

CSL_ROOT <- find_codespringlab_root()
SCRIPTS_DIR <- file.path(CSL_ROOT, "scripts_DoNotTouch")
APP_HOME <- path.expand(Sys.getenv("CSL_WEB_HOME", unset = "~/.codespringweb"))
dir.create(APP_HOME, recursive = TRUE, showWarnings = FALSE)
JOBS_PATH <- file.path(APP_HOME, "jobs.tsv")

analysis_label <- function(x) {
  x <- tolower(as.character(x %||% "rna"))
  if (grepl("atac", x)) return("ATAC-seq")
  if (grepl("chip", x)) return("ChIP-seq")
  "RNA-seq"
}

analysis_key <- function(x) {
  x <- tolower(as.character(x %||% "rna"))
  if (grepl("atac", x)) return("atac")
  if (grepl("chip", x)) return("chip")
  "rna"
}

analysis_notebook_dir <- function(key) {
  switch(analysis_key(key), atac = "bulkATACseq", chip = "bulkChIPseq", rna = "bulkRNAseq")
}

parse_py_config <- function(path) {
  values <- list()
  if (!file.exists(path)) return(values)
  lines <- readLines(path, warn = FALSE)
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line) || grepl("^#", line)) next
    m <- regexec("^([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*(.*)$", line)
    hit <- regmatches(line, m)[[1]]
    if (length(hit) != 3) next
    key <- hit[2]
    val <- trimws(hit[3])
    val <- sub("\\s+#.*$", "", val)
    if ((startsWith(val, "'") && endsWith(val, "'")) || (startsWith(val, "\"") && endsWith(val, "\""))) {
      val <- substr(val, 2, nchar(val) - 1)
    }
    values[[key]] <- val
  }
  values
}

resolve_legacy_path <- function(value, key = "rna") {
  value <- trimws(as.character(value %||% ""))
  if (!nzchar(value)) return("")
  value <- path.expand(value)
  if (startsWith(value, "/")) return(normalizePath(value, winslash = "/", mustWork = FALSE))
  base <- file.path(CSL_ROOT, analysis_notebook_dir(key))
  normalizePath(file.path(base, value), winslash = "/", mustWork = FALSE)
}

with_slash <- function(path) {
  path <- trimws(as.character(path %||% ""))
  if (!nzchar(path)) return(path)
  paste0(sub("/+$", "", path), "/")
}

design_path_from_dir <- function(path) {
  path <- trimws(as.character(path %||% ""))
  if (!nzchar(path)) return("")
  if (basename(path) == "design_matrix.txt") return(path)
  file.path(path, "design_matrix.txt")
}

legacy_project_from_config <- function(path) {
  vals <- parse_py_config(path)
  if (!length(vals) && basename(path) != "config.py") return(NULL)
  key <- analysis_key(vals$analysis_type %||% basename(dirname(path)))
  project_name <- vals$project_name %||% tools::file_path_sans_ext(basename(path))
  if (!nzchar(project_name)) return(NULL)
  results_root <- resolve_legacy_path(vals$results_directory %||% "../../csl_results/", key)
  visualizer_data_dir <- resolve_legacy_path(vals$visualizer_data_dir %||% "", key)
  if (nzchar(visualizer_data_dir) && basename(visualizer_data_dir) == "data" && basename(dirname(visualizer_data_dir)) == project_name) {
    results_root <- dirname(dirname(visualizer_data_dir))
  }
  inpath_design <- resolve_legacy_path(vals$inpath_design %||% "", key)
  fastq_dir <- resolve_legacy_path(vals$read_path_destination %||% vals$read_path_original %||% "", key)
  pairing <- tolower(vals$pairing %||% "y")
  data_dir <- if (nzchar(visualizer_data_dir)) visualizer_data_dir else file.path(results_root, project_name, "data")
  list(
    id = paste0(key, "/", clean_name(project_name, "project")),
    name = clean_name(project_name, "project"),
    label = project_name,
    analysis = analysis_label(key),
    analysis_key = key,
    genome = tolower(vals$genome %||% "mouse"),
    paired_end = !(pairing %in% c("n", "no", "false", "single", "se")),
    results_root = results_root,
    data_dir = data_dir,
    fastq_dir = fastq_dir,
    design_matrix_path = design_path_from_dir(inpath_design),
    source_config = normalizePath(path, winslash = "/", mustWork = FALSE),
    source = "CodeSpringLab config"
  )
}

discover_projects <- function() {
  roots <- c(file.path(SCRIPTS_DIR, "project_configs"), file.path(CSL_ROOT, "project_configs"))
  files <- character(0)
  for (root in roots) {
    if (dir.exists(root)) files <- c(files, list.files(root, pattern = "\\.py$", recursive = TRUE, full.names = TRUE))
  }
  active <- file.path(SCRIPTS_DIR, "config.py")
  if (file.exists(active)) files <- c(files, active)
  files <- unique(normalizePath(files, winslash = "/", mustWork = FALSE))
  projects <- Filter(Negate(is.null), lapply(files, legacy_project_from_config))
  if (!length(projects)) {
    projects <- list(list(
      id = "rna/example_dataset",
      name = "example_dataset",
      label = "example_dataset",
      analysis = "RNA-seq",
      analysis_key = "rna",
      genome = "mouse",
      paired_end = TRUE,
      results_root = normalizePath(path.expand("~/csl_results"), winslash = "/", mustWork = FALSE),
      data_dir = normalizePath(path.expand("~/csl_results/example_dataset/data"), winslash = "/", mustWork = FALSE),
      fastq_dir = "",
      design_matrix_path = normalizePath(path.expand("~/csl_results/example_dataset/data/manifest/design_matrix.txt"), winslash = "/", mustWork = FALSE),
      source_config = "",
      source = "default"
    ))
  }
  names(projects) <- vapply(projects, `[[`, character(1), "id")
  projects
}

safe_read_table <- function(path, n = Inf) {
  if (!file.exists(path)) return(data.frame())
  ext <- tolower(tools::file_ext(path))
  sep <- if (ext == "csv") "," else "\t"
  tryCatch({
    utils::read.table(
      path,
      sep = sep,
      header = TRUE,
      quote = "\"",
      comment.char = "",
      check.names = FALSE,
      nrows = if (is.finite(n)) n else -1
    )
  }, error = function(e) {
    tryCatch({
      utils::read.table(
        path,
        sep = "",
        header = TRUE,
        quote = "\"",
        comment.char = "",
        check.names = FALSE,
        nrows = if (is.finite(n)) n else -1
      )
    }, error = function(e2) data.frame())
  })
}

render_data_table <- function(df, page_length = 25, height = NULL) {
  if (!NROW(df)) return(tags$div(class = "empty-box", "No rows available."))
  if (DT_AVAILABLE) {
    DT::datatable(
      df,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = page_length, scrollX = TRUE, scrollY = height %||% "520px")
    )
  } else {
    tableOutput(NULL)
  }
}

download_table <- function(df, file) {
  utils::write.csv(df, file, row.names = FALSE)
}

fastq_suffix_regex <- "\\.(fastq\\.gz|fq\\.gz|fastq|fq)$"

fastq_files <- function(folder) {
  if (!dir.exists(folder)) return(character(0))
  files <- list.files(folder, full.names = FALSE)
  files[grepl(fastq_suffix_regex, tolower(files))]
}

mate_name <- function(x, mate = 2) {
  stem <- sub(fastq_suffix_regex, "", x, ignore.case = TRUE)
  suffix <- regmatches(x, regexpr(fastq_suffix_regex, x, ignore.case = TRUE))
  if (!length(suffix) || suffix == "-1") suffix <- ""
  if (mate == 2) {
    out <- sub("([._-]R)1([._-]?[0-9]*)$", "\\12\\2", stem, ignore.case = TRUE)
    if (identical(out, stem)) out <- sub("([._-])1$", "\\12", stem)
  } else {
    out <- sub("([._-]R)2([._-]?[0-9]*)$", "\\11\\2", stem, ignore.case = TRUE)
    if (identical(out, stem)) out <- sub("([._-])2$", "\\11", stem)
  }
  if (identical(out, stem)) return(NA_character_)
  paste0(out, suffix)
}

infer_sample <- function(x) {
  stem <- sub(fastq_suffix_regex, "", basename(x), ignore.case = TRUE)
  stem <- sub("([._-]R)[12]([._-]?[0-9]*)$", "", stem, ignore.case = TRUE)
  stem <- sub("([._-])[12]$", "", stem)
  clean_name(stem)
}

scan_fastqs <- function(folder, paired = TRUE, metadata_cols = "treatment") {
  files <- fastq_files(folder)
  rows <- list()
  used <- character(0)
  if (paired) {
    for (r1 in files) {
      if (r1 %in% used) next
      r2 <- mate_name(r1, 2)
      if (!is.na(r2) && r2 %in% files) {
        rows[[length(rows) + 1]] <- data.frame(include = TRUE, sample = infer_sample(r1), filename = paste(r1, r2, sep = ","), status = "paired")
        used <- c(used, r1, r2)
      } else if (grepl("([._-]R)1|([._-])1", r1, ignore.case = TRUE)) {
        rows[[length(rows) + 1]] <- data.frame(include = FALSE, sample = infer_sample(r1), filename = r1, status = "missing R2")
        used <- c(used, r1)
      }
    }
  } else {
    for (f in files) rows[[length(rows) + 1]] <- data.frame(include = TRUE, sample = infer_sample(f), filename = f, status = "single")
  }
  df <- if (length(rows)) do.call(rbind, rows) else data.frame(include = logical(), sample = character(), filename = character(), status = character())
  for (col in metadata_cols) if (!col %in% names(df)) df[[col]] <- ""
  df[, c("include", "sample", metadata_cols, "filename", "status"), drop = FALSE]
}

write_design_matrix <- function(project, df, metadata_cols) {
  keep <- df[isTRUE(df$include) | df$include %in% c(TRUE, "TRUE", "true", "1"), , drop = FALSE]
  if (!NROW(keep)) stop("No samples are included.")
  keep$sample <- clean_name(keep$sample)
  out <- project$design_matrix_path
  if (!nzchar(out) || basename(out) != "design_matrix.txt") {
    out <- file.path(project$data_dir, "manifest", "design_matrix.txt")
  }
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  keep <- keep[, c("sample", metadata_cols, "filename"), drop = FALSE]
  utils::write.table(keep, out, sep = "\t", row.names = FALSE, quote = FALSE)
  out
}

count_files <- function(path, pattern) {
  if (!dir.exists(path)) return(0)
  length(list.files(path, pattern = pattern, recursive = TRUE, full.names = TRUE))
}

project_status <- function(project) {
  data_dir <- project$data_dir
  design <- project$design_matrix_path
  data.frame(
    step = c("Setup", "Design matrix", "FASTQ reads", "FastQC", "Cutadapt", "STAR", "Kallisto", "featureCounts", "Count matrix", "DESeq2", "GSEA"),
    status = c(
      if (nzchar(project$name)) "Complete" else "Needs attention",
      if (file.exists(design)) "Complete" else "Missing",
      if (dir.exists(project$fastq_dir) && length(fastq_files(project$fastq_dir))) "Complete" else "Optional/missing",
      if (count_files(file.path(data_dir, "fastqc"), "\\.html$") + count_files(file.path(data_dir, "fastqc_cutadapt"), "\\.html$") > 0) "Complete" else "Not found",
      if (count_files(file.path(data_dir, "cutadapt"), fastq_suffix_regex) > 0) "Complete" else "Not found",
      if (count_files(file.path(data_dir, "star"), "Aligned\\.sortedByCoord\\.out\\.bam$") > 0) "Complete" else "Not found",
      if (count_files(file.path(data_dir, "kallisto"), "abundance\\.tsv$") > 0) "Complete" else "Not found",
      if (count_files(file.path(data_dir, "featurecounts"), "_counts\\.txt$") > 0) "Complete" else "Not found",
      if (file.exists(file.path(data_dir, "counts", "count_matrix.txt"))) "Complete" else "Not found",
      if (count_files(file.path(data_dir, "deseq2"), "DEG|normalized") > 0) "Complete" else "Not found",
      if (count_files(file.path(data_dir, "gseapy"), "\\.(csv|txt|png|pdf)$") > 0) "Complete" else "Not found"
    ),
    path = c(
      dirname(data_dir),
      design,
      project$fastq_dir,
      file.path(data_dir, "fastqc"),
      file.path(data_dir, "cutadapt"),
      file.path(data_dir, "star"),
      file.path(data_dir, "kallisto"),
      file.path(data_dir, "featurecounts"),
      file.path(data_dir, "counts", "count_matrix.txt"),
      file.path(data_dir, "deseq2"),
      file.path(data_dir, "gseapy")
    ),
    stringsAsFactors = FALSE
  )
}

sample_progress <- function(project) {
  design <- safe_read_table(project$design_matrix_path)
  if (!NROW(design) || !"sample" %in% names(design)) return(data.frame())
  data_dir <- project$data_dir
  rows <- lapply(seq_len(NROW(design)), function(i) {
    sample <- as.character(design$sample[i])
    data.frame(
      sample = sample,
      FastQC = if (count_files(file.path(data_dir, "fastqc"), paste0(sample, ".*\\.html$")) > 0 || count_files(file.path(data_dir, "fastqc_cutadapt"), paste0(sample, ".*\\.html$")) > 0) "ready" else "missing",
      Trim = if (count_files(file.path(data_dir, "cutadapt"), paste0(sample, ".*", fastq_suffix_regex)) > 0) "ready" else "missing",
      STAR = if (file.exists(file.path(data_dir, "star", sample, paste0(sample, "Aligned.sortedByCoord.out.bam")))) "ready" else "missing",
      Kallisto = if (file.exists(file.path(data_dir, "kallisto", sample, "abundance.tsv"))) "ready" else "missing",
      featureCounts = if (file.exists(file.path(data_dir, "featurecounts", sample, paste0(sample, "_counts.txt")))) "ready" else "missing",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

save_job <- function(project, step, command, output = "") {
  row <- data.frame(
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    project = project$name,
    analysis = project$analysis,
    step = step,
    command = paste(command, collapse = " "),
    output = output,
    stringsAsFactors = FALSE
  )
  utils::write.table(row, JOBS_PATH, sep = "\t", row.names = FALSE, quote = TRUE, append = file.exists(JOBS_PATH), col.names = !file.exists(JOBS_PATH))
}

rna_workdir <- function(project) {
  normalizePath(file.path(CSL_ROOT, analysis_notebook_dir(project$analysis_key)), winslash = "/", mustWork = FALSE)
}

genome_resources <- function(project) {
  genome <- tolower(project$genome %||% "mouse")
  if (genome == "human") {
    list(
      star_index = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/hg38_p13_gencode_rel42_all_starindex",
      kallisto_index = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v45.transcripts.idx",
      gtf = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.chr_patch_hapl_scaff.annotation.gtf",
      strand_bed = "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.chr_patch_hapl_scaff.annotation_forStrandDetect_geneID.bed"
    )
  } else {
    list(
      star_index = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/GRCm39_M29_gencode_starindex",
      kallisto_index = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.transcripts.idx",
      gtf = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.annotation.gtf",
      strand_bed = "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.annotation_forStrandDetect_geneID.bed"
    )
  }
}

resolve_read_path <- function(base, value) {
  value <- trimws(as.character(value %||% ""))
  if (!nzchar(value)) return("")
  if (startsWith(path.expand(value), "/")) return(path.expand(value))
  file.path(base, basename(value))
}

sample_fastq_pairs <- function(project, trimmed = FALSE) {
  design <- safe_read_table(project$design_matrix_path)
  if (!NROW(design) || !"sample" %in% names(design) || !"filename" %in% names(design)) return(data.frame())
  base <- if (trimmed) file.path(project$data_dir, "cutadapt") else project$fastq_dir
  rows <- lapply(seq_len(NROW(design)), function(i) {
    parts <- trimws(unlist(strsplit(as.character(design$filename[i]), ",")))
    parts <- parts[nzchar(parts)]
    if (!length(parts)) return(NULL)
    r1 <- resolve_read_path(base, parts[1])
    r2 <- if (project$paired_end && length(parts) >= 2) resolve_read_path(base, parts[2]) else r1
    data.frame(sample = as.character(design$sample[i]), r1 = r1, r2 = r2, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(out)) data.frame() else out
}

parse_sbatch_job_id <- function(output) {
  m <- regexpr("Submitted batch job[[:space:]]+[0-9]+", output)
  if (m < 0) return("")
  sub(".*Submitted batch job[[:space:]]+", "", regmatches(output, m))
}

submit_sbatch <- function(project, step, script, args, log_name) {
  log_dir <- file.path(dirname(project$data_dir), "log")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  stdout <- file.path(log_dir, paste0("output_", log_name, "_", stamp, ".txt"))
  stderr <- file.path(log_dir, paste0("error_", log_name, "_", stamp, ".txt"))
  cmd <- c("sbatch", "-e", stderr, "-o", stdout, script, args)
  if (Sys.which("sbatch") == "") {
    msg <- "sbatch was not found. Run on the server to submit jobs."
    save_job(project, step, cmd, msg)
    return(msg)
  }
  wd <- rna_workdir(project)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  setwd(wd)
  out <- tryCatch(system2(cmd[1], cmd[-1], stdout = TRUE, stderr = TRUE), error = function(e) conditionMessage(e))
  out_text <- paste(out, collapse = "\n")
  job_id <- parse_sbatch_job_id(out_text)
  save_job(project, step, cmd, paste(c(out_text, if (nzchar(job_id)) paste("job_id:", job_id), paste("stdout:", stdout), paste("stderr:", stderr)), collapse = "\n"))
  out_text
}

submit_fastqc_jobs <- function(project, trimmed = FALSE) {
  outdir <- file.path(project$data_dir, if (trimmed) "fastqc_cutadapt" else "fastqc")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  pairs <- sample_fastq_pairs(project, trimmed)
  if (!NROW(pairs)) return("No samples/read files found in design matrix.")
  reads <- unique(c(pairs$r1, if (project$paired_end) pairs$r2 else character(0)))
  script <- file.path(SCRIPTS_DIR, "FastQC", "qsub_fastqc.sh")
  paste(vapply(reads, function(read) submit_sbatch(project, "FastQC", script, c(read, outdir, project$name), "fastQC"), character(1)), collapse = "\n")
}

submit_cutadapt_jobs <- function(project, adapter1, adapter2, min_length) {
  outdir <- file.path(project$data_dir, "cutadapt")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  pairs <- sample_fastq_pairs(project, FALSE)
  if (!NROW(pairs)) return("No samples/read files found in design matrix.")
  script <- file.path(SCRIPTS_DIR, if (project$paired_end) "cutadapt_PE/qsub_cutadapt_PE.sh" else "cutadapt_SE/qsub_cutadapt_SE.sh")
  paste(apply(pairs, 1, function(row) {
    trimmed1 <- file.path(outdir, basename(row[["r1"]]))
    trimmed2 <- if (project$paired_end) file.path(outdir, basename(row[["r2"]])) else trimmed1
    submit_sbatch(project, "Cutadapt", script, c(min_length, adapter1, adapter2, trimmed1, trimmed2, row[["r1"]], row[["r2"]], project$name), "cutadapt")
  }), collapse = "\n")
}

submit_star_jobs <- function(project, trimmed = FALSE) {
  res <- genome_resources(project)
  outdir <- file.path(project$data_dir, "star")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  pairs <- sample_fastq_pairs(project, trimmed)
  if (!NROW(pairs)) return("No samples/read files found in design matrix.")
  script <- file.path(SCRIPTS_DIR, "STAR", if (project$paired_end) "qsub_star_PE.sh" else "qsub_star_SE.sh")
  paste(apply(pairs, 1, function(row) {
    sample_dir <- file.path(outdir, row[["sample"]])
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    out_prefix <- file.path(sample_dir, row[["sample"]])
    submit_sbatch(project, "STAR", script, c(out_prefix, res$star_index, row[["r1"]], row[["r2"]], project$name), "star")
  }), collapse = "\n")
}

submit_kallisto_jobs <- function(project, trimmed = FALSE) {
  res <- genome_resources(project)
  outdir <- file.path(project$data_dir, "kallisto")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  pairs <- sample_fastq_pairs(project, trimmed)
  if (!NROW(pairs)) return("No samples/read files found in design matrix.")
  script <- file.path(SCRIPTS_DIR, "Kallisto", if (project$paired_end) "qsub_kallisto_PE.sh" else "qsub_kallisto_SE.sh")
  paste(apply(pairs, 1, function(row) {
    sample_dir <- file.path(outdir, row[["sample"]])
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    submit_sbatch(project, "Kallisto", script, c(sample_dir, res$kallisto_index, row[["r1"]], row[["r2"]], project$name), "kallisto")
  }), collapse = "\n")
}

submit_featurecounts_jobs <- function(project, feature = "gene_id") {
  res <- genome_resources(project)
  design <- safe_read_table(project$design_matrix_path)
  if (!NROW(design) || !"sample" %in% names(design)) return("No samples found in design matrix.")
  outdir <- file.path(project$data_dir, "featurecounts")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  script <- file.path(SCRIPTS_DIR, "featureCounts", if (project$paired_end) "qsub_featurecounts_PE.sh" else "qsub_featurecounts_SE.sh")
  paste(vapply(as.character(design$sample), function(sample) {
    sample_dir <- file.path(outdir, sample)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    bam <- file.path(project$data_dir, "star", sample, paste0(sample, "Aligned.sortedByCoord.out.bam"))
    count_prefix <- file.path(sample_dir, sample)
    submit_sbatch(project, "featureCounts", script, c(bam, res$gtf, feature, count_prefix, res$strand_bed, project$name), "featurecounts")
  }, character(1)), collapse = "\n")
}

submit_deseq2_job <- function(project, reference, comparison, redundant = "NoRedundant") {
  outdir <- file.path(project$data_dir, "deseq2")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  script <- file.path(SCRIPTS_DIR, "DESeq2", "qsub_deseq2.sh")
  rscript <- file.path(SCRIPTS_DIR, "DESeq2", "DESeq2.R")
  count_matrix <- file.path(project$data_dir, "counts", "count_matrix.txt")
  submit_sbatch(project, "DESeq2", script, c(rscript, count_matrix, project$design_matrix_path, outdir, reference, comparison, redundant, project$name), "deseq2")
}

list_result_files <- function(project, pattern = "\\.(txt|csv|tsv|html|png|pdf)$") {
  if (!dir.exists(project$data_dir)) return(character(0))
  list.files(project$data_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
}

image_or_file_ui <- function(path, height = "900px") {
  if (!file.exists(path)) return(tags$div(class = "empty-box", "File not found."))
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("png", "jpg", "jpeg", "webp") && BASE64_AVAILABLE) {
    mime <- if (ext == "png") "image/png" else "image/jpeg"
    tags$img(src = paste0("data:", mime, ";base64,", base64enc::base64encode(path)), style = "max-width:100%; border:1px solid #d8dde8; border-radius:8px;")
  } else if (ext == "html") {
    html <- paste(readLines(path, warn = FALSE), collapse = "\n")
    tags$iframe(srcdoc = htmltools::HTML(html), style = paste0("width:100%; height:", height, "; border:1px solid #d8dde8; border-radius:8px;"))
  } else if (ext == "pdf" && BASE64_AVAILABLE) {
    tags$iframe(src = paste0("data:application/pdf;base64,", base64enc::base64encode(path)), style = paste0("width:100%; height:", height, "; border:1px solid #d8dde8; border-radius:8px;"))
  } else {
    tags$div(class = "empty-box", tags$p(basename(path)), tags$p(path))
  }
}

app_css <- "
body { background:#f5f7fb; color:#17202f; }
.navbar, .navbar-default { background:#0f1724 !important; border:0; }
.navbar-default .navbar-nav > li > a, .navbar-default .navbar-brand { color:#f8fafc !important; }
.well, .panel, .tab-content { border-radius:8px; border-color:#d8dde8; }
.csl-header { background:white; border:1px solid #d8dde8; border-radius:8px; padding:16px 18px; margin-bottom:14px; }
.csl-header h2 { margin:0 0 6px 0; font-weight:700; }
.muted { color:#657084; }
.empty-box { background:white; border:1px solid #d8dde8; border-radius:8px; padding:18px; color:#657084; }
.btn-primary { background:#1f5eff; border-color:#1f5eff; }
"

ui <- fluidPage(
  tags$head(tags$style(HTML(app_css))),
  div(class = "csl-header",
      h2("CodeSpringWeb"),
      div(class = "muted", "Shiny control center for CodeSpringLab projects: configure, run, track, and visualize results from one port.")
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("analysis", "Analysis", choices = c("RNA-seq", "ATAC-seq", "ChIP-seq", "All analyses"), selected = "RNA-seq"),
      uiOutput("project_ui"),
      tags$hr(),
      div(class = "muted", sprintf("CodeSpringLab root: %s", CSL_ROOT)),
      tags$hr(),
      h4("Selected Project"),
      verbatimTextOutput("project_paths")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Setup", br(), h3("Project Setup"), tableOutput("setup_table"), uiOutput("source_config_ui")),
        tabPanel("Design Matrix", br(), h3("Design Matrix Builder"),
                 fluidRow(
                   column(8, textInput("metadata_cols", "Metadata columns", value = "treatment", placeholder = "treatment, batch, replicate")),
                   column(4, br(), actionButton("scan_fastqs", "Scan FASTQ folder", class = "btn-primary"))
                 ),
                 table_output("design_editor"),
                 br(),
                 actionButton("save_design", "Save design_matrix.txt", class = "btn-primary"),
                 verbatimTextOutput("design_save_status")),
        tabPanel("Progress", br(), h3("Pipeline Progress"), table_output("status_table"), br(), h4("Sample Progress"), table_output("sample_progress_table")),
        tabPanel("Run Pipeline", br(), h3("Run Pipeline"),
                 tags$p(class = "muted", "Buttons submit SLURM sbatch jobs. After submission, jobs keep running even if this app or browser is closed."),
                 fluidRow(
                   column(3, checkboxInput("use_trimmed_reads", "Use trimmed reads", value = FALSE)),
                   column(3, selectInput("feature_attr", "featureCounts attribute", choices = c("gene_id", "gene_name"), selected = "gene_id")),
                   column(3, textInput("deseq_reference", "DESeq2 reference", value = "control")),
                   column(3, textInput("deseq_comparison", "DESeq2 comparison", value = "treated"))
                 ),
                 fluidRow(
                   column(4, textInput("adapter1", "R1 adapter", value = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA")),
                   column(4, textInput("adapter2", "R2 adapter", value = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")),
                   column(4, textInput("min_length", "Minimum read length", value = "20"))
                 ),
                 fluidRow(
                   column(2, actionButton("run_fastqc", "FastQC", class = "btn-primary")),
                   column(2, actionButton("run_cutadapt", "cutadapt", class = "btn-primary")),
                   column(2, actionButton("run_star", "STAR", class = "btn-primary")),
                   column(2, actionButton("run_kallisto", "Kallisto", class = "btn-primary")),
                   column(2, actionButton("run_featurecounts", "featureCounts", class = "btn-primary")),
                   column(2, actionButton("run_deseq2", "DESeq2", class = "btn-primary"))
                 ),
                 br(),
                 verbatimTextOutput("run_output")),
        tabPanel("Results Explorer", br(),
                 tabsetPanel(
                   tabPanel("Overview", br(), table_output("results_overview"), br(), h4("Design Matrix"), table_output("design_table")),
                   tabPanel("QC", br(), uiOutput("fastqc_select_ui"), uiOutput("fastqc_view")),
                   tabPanel("Alignment QC", br(), h4("STAR Summary"), table_output("star_summary"), br(), h4("featureCounts Summary"), table_output("featurecounts_summary")),
                   tabPanel("Counts", br(), tabsetPanel(
                     tabPanel("Raw Counts", br(), table_output("count_matrix")),
                     tabPanel("RSEM", br(), uiOutput("rsem_file_ui"), table_output("rsem_table")),
                     tabPanel("Kallisto", br(), uiOutput("kallisto_file_ui"), table_output("kallisto_table")),
                     tabPanel("DESeq2 Normalized", br(), uiOutput("norm_file_ui"), table_output("norm_table"))
                   )),
                   tabPanel("DESeq2", br(), uiOutput("deseq_file_ui"), uiOutput("deseq_file_view")),
                   tabPanel("GSEA", br(), uiOutput("gsea_file_ui"), uiOutput("gsea_file_view")),
                   tabPanel("Files", br(), uiOutput("all_file_ui"), uiOutput("all_file_view"))
                 )),
        tabPanel("Logs", br(), h3("Submitted Jobs"), table_output("jobs_table"))
      )
    )
  )
)

server <- function(input, output, session) {
  projects <- reactiveVal(discover_projects())
  design_state <- reactiveVal(data.frame())
  run_message <- reactiveVal("")

  filtered_projects <- reactive({
    p <- projects()
    analysis <- input$analysis
    if (!length(analysis) || is.null(analysis) || !nzchar(analysis) || identical(analysis, "All analyses")) return(p)
    p[vapply(p, function(x) identical(x$analysis, analysis), logical(1))]
  })

  output$project_ui <- renderUI({
    p <- filtered_projects()
    if (!length(p)) {
      return(div(class = "muted", "No project configs were found. Check CSL_CODESPRINGLAB_ROOT or create a new project config."))
    }
    labels <- vapply(p, function(x) paste0(x$label, " (", x$analysis, if (nzchar(x$source_config)) " · CSL config" else "", ")"), character(1))
    selectInput("project_id", "Project config", choices = labels, selected = labels[[1]])
  })

  current_project <- reactive({
    p <- filtered_projects()
    req(length(p) > 0)
    labels <- vapply(p, function(x) paste0(x$label, " (", x$analysis, if (nzchar(x$source_config)) " · CSL config" else "", ")"), character(1))
    selected <- input$project_id
    if (!length(selected) || is.null(selected) || !nzchar(selected)) {
      idx <- 1
    } else {
      idx <- match(selected, labels)
      if (!length(idx) || is.na(idx)) idx <- 1
    }
    p[[idx]]
  })

  output$project_paths <- renderText({
    p <- current_project()
    paste(c(
      paste("Project:", p$label),
      paste("Analysis:", p$analysis),
      paste("Genome:", p$genome),
      paste("Data:", p$data_dir),
      paste("Design:", p$design_matrix_path),
      paste("FASTQ:", p$fastq_dir)
    ), collapse = "\n")
  })

  output$setup_table <- renderTable({
    p <- current_project()
    data.frame(
      field = c("Project", "Analysis", "Genome", "Paired-end", "Results root", "Data folder", "FASTQ folder", "Design matrix"),
      value = c(p$label, p$analysis, p$genome, as.character(p$paired_end), p$results_root, p$data_dir, p$fastq_dir, p$design_matrix_path),
      stringsAsFactors = FALSE
    )
  })

  output$source_config_ui <- renderUI({
    p <- current_project()
    if (!nzchar(p$source_config)) return(NULL)
    tagList(h4("Imported CodeSpringLab Config"), tags$pre(p$source_config))
  })

  observeEvent(input$scan_fastqs, {
    p <- current_project()
    cols <- clean_name(unlist(strsplit(input$metadata_cols, ",")))
    cols <- cols[nzchar(cols) & !cols %in% c("sample", "filename", "include", "status")]
    if (!length(cols)) cols <- "treatment"
    design_state(scan_fastqs(p$fastq_dir, p$paired_end, cols))
  })

  observeEvent(current_project(), {
    p <- current_project()
    if (file.exists(p$design_matrix_path)) {
      df <- safe_read_table(p$design_matrix_path)
      if (NROW(df)) {
        df$include <- TRUE
        df$status <- "saved"
        df <- df[, c("include", setdiff(names(df), c("include", "status")), "status"), drop = FALSE]
        design_state(df)
      }
    }
  }, ignoreInit = FALSE)

  output$design_editor <- render_csl_table({
    df <- design_state()
    if (!NROW(df)) df <- data.frame(include = logical(), sample = character(), treatment = character(), filename = character(), status = character())
    df
  }, page_length = 25, editable = TRUE)

  observeEvent(input$design_editor_cell_edit, {
    info <- input$design_editor_cell_edit
    df <- design_state()
    if (NROW(df)) {
      df[info$row, info$col + 1] <- info$value
      design_state(df)
    }
  })

  output$design_save_status <- renderText("")
  observeEvent(input$save_design, {
    p <- current_project()
    df <- design_state()
    metadata <- setdiff(names(df), c("include", "sample", "filename", "status"))
    msg <- tryCatch(write_design_matrix(p, df, metadata), error = function(e) paste("ERROR:", conditionMessage(e)))
    output$design_save_status <- renderText(msg)
  })

  output$status_table <- render_csl_table(project_status(current_project()), page_length = 20)

  output$sample_progress_table <- render_csl_table(sample_progress(current_project()), page_length = 25)

  observeEvent(input$run_fastqc, {
    run_message(submit_fastqc_jobs(current_project(), isTRUE(input$use_trimmed_reads)))
  })
  observeEvent(input$run_cutadapt, {
    run_message(submit_cutadapt_jobs(current_project(), input$adapter1, input$adapter2, input$min_length))
  })
  observeEvent(input$run_star, {
    run_message(submit_star_jobs(current_project(), isTRUE(input$use_trimmed_reads)))
  })
  observeEvent(input$run_kallisto, {
    run_message(submit_kallisto_jobs(current_project(), isTRUE(input$use_trimmed_reads)))
  })
  observeEvent(input$run_featurecounts, {
    run_message(submit_featurecounts_jobs(current_project(), input$feature_attr))
  })
  observeEvent(input$run_deseq2, {
    if (identical(input$deseq_reference, input$deseq_comparison)) {
      run_message("Reference and comparison must be different.")
    } else {
      run_message(submit_deseq2_job(current_project(), input$deseq_reference, input$deseq_comparison, "NoRedundant"))
    }
  })
  output$run_output <- renderText(run_message())

  output$results_overview <- render_csl_table(project_status(current_project()), page_length = 20)
  output$design_table <- render_csl_table(safe_read_table(current_project()$design_matrix_path), page_length = 25)
  output$fastqc_select_ui <- renderUI({
    p <- current_project()
    files <- c(list.files(file.path(p$data_dir, "fastqc"), pattern = "\\.html$", full.names = TRUE),
               list.files(file.path(p$data_dir, "fastqc_cutadapt"), pattern = "\\.html$", full.names = TRUE))
    selectInput("fastqc_file", "FastQC report", choices = files, selected = files[1] %||% character(0))
  })
  output$fastqc_view <- renderUI({ req(input$fastqc_file); image_or_file_ui(input$fastqc_file, "1050px") })
  output$star_summary <- render_csl_table(safe_read_table(file.path(current_project()$data_dir, "star_summary", "summary_matrix.txt")), page_length = 25)
  output$featurecounts_summary <- render_csl_table(safe_read_table(file.path(current_project()$data_dir, "counts", "featurecounts_summary.txt")), page_length = 25)
  output$count_matrix <- render_csl_table(safe_read_table(file.path(current_project()$data_dir, "counts", "count_matrix.txt"), 5000), page_length = 25)

  file_select <- function(id, label, dir, pattern) {
    files <- if (dir.exists(dir)) list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE) else character(0)
    selectInput(id, label, choices = files, selected = files[1] %||% character(0))
  }
  output$rsem_file_ui <- renderUI({ file_select("rsem_file", "RSEM table", file.path(current_project()$data_dir, "rsem"), "\\.(txt|csv|results)$") })
  output$rsem_table <- render_csl_table({ req(input$rsem_file); safe_read_table(input$rsem_file, 5000) }, page_length = 25)
  output$kallisto_file_ui <- renderUI({ file_select("kallisto_file", "Kallisto table", file.path(current_project()$data_dir, "kallisto"), "\\.(tsv|txt|csv)$") })
  output$kallisto_table <- render_csl_table({ req(input$kallisto_file); safe_read_table(input$kallisto_file, 5000) }, page_length = 25)
  output$norm_file_ui <- renderUI({ file_select("norm_file", "DESeq2 normalized counts", file.path(current_project()$data_dir, "deseq2"), "normalized.*\\.(txt|csv)$") })
  output$norm_table <- render_csl_table({ req(input$norm_file); safe_read_table(input$norm_file, 5000) }, page_length = 25)
  output$deseq_file_ui <- renderUI({ file_select("deseq_file", "DESeq2 file", file.path(current_project()$data_dir, "deseq2"), "\\.(txt|csv|png|pdf)$") })
  output$deseq_file_view <- renderUI({
    req(input$deseq_file)
    if (tolower(tools::file_ext(input$deseq_file)) %in% c("txt", "csv", "tsv")) {
      table_output("deseq_selected_table")
    } else image_or_file_ui(input$deseq_file)
  })
  output$deseq_selected_table <- render_csl_table({ req(input$deseq_file); safe_read_table(input$deseq_file, 5000) }, page_length = 25)
  output$gsea_file_ui <- renderUI({ file_select("gsea_file", "GSEA file", file.path(current_project()$data_dir, "gseapy"), "\\.(txt|csv|png|pdf)$") })
  output$gsea_file_view <- renderUI({
    req(input$gsea_file)
    if (tolower(tools::file_ext(input$gsea_file)) %in% c("txt", "csv", "tsv")) {
      table_output("gsea_selected_table")
    } else image_or_file_ui(input$gsea_file, "950px")
  })
  output$gsea_selected_table <- render_csl_table({ req(input$gsea_file); safe_read_table(input$gsea_file, 5000) }, page_length = 25)
  output$all_file_ui <- renderUI({ file_select("all_file", "Result file", current_project()$data_dir, "\\.(txt|csv|tsv|html|png|pdf)$") })
  output$all_file_view <- renderUI({ req(input$all_file); image_or_file_ui(input$all_file) })
  output$jobs_table <- render_csl_table({
    if (!file.exists(JOBS_PATH)) return(data.frame())
    utils::read.delim(JOBS_PATH, check.names = FALSE)
  }, page_length = 25)
}

shinyApp(ui, server)
