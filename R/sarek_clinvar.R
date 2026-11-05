# Optional, offline ClinVar annotation of one completed Sarek VCF. ClinVar
# references are supplied by the site/user; this module never downloads data.
sarek_clinvar_header <- function(path, max_lines = 10000L) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  header <- character()
  contig <- ""
  found_columns <- FALSE
  for (i in seq_len(max_lines)) {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    if (startsWith(line, "##")) { header <- c(header, line); next }
    if (startsWith(line, "#CHROM")) {
      found_columns <- TRUE
      record <- readLines(con, n = 1L, warn = FALSE)
      if (length(record)) contig <- strsplit(record[[1]], "\t", fixed = TRUE)[[1]][[1]]
      break
    }
    stop("Expected a VCF header in ", basename(path), ".", call. = FALSE)
  }
  if (!length(header) || !found_columns) stop("The VCF header is incomplete in ", basename(path), ".", call. = FALSE)
  if (!nzchar(contig)) {
    ids <- sub('^##contig=<ID=([^,>]+).*$', '\\1', grep('^##contig=<ID=', header, value = TRUE))
    if (length(ids)) contig <- ids[[1]]
  }
  list(lines = header, contig = contig)
}

sarek_clinvar_assembly <- function(manifest, declared, reference_header) {
  run_assembly <- sarek_text(manifest$reference$assembly)
  if (!run_assembly %in% c("GRCh37", "GRCh38")) {
    stop("ClinVar VCF annotation requires a Sarek run with a known GRCh37 or GRCh38 assembly.", call. = FALSE)
  }
  if (!identical(sarek_text(declared), run_assembly)) {
    stop("Select the ClinVar VCF assembly that matches this run: ", run_assembly, ".", call. = FALSE)
  }
  build_line <- grep('^##(reference|assembly)=', reference_header, value = TRUE, ignore.case = TRUE)
  known <- unique(c(if (any(grepl('GRCh37|hg19', build_line, ignore.case = TRUE))) "GRCh37",
                    if (any(grepl('GRCh38|hg38', build_line, ignore.case = TRUE))) "GRCh38"))
  if (length(known) && !identical(known, run_assembly)) {
    stop("ClinVar VCF header reports a different genome assembly from this run.", call. = FALSE)
  }
  run_assembly
}

sarek_clinvar_contigs <- function(header) {
  lines <- grep("^##contig=<ID=", header$lines, value = TRUE)
  ids <- sub("^##contig=<ID=([^,>]+).*$", "\\1", lines)
  unique(c(ids, header$contig)[nzchar(c(ids, header$contig))])
}

# Translate the source VCF into ClinVar's chromosome convention for annotation.
# The annotation output is translated back to the original source convention.
sarek_clinvar_chromosome_map <- function(source_header, reference_header) {
  empty <- data.frame(
    source = character(),
    clinvar = character(),
    stringsAsFactors = FALSE
  )

  source_chr <- startsWith(source_header$contig, "chr")
  reference_chr <- startsWith(reference_header$contig, "chr")

  if (!nzchar(source_header$contig) ||
      !nzchar(reference_header$contig) ||
      identical(source_chr, reference_chr)) {
    return(empty)
  }

  source_ids <- sarek_clinvar_contigs(source_header)
  reference_ids <- sarek_clinvar_contigs(reference_header)

  if (!any(startsWith(source_header$lines, "##contig=<ID="))) {
    inferred <- if (reference_chr) {
      sub("^chr", "", reference_ids)
    } else {
      paste0("chr", reference_ids)
    }

    if (length(inferred) == 1L) {
      inferred <- if (source_chr) {
        paste0("chr", c(1:22, "X", "Y", "M"))
      } else {
        c(1:22, "X", "Y", "MT")
      }
    }

    if (source_chr && "MT" %in% reference_ids) {
      inferred[inferred == "chrMT"] <- "chrM"
    }

    if (!source_chr && "chrM" %in% reference_ids) {
      inferred[inferred == "M"] <- "MT"
    }

    source_ids <- unique(c(source_ids, inferred))
  }

  changing <- if (reference_chr) {
    !startsWith(source_ids, "chr")
  } else {
    startsWith(source_ids, "chr")
  }

  old <- source_ids[changing]

  new <- if (reference_chr) {
    paste0("chr", old)
  } else {
    substring(old, 4L)
  }

  if (reference_chr) {
    new[new == "chrMT" & "chrM" %in% reference_ids] <- "chrM"
  }

  if (!reference_chr) {
    new[new == "M" & "MT" %in% reference_ids] <- "MT"
  }

  if (any(!nzchar(new)) ||
      anyDuplicated(c(new, source_ids[!changing]))) {
    stop(
      "Cannot safely reconcile VCF chromosome names: renaming would create duplicate contigs.",
      call. = FALSE
    )
  }

  data.frame(
    source = old,
    clinvar = new,
    stringsAsFactors = FALSE
  )
}


sarek_clinvar_sources <- function(catalog) {
  if (!NROW(catalog)) return(catalog)
  keep <- grepl("\\.vcf(\\.gz)?$", catalog$file, ignore.case = TRUE) &
    !grepl("\\.g\\.vcf(\\.gz)?$|\\.clinvar\\.vcf(\\.gz)?$", catalog$file, ignore.case = TRUE) &
    !grepl("(^|/)clinvar(/|$)", as.character(catalog$folder), ignore.case = TRUE)
  catalog[which(keep), , drop = FALSE]
}

sarek_clinvar_outputs <- function(catalog) {
  if (!NROW(catalog)) return(catalog)
  catalog[which(grepl("\\.clinvar\\.vcf\\.gz$", catalog$file, ignore.case = TRUE) &
    grepl("(^|/)annotation/clinvar(/|$)", paste0("/", catalog$folder), ignore.case = TRUE)), , drop = FALSE]
}

sarek_clinvar_log_file <- function(root, record, suffix) {
  prefix <- sub("\\.tsv$", "_", basename(record))
  candidates <- list.files(
    root,
    pattern = paste0("^", prefix, ".*\\.", suffix, "$"),
    full.names = TRUE
  )

  if (!length(candidates)) return("")

  info <- file.info(candidates)
  candidates[[order(info$mtime, decreasing = TRUE, na.last = TRUE)[[1]]]]
}

sarek_clinvar_log_tail <- function(path, max_lines = 120L) {
  if (!nzchar(path) || !file.exists(path)) {
    return("Log file not available yet.")
  }
  if (file.access(path, 4L) != 0L) {
    return("Log file exists but cannot be read.")
  }

  tryCatch({
    lines <- sarek_tail_text_lines(
      path,
      max_lines = max_lines,
      max_bytes = 262144L
    )
    if (!length(lines)) return("No log text available.")
    paste(lines, collapse = "\n")
  }, error = function(error) {
    paste("Could not read log:", conditionMessage(error))
  })
}

sarek_clinvar_slurm_status <- function(job) {
  status <- tryCatch(
    suppressWarnings(sarek_query_slurm_job(job)),
    error = function(error) list()
  )

  state <- if (!is.null(status$state)) sarek_text(status$state) else ""
  exit_code <- if (!is.null(status$exit_code)) sarek_text(status$exit_code) else ""

  sacct <- Sys.which("sacct")

  if (nzchar(sacct) && (!nzchar(state) || !nzchar(exit_code))) {
    raw <- tryCatch(
      system2(
        sacct,
        c(
          "-j", job,
          "-X",
          "--noheader",
          "--parsable2",
          "--format=State,ExitCode"
        ),
        stdout = TRUE,
        stderr = FALSE
      ),
      error = function(error) character()
    )

    raw <- raw[nzchar(trimws(raw))]

    if (length(raw)) {
      fields <- strsplit(raw[[length(raw)]], "\\|")[[1]]

      if (!nzchar(state) && length(fields) >= 1L) {
        state <- trimws(fields[[1]])
      }

      if (!nzchar(exit_code) && length(fields) >= 2L) {
        exit_code <- trimws(fields[[2]])
      }
    }
  }

  list(
    state = sarek_text(state, "SUBMITTED / STATUS UNAVAILABLE"),
    exit_code = exit_code
  )
}

sarek_clinvar_jobs <- function(run) {
  root <- file.path(sarek_text(run$run_dir), ".codespring", "clinvar")

  records <- if (dir.exists(root)) {
    list.files(root, pattern = "\\.tsv$", full.names = TRUE)
  } else {
    character()
  }

  records <- records[!grepl("/status\\.tsv$", records)]

  empty <- data.frame(
    Source = character(),
    Output = character(),
    Job = character(),
    Status = character(),
    ExitCode = character(),
    ErrorLog = character(),
    OutputLog = character(),
    LogText = character(),
    Path = character(),
    stringsAsFactors = FALSE
  )

  if (!length(records)) return(empty)

  rows <- lapply(records, function(record) {
    info <- sarek_read_key_value_file(record)

    output_path <- sarek_text(info["output"])
    source_path <- sarek_text(info["source"])
    job <- sarek_text(info["job_id"])

    error_log <- sarek_clinvar_log_file(root, record, "err")
    output_log <- sarek_clinvar_log_file(root, record, "out")

    slurm <- sarek_clinvar_slurm_status(job)

    output_complete <- nzchar(output_path) &&
      file.exists(output_path) &&
      file.exists(paste0(output_path, ".csi"))

    # Each attempt retains its own Slurm state, even after a retry succeeds.
    state <- slurm$state
    if (identical(state, "COMPLETED") && !output_complete) {
      state <- "COMPLETED (OUTPUT INCOMPLETE)"
    }

    error_text <- sarek_clinvar_log_tail(error_log)
    output_text <- sarek_clinvar_log_tail(output_log)

    log_text <- paste0(
      "[stderr — latest 120 lines, up to 256 KB]\n",
      error_text,
      "\n\n[stdout — latest 120 lines, up to 256 KB]\n",
      output_text
    )

    data.frame(
      Source = basename(source_path),
      Output = basename(output_path),
      Job = job,
      Status = state,
      ExitCode = slurm$exit_code,
      ErrorLog = if (nzchar(error_log)) basename(error_log) else "not available",
      OutputLog = if (nzchar(output_log)) basename(output_log) else "not available",
      LogText = log_text,
      Path = output_path,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}


sarek_clinvar_submit <- function(run, manifest, source, reference, declared_assembly,
                                 bcftools = Sys.getenv("CSL_SAREK_BCFTOOLS", unset = "bcftools"),
                                 queue = Sys.getenv("CSL_SAREK_QUEUE", unset = "cpuq"),
                                 sbatch = "sbatch", runner = NULL) {
  run_dir <- normalizePath(sarek_text(run$run_dir), winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(sarek_text(run$output_dir), winslash = "/", mustWork = TRUE)
  if (!sarek_path_is_within(output_dir, run_dir) ||
      !file.exists(file.path(run_dir, ".codespring", "submission.tsv"))) {
    stop("The selected results directory is not a recorded Sarek run.", call. = FALSE)
  }
  source <- sarek_explorer_file(source, output_dir)
  if (!nzchar(source) || !grepl("\\.vcf(\\.gz)?$", source, ignore.case = TRUE) ||
      grepl("\\.g\\.vcf(\\.gz)?$|\\.clinvar\\.vcf(\\.gz)?$", source, ignore.case = TRUE)) {
    stop("Choose an original Sarek variant VCF for ClinVar annotation.", call. = FALSE)
  }
  reference <- sarek_text(reference)
  if (!sarek_is_absolute_path(reference) || !grepl("\\.vcf\\.gz$", reference, ignore.case = TRUE) ||
      !file.exists(reference) || file.access(reference, 4L) != 0L ||
      !any(file.exists(paste0(reference, c(".tbi", ".csi"))))) {
    stop("Provide a readable, indexed ClinVar .vcf.gz file on the server (.tbi or .csi required).", call. = FALSE)
  }
  reference <- normalizePath(reference, winslash = "/", mustWork = TRUE)
  source_header <- sarek_clinvar_header(source)
  reference_header <- sarek_clinvar_header(reference)
  assembly <- sarek_clinvar_assembly(manifest, declared_assembly, reference_header$lines)
  needed <- c("CLNSIG", "CLNREVSTAT", "CLNDN")
  missing <- needed[!vapply(needed, function(tag) any(grepl(
    paste0("^##INFO=<ID=", tag, "[,>]"), reference_header$lines)), logical(1))]
  if (length(missing)) stop("ClinVar VCF lacks expected INFO tags: ", paste(missing, collapse = ", "), call. = FALSE)
  existing <- needed[vapply(needed, function(tag) any(grepl(
    paste0("^##INFO=<ID=", tag, "[,>]"), source_header$lines)), logical(1))]
  if (length(existing)) stop("Source VCF already defines ClinVar INFO fields: ",
                             paste(existing, collapse = ", "), ". Choose an unannotated VCF.", call. = FALSE)
  chromosome_map <- sarek_clinvar_chromosome_map(
    source_header,
    reference_header
  )
  tool <- sarek_text(bcftools, "bcftools")
  if (sarek_is_absolute_path(tool) && file.access(tool, 1L) != 0L) {
    stop("The configured CSL_SAREK_BCFTOOLS path is not executable.", call. = FALSE)
  }
  if (!grepl("^[A-Za-z0-9_.-]+$", queue)) stop("Invalid Slurm queue.", call. = FALSE)
  relative <- substring(dirname(source), nchar(output_dir) + 2L)
  destination <- file.path(output_dir, "annotation", "clinvar", relative)
  target <- file.path(destination, paste0(sub("\\.vcf(\\.gz)?$", "", basename(source), ignore.case = TRUE),
                                          ".clinvar.vcf.gz"))
  if (file.exists(target)) stop("This VCF already has a ClinVar output: ", basename(target), call. = FALSE)
  previous <- sarek_clinvar_jobs(run)
  if (NROW(previous) && any(previous$Path == target &
                           !previous$Status %in% c("FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY",
                                                        "NODE_FAIL", "BOOT_FAIL", "DEADLINE",
                                                        "Incomplete", "COMPLETED (OUTPUT INCOMPLETE)"))) {
    stop("A ClinVar job has already been submitted for this VCF. Check its status before retrying.", call. = FALSE)
  }
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  internal <- file.path(run_dir, ".codespring", "clinvar")
  dir.create(internal, recursive = TRUE, showWarnings = FALSE)
  attempt <- paste0(format(Sys.time(), "%Y%m%d%H%M%S", tz = "UTC"), "_", Sys.getpid())
  script <- file.path(internal, paste0("clinvar_", attempt, ".sh"))
  temp <- file.path(destination, paste0(".", basename(target), ".", attempt, ".partial.vcf.gz"))
  lines <- c("#!/bin/bash -l", "set -euo pipefail")
  if (!sarek_is_absolute_path(tool)) lines <- c(lines,
    "if ! type module >/dev/null 2>&1; then",
    "  for init in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh; do",
    "    if [[ -r \"$init\" ]]; then source \"$init\"; break; fi",
    "  done",
    "fi",
    "type module >/dev/null 2>&1 || { echo 'Environment Modules is unavailable in this Slurm job.' >&2; exit 127; }",
    "if ! module load EB5 >/dev/null 2>&1; then module load EB5Modules; fi",
    "module load BCFtools/1.23.1-GCC-15.2.0",
    "command -v bcftools >/dev/null 2>&1 || { echo 'BCFtools module loaded, but bcftools is unavailable.' >&2; exit 127; }")
  lines <- c(lines,
    paste0("source_vcf=", shQuote(source)),
    paste0("clinvar_vcf=", shQuote(reference)),
    paste0("out=", shQuote(target)),
    paste0("tmp=", shQuote(temp)),
    "trap 'rm -f \"$tmp\" \"$tmp.csi\"' EXIT")

  annotate <- paste(
    shQuote(tool),
    "annotate --pair-logic exact -a \"$clinvar_vcf\"",
    "-c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN"
  )

  if (NROW(chromosome_map)) {
    forward <- file.path(
      internal,
      paste0("clinvar_", attempt, ".forward.map")
    )

    reverse <- file.path(
      internal,
      paste0("clinvar_", attempt, ".reverse.map")
    )

    normalized <- file.path(
      destination,
      paste0(".", basename(target), ".", attempt, ".normalized.vcf.gz")
    )

    annotated <- file.path(
      destination,
      paste0(".", basename(target), ".", attempt, ".annotated.vcf.gz")
    )

    sarek_write_text_atomic(
      paste(chromosome_map$source, chromosome_map$clinvar, sep = "\t"),
      forward
    )

    sarek_write_text_atomic(
      paste(chromosome_map$clinvar, chromosome_map$source, sep = "\t"),
      reverse
    )

    lines <- c(
      lines,
      paste0("normalized_vcf=", shQuote(normalized)),
      paste0("annotated_vcf=", shQuote(annotated)),
      "trap 'rm -f \"$tmp\" \"$tmp.csi\" \"$normalized_vcf\" \"$normalized_vcf.csi\" \"$annotated_vcf\" \"$annotated_vcf.csi\"' EXIT",
      paste(
        shQuote(tool),
        "annotate --rename-chrs",
        shQuote(forward),
        "-Oz -o \"$normalized_vcf\" \"$source_vcf\""
      ),
      paste(
        shQuote(tool),
        "index --csi --force \"$normalized_vcf\""
      ),
      paste(
        annotate,
        "-Oz -o \"$annotated_vcf\" \"$normalized_vcf\""
      ),
      paste(
        shQuote(tool),
        "index --csi --force \"$annotated_vcf\""
      ),
      paste(
        shQuote(tool),
        "annotate --rename-chrs",
        shQuote(reverse),
        "-Oz -o \"$tmp\" \"$annotated_vcf\""
      )
    )
  } else {
    lines <- c(
      lines,
      paste(annotate, "-Oz -o \"$tmp\" \"$source_vcf\"")
    )
  }

  lines <- c(
    lines,
    paste(shQuote(tool), "index --csi --force \"$tmp\""),
    "mv -f \"$tmp.csi\" \"$out.csi\"",
    "mv -f \"$tmp\" \"$out\""
  )

  sarek_write_text_atomic(lines, script)
  Sys.chmod(script, "0700")
  arguments <- c("--parsable", paste0("--job-name=clinvar_", substr(sarek_identifier(run$run_id, "sarek"), 1L, 40L)),
                 paste0("--partition=", queue), "--time=08:00:00", "--cpus-per-task=2", "--mem=8G",
                 paste0("--output=", file.path(internal, paste0("clinvar_", attempt, "_%j.out"))),
                 paste0("--error=", file.path(internal, paste0("clinvar_", attempt, "_%j.err"))), script)
  result <- if (is.function(runner)) runner(sbatch, arguments) else
    suppressWarnings(system2(sbatch, shQuote(arguments), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0L) {
    stop("Slurm rejected the ClinVar job: ", paste(result, collapse = " "), call. = FALSE)
  }
  job_id <- sarek_parse_slurm_job_id(result)
  if (!nzchar(job_id)) stop("Slurm did not return a ClinVar job ID: ", paste(result, collapse = " "), call. = FALSE)
  sarek_submission_write_values(file.path(internal, paste0("clinvar_", attempt, ".tsv")), c(
    source = source, output = target, reference = reference, assembly = assembly,
    submitted_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), job_id = job_id))
  list(job_id = job_id, output = target)
}
