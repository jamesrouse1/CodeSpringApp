# Lightweight BAM inspection helpers for the Sarek manifest workflow.
# These functions read only the BAM header and index metadata. They do not scan
# alignments or call variants.

SAREK_BAM_INSPECTOR_VERSION <- "0.1.0"

sarek_bam_command <- function(samtools, args, timeout_seconds = 30L) {
  timeout_seconds <- suppressWarnings(as.integer(timeout_seconds))
  if (!is.finite(timeout_seconds) || timeout_seconds < 1L) timeout_seconds <- 30L

  output <- tryCatch(
    system2(
      samtools,
      args = args,
      stdout = TRUE,
      stderr = TRUE,
      timeout = timeout_seconds
    ),
    error = function(error) structure(conditionMessage(error), status = 127L)
  )
  status <- attr(output, "status")
  if (is.null(status) || !length(status) || is.na(status[[1]])) status <- 0L
  list(status = as.integer(status[[1]]), output = as.character(output))
}

sarek_resolve_samtools <- function(samtools = "samtools") {
  samtools <- sarek_text(samtools, "samtools")
  resolved <- if (grepl("/", samtools, fixed = TRUE)) samtools else Sys.which(samtools)[[1]]
  if (!nzchar(resolved) || !file.exists(resolved) || dir.exists(resolved) || file.access(resolved, mode = 1) != 0) {
    stop(
      "samtools is unavailable. Configure an executable with CSL_SAREK_SAMTOOLS; received: ",
      samtools,
      call. = FALSE
    )
  }
  normalizePath(resolved, winslash = "/", mustWork = TRUE)
}

sarek_bam_header_tag <- function(lines, record, tag) {
  records <- lines[startsWith(lines, paste0(record, "\t"))]
  if (!length(records)) return(character(0))
  values <- unlist(strsplit(records, "\t", fixed = TRUE), use.names = FALSE)
  prefix <- paste0(tag, ":")
  unique(sub(prefix, "", values[startsWith(values, prefix)], fixed = TRUE))
}

sarek_bam_reference_evidence <- function(header_lines) {
  sq <- header_lines[startsWith(header_lines, "@SQ\t")]
  if (!length(sq)) {
    return(list(
      grch38_compatible = NA,
      contig_style = "unknown",
      summary = "No sequence dictionary was present in the BAM header."
    ))
  }

  parsed <- lapply(sq, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    sn <- sub("^SN:", "", fields[startsWith(fields, "SN:")][1])
    ln <- suppressWarnings(as.numeric(sub("^LN:", "", fields[startsWith(fields, "LN:")][1])))
    list(name = if (length(sn)) sn else NA_character_, length = if (length(ln)) ln else NA_real_)
  })
  names <- vapply(parsed, function(item) item$name, character(1))
  lengths <- vapply(parsed, function(item) item$length, numeric(1))

  chr1_index <- match("chr1", names)
  plain1_index <- match("1", names)
  if (!is.na(chr1_index)) {
    compatible <- isTRUE(lengths[[chr1_index]] == 248956422)
    return(list(
      grch38_compatible = compatible,
      contig_style = "chr-prefixed",
      summary = if (compatible) {
        "chr1 length matches GRCh38 and the contig naming matches GATK.GRCh38."
      } else {
        paste0("chr1 length is ", lengths[[chr1_index]], ", not the GRCh38 length 248956422.")
      }
    ))
  }
  if (!is.na(plain1_index)) {
    compatible_length <- isTRUE(lengths[[plain1_index]] == 248956422)
    return(list(
      grch38_compatible = if (compatible_length) NA else FALSE,
      contig_style = "non-chr-prefixed",
      summary = if (compatible_length) {
        "Chromosome 1 length matches GRCh38, but non-chr contig names may not match GATK.GRCh38 resources."
      } else {
        paste0("Chromosome 1 length is ", lengths[[plain1_index]], ", not the GRCh38 length 248956422.")
      }
    ))
  }

  list(
    grch38_compatible = NA,
    contig_style = "other",
    summary = "The BAM header has no chr1 or 1 contig, so GRCh38 compatibility could not be inferred."
  )
}

sarek_infer_bam_processing_state <- function(header_lines) {
  hd_sort <- sarek_bam_header_tag(header_lines, "@HD", "SO")
  sort_order <- if (length(hd_sort)) hd_sort[[1]] else "unknown"
  pg_lines <- header_lines[startsWith(header_lines, "@PG\t")]
  pg_text <- tolower(paste(pg_lines, collapse = " "))

  applied_bqsr <- grepl(
    "applybqsr|apply_bqsr|printreads|recalibrated[.]bam|recalibrated[.]cram",
    pg_text,
    perl = TRUE
  )
  duplicate_marked <- grepl(
    "markduplicates|mark_duplicates|samtools[[:space:]]+markdup|picard[^[:space:]]*[[:space:]]+markduplicates",
    pg_text,
    perl = TRUE
  )
  has_dictionary <- any(startsWith(header_lines, "@SQ\t"))

  if (applied_bqsr) {
    return(list(
      recommendation = "recalibrated",
      confidence = "high",
      sort_order = sort_order,
      evidence = "The BAM program records show that base-quality recalibration was applied."
    ))
  }
  if (duplicate_marked) {
    return(list(
      recommendation = "duplicate_marked",
      confidence = "medium",
      sort_order = sort_order,
      evidence = "The BAM program records show duplicate marking but do not prove that BQSR was applied."
    ))
  }
  if (identical(sort_order, "coordinate") && has_dictionary) {
    return(list(
      recommendation = "aligned",
      confidence = "medium",
      sort_order = sort_order,
      evidence = "The BAM is coordinate-sorted and has a sequence dictionary, but later processing was not proven."
    ))
  }

  list(
    recommendation = "unknown",
    confidence = "none",
    sort_order = sort_order,
    evidence = "The BAM header does not contain enough provenance to recommend a safe processing state."
  )
}

sarek_bam_inspection_result <- function(
  path,
  index,
  status,
  summary,
  quickcheck_ok = FALSE,
  index_ok = FALSE,
  sample_ids = character(0),
  sort_order = "unknown",
  grch38_compatible = NA,
  reference_summary = "Not inspected.",
  recommendation = "unknown",
  confidence = "none",
  evidence = character(0)
) {
  list(
    inspector_version = SAREK_BAM_INSPECTOR_VERSION,
    inspected_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    path = sarek_text(path),
    index = sarek_text(index),
    status = status,
    summary = summary,
    quickcheck_ok = isTRUE(quickcheck_ok),
    index_ok = isTRUE(index_ok),
    sample_ids = unique(sample_ids[nzchar(sample_ids)]),
    sort_order = sarek_text(sort_order, "unknown"),
    grch38_compatible = grch38_compatible,
    reference_summary = sarek_text(reference_summary, "Not inspected."),
    recommendation = sarek_text(recommendation, "unknown"),
    confidence = sarek_text(confidence, "none"),
    evidence = unique(evidence[nzchar(evidence)])
  )
}

sarek_inspect_bam <- function(path, index = "", samtools = "samtools", timeout_seconds = 30L) {
  path <- sarek_text(path)
  index <- sarek_text(index)
  if (!nzchar(path) || !file.exists(path) || dir.exists(path) || file.access(path, mode = 4) != 0) {
    return(sarek_bam_inspection_result(
      path, index, "failed", "The BAM file is missing or unreadable.",
      evidence = "BAM inspection stopped before running samtools."
    ))
  }
  if (!identical(sarek_detect_input_format(path), "bam")) {
    return(sarek_bam_inspection_result(
      path, index, "failed", "The selected file is not a BAM file.",
      evidence = "Only BAM inspection is implemented in this lightweight inspector."
    ))
  }
  if (!nzchar(index) || !file.exists(index) || dir.exists(index) || file.access(index, mode = 4) != 0) {
    return(sarek_bam_inspection_result(
      path, index, "failed", "A readable BAM index was not detected.",
      evidence = "Expected a readable .bam.bai or sibling .bai file."
    ))
  }

  resolved_samtools <- tryCatch(sarek_resolve_samtools(samtools), error = function(error) error)
  if (inherits(resolved_samtools, "error")) {
    return(sarek_bam_inspection_result(
      path, index, "unavailable", conditionMessage(resolved_samtools),
      evidence = "No BAM content was inspected."
    ))
  }

  quoted_path <- shQuote(path)
  quickcheck <- sarek_bam_command(
    resolved_samtools,
    c("quickcheck", "-v", quoted_path),
    timeout_seconds = timeout_seconds
  )
  if (quickcheck$status != 0L) {
    details <- trimws(paste(quickcheck$output, collapse = " "))
    return(sarek_bam_inspection_result(
      path, index, "failed",
      paste0("samtools quickcheck failed", if (nzchar(details)) paste0(": ", details) else "."),
      evidence = "The BAM failed the lightweight integrity check."
    ))
  }

  header <- sarek_bam_command(
    resolved_samtools,
    c("view", "-H", quoted_path),
    timeout_seconds = timeout_seconds
  )
  if (header$status != 0L || !length(header$output)) {
    details <- trimws(paste(header$output, collapse = " "))
    return(sarek_bam_inspection_result(
      path, index, "failed",
      paste0("The BAM header could not be read", if (nzchar(details)) paste0(": ", details) else "."),
      quickcheck_ok = TRUE,
      evidence = "samtools quickcheck passed, but header inspection failed."
    ))
  }

  idxstats <- sarek_bam_command(
    resolved_samtools,
    c("idxstats", quoted_path),
    timeout_seconds = timeout_seconds
  )
  if (idxstats$status != 0L) {
    details <- trimws(paste(idxstats$output, collapse = " "))
    return(sarek_bam_inspection_result(
      path, index, "failed",
      paste0("The BAM index could not be used", if (nzchar(details)) paste0(": ", details) else "."),
      quickcheck_ok = TRUE,
      evidence = "The index file exists, but samtools idxstats could not use it."
    ))
  }

  processing <- sarek_infer_bam_processing_state(header$output)
  reference <- sarek_bam_reference_evidence(header$output)
  sample_ids <- sarek_bam_header_tag(header$output, "@RG", "SM")
  evidence <- c(
    "samtools quickcheck passed.",
    "samtools idxstats confirmed that the BAM index is usable.",
    processing$evidence,
    reference$summary
  )
  summary <- paste0(
    "Inspection passed; suggested processing state: ", processing$recommendation,
    " (", processing$confidence, " confidence)."
  )
  sarek_bam_inspection_result(
    path = path,
    index = index,
    status = "passed",
    summary = summary,
    quickcheck_ok = TRUE,
    index_ok = TRUE,
    sample_ids = sample_ids,
    sort_order = processing$sort_order,
    grch38_compatible = reference$grch38_compatible,
    reference_summary = reference$summary,
    recommendation = processing$recommendation,
    confidence = processing$confidence,
    evidence = evidence
  )
}

sarek_apply_bam_inspection <- function(table, path, inspection) {
  if (!is.data.frame(table) || !NROW(table)) stop("A non-empty discovery table is required.", call. = FALSE)
  path <- sarek_text(path)
  rows <- which(as.character(table$path) == path)
  if (length(rows) != 1L) stop("The inspected BAM must match exactly one discovery row.", call. = FALSE)
  if (!identical(as.character(table$input_format[[rows]]), "bam")) {
    stop("Inspection results can only be attached to a BAM row.", call. = FALSE)
  }

  fields <- list(
    inspection_status = sarek_text(inspection$status, "failed"),
    processing_recommendation = sarek_text(inspection$recommendation, "unknown"),
    processing_confidence = sarek_text(inspection$confidence, "none"),
    inspection_summary = sarek_text(inspection$summary),
    header_sample_ids = paste(inspection$sample_ids, collapse = ", "),
    sort_order = sarek_text(inspection$sort_order, "unknown"),
    reference_compatibility = sarek_text(inspection$reference_summary),
    inspection_evidence = paste(inspection$evidence, collapse = " ")
  )
  for (field in names(fields)) {
    if (!field %in% names(table)) table[[field]] <- ""
    table[[field]][rows] <- fields[[field]]
  }
  table
}

sarek_auto_inspect_bams <- function(
  table,
  inspector = sarek_inspect_bam,
  samtools = "samtools",
  max_bams = 20L
) {
  if (!is.data.frame(table) || !NROW(table)) return(table)
  bam_rows <- which(as.character(table$input_format) == "bam" & suppressWarnings(as.logical(table$include)))
  if (!length(bam_rows)) return(table)
  max_bams <- suppressWarnings(as.integer(max_bams))
  if (!is.finite(max_bams) || max_bams < 0L) max_bams <- 0L
  inspect_rows <- utils::head(bam_rows, max_bams)

  for (row in inspect_rows) {
    inspection <- tryCatch(
      inspector(
        path = as.character(table$path[[row]]),
        index = as.character(table$index[[row]]),
        samtools = samtools
      ),
      error = function(error) sarek_bam_inspection_result(
        table$path[[row]], table$index[[row]], "failed", conditionMessage(error),
        evidence = "The inspector returned an unexpected error."
      )
    )
    table <- sarek_apply_bam_inspection(table, table$path[[row]], inspection)
  }
  if (length(bam_rows) > length(inspect_rows)) {
    deferred <- setdiff(bam_rows, inspect_rows)
    table$inspection_status[deferred] <- "deferred"
    table$inspection_summary[deferred] <- paste0(
      "Automatic inspection was deferred because discovery exceeded the ", max_bams, "-BAM safety limit."
    )
  }
  table
}
