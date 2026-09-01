# Pure conversion helpers for turning a confirmed CodeSpring Sarek manifest
# into an nf-core/sarek 3.9.0 samplesheet. This module writes launch inputs;
# it does not invoke Nextflow or submit work to Slurm.

SAREK_NEXTFLOW_INPUT_VERSION <- "0.1.0"
SAREK_NEXTFLOW_PIPELINE <- "nf-core/sarek"
SAREK_NEXTFLOW_PIPELINE_VERSION <- "3.9.0"

sarek_nextflow_require_helpers <- function() {
  required <- c("sarek_text", "sarek_is_absolute_path", "sarek_normalize_sex")
  missing <- required[!vapply(required, exists, logical(1), mode = "function")]
  if (length(missing)) {
    stop(
      "Source R/sarek_manifest.R before R/sarek_nextflow_input.R. Missing helpers: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  constants <- c("SAREK_MANIFEST_SCHEMA_VERSION", "SAREK_SEX_CHROMOSOMES")
  missing_constants <- constants[!vapply(constants, exists, logical(1), inherits = TRUE)]
  if (length(missing_constants)) {
    stop(
      "Source R/sarek_manifest.R before R/sarek_nextflow_input.R. Missing constants: ",
      paste(missing_constants, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

sarek_nextflow_read_manifest <- function(path) {
  sarek_nextflow_require_helpers()
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The jsonlite R package is required to read a Sarek manifest.", call. = FALSE)
  }
  path <- sarek_text(path)
  if (!nzchar(path) || !file.exists(path) || dir.exists(path)) {
    stop("Manifest JSON does not exist: ", path, call. = FALSE)
  }
  jsonlite::read_json(path, simplifyVector = FALSE)
}

sarek_nextflow_manifest_samples <- function(manifest) {
  if (!is.list(manifest)) stop("The confirmed manifest must be an R list.", call. = FALSE)
  if (!identical(sarek_text(manifest$schema_version), SAREK_MANIFEST_SCHEMA_VERSION)) {
    stop(
      "Unsupported manifest schema version: ", sarek_text(manifest$schema_version, "missing"),
      ". Expected ", SAREK_MANIFEST_SCHEMA_VERSION, ".",
      call. = FALSE
    )
  }
  if (!identical(sarek_text(manifest$status), "confirmed")) {
    stop("Only a confirmed manifest can be converted to a Nextflow input.", call. = FALSE)
  }
  if (!length(manifest$patients)) stop("The confirmed manifest contains no patients.", call. = FALSE)

  records <- list()
  for (patient in manifest$patients) {
    patient_id <- sarek_text(patient$patient_id)
    if (!nzchar(patient_id)) stop("A manifest patient is missing patient_id.", call. = FALSE)
    patient_sex <- sarek_normalize_sex(patient$sex, "NA")
    if (!patient_sex %in% SAREK_SEX_CHROMOSOMES) {
      stop(
        "Patient ", patient_id, " has unsupported sex chromosomes: ", patient_sex,
        ". Expected XX, XY, or NA.",
        call. = FALSE
      )
    }
    if (!length(patient$samples)) {
      stop("Patient ", patient_id, " contains no samples.", call. = FALSE)
    }
    for (sample in patient$samples) {
      sample_id <- sarek_text(sample$sample_id)
      if (!nzchar(sample_id)) stop("Patient ", patient_id, " contains a sample without sample_id.", call. = FALSE)
      files <- sample$files
      if (!length(files)) {
        stop("Sample ", patient_id, "::", sample_id, " contains no input files.", call. = FALSE)
      }
      records[[length(records) + 1L]] <- list(
        patient = patient_id,
        sex = patient_sex,
        sample = sample_id,
        role = sarek_text(sample$role),
        input_format = sarek_text(sample$input_format),
        processing_state = sarek_text(sample$processing_state),
        files = files
      )
    }
  }
  records
}

sarek_nextflow_validate_identifier <- function(value, label) {
  if (!grepl("^\\S+$", value, perl = TRUE)) {
    stop(label, " cannot be empty or contain whitespace: ", value, call. = FALSE)
  }
  invisible(TRUE)
}

sarek_nextflow_file_path <- function(file, sample_key) {
  path <- sarek_text(file$path)
  if (!nzchar(path) || !sarek_is_absolute_path(path)) {
    stop("Input path must be absolute for sample ", sample_key, ": ", path, call. = FALSE)
  }
  if (grepl("\\s", path, perl = TRUE)) {
    stop("nf-core/sarek input paths cannot contain whitespace: ", path, call. = FALSE)
  }
  if (!file.exists(path) || dir.exists(path)) {
    stop("Input file no longer exists for sample ", sample_key, ": ", path, call. = FALSE)
  }
  if (file.access(path, mode = 4) != 0) {
    stop("Input file is not readable for sample ", sample_key, ": ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

sarek_nextflow_index_path <- function(file, sample_key, extension) {
  path <- sarek_text(file$index)
  if (!nzchar(path)) {
    stop("A ", extension, " index is required for sample ", sample_key, ".", call. = FALSE)
  }
  if (!sarek_is_absolute_path(path) || grepl("\\s", path, perl = TRUE)) {
    stop("The index path must be absolute and contain no whitespace: ", path, call. = FALSE)
  }
  if (!file.exists(path) || dir.exists(path) || file.access(path, mode = 4) != 0) {
    stop("The index file is missing or unreadable for sample ", sample_key, ": ", path, call. = FALSE)
  }
  expected <- paste0("\\.", extension, "$")
  if (!grepl(expected, path, ignore.case = TRUE, perl = TRUE)) {
    stop("Expected a .", extension, " index for sample ", sample_key, ": ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

sarek_nextflow_step_for_sample <- function(record) {
  format <- record$input_format
  state <- record$processing_state
  sample_key <- paste(record$patient, record$sample, sep = "::")

  if (identical(format, "fastq")) {
    if (!identical(state, "unmapped")) {
      stop("FASTQ sample ", sample_key, " must use processing_state 'unmapped'.", call. = FALSE)
    }
    return("mapping")
  }
  if (identical(format, "ubam")) {
    stop(
      "Sample ", sample_key,
      " uses a .ubam filename. nf-core/sarek 3.9.0 accepts unmapped BAM content in the bam column, ",
      "but its input schema requires a .bam filename. A .bam alias must be created before conversion.",
      call. = FALSE
    )
  }
  if (identical(format, "bam")) {
    return(switch(
      state,
      unmapped = "mapping",
      aligned = "markduplicates",
      duplicate_marked = "prepare_recalibration",
      recalibrated = "variant_calling",
      analysis_ready = "variant_calling",
      stop("BAM sample ", sample_key, " has an unsupported processing state: ", state, call. = FALSE)
    ))
  }
  if (identical(format, "cram")) {
    return(switch(
      state,
      aligned = "markduplicates",
      duplicate_marked = "prepare_recalibration",
      recalibrated = "variant_calling",
      analysis_ready = "variant_calling",
      stop("CRAM sample ", sample_key, " has an unsupported processing state: ", state, call. = FALSE)
    ))
  }
  if (identical(format, "vcf")) {
    if (!identical(state, "variant_calls")) {
      stop("VCF sample ", sample_key, " must use processing_state 'variant_calls'.", call. = FALSE)
    }
    return("annotate")
  }
  if (identical(format, "bcf")) {
    stop(
      "Sample ", sample_key,
      " uses BCF, but nf-core/sarek 3.9.0 annotation input accepts VCF only. Convert it to sorted .vcf.gz first.",
      call. = FALSE
    )
  }
  stop("Sample ", sample_key, " has an unsupported input format: ", format, call. = FALSE)
}

sarek_nextflow_status <- function(role, sample_key) {
  if (identical(role, "tumor")) return(1L)
  if (role %in% c("normal", "germline")) return(0L)
  stop("Sample role must be confirmed before Nextflow conversion: ", sample_key, call. = FALSE)
}

sarek_nextflow_validate_mode <- function(records, analysis_mode, step) {
  roles <- vapply(records, function(record) record$role, character(1))
  patients <- vapply(records, function(record) record$patient, character(1))

  if (identical(analysis_mode, "annotation_only")) {
    if (!identical(step, "annotate")) {
      stop("Annotation-only manifests must contain VCF inputs at the annotation step.", call. = FALSE)
    }
    return(invisible(TRUE))
  }
  if (identical(step, "annotate")) {
    stop("VCF inputs require analysis mode 'annotation_only'.", call. = FALSE)
  }
  if (any(!roles %in% c("germline", "tumor", "normal"))) {
    stop("Every sample role must be confirmed before Nextflow conversion.", call. = FALSE)
  }
  if (identical(analysis_mode, "germline") && any(roles != "germline")) {
    stop("Germline mode can contain only germline samples.", call. = FALSE)
  }
  if (identical(analysis_mode, "tumor_only") && any(roles != "tumor")) {
    stop("Tumor-only mode can contain only tumor samples.", call. = FALSE)
  }
  if (identical(analysis_mode, "matched_tumor_normal")) {
    by_patient <- split(roles, patients)
    invalid <- names(by_patient)[!vapply(by_patient, function(x) {
      any(x == "tumor") && any(x == "normal")
    }, logical(1))]
    if (length(invalid)) {
      stop(
        "Matched tumor-normal conversion requires tumor and normal samples for every patient: ",
        paste(invalid, collapse = ", "),
        call. = FALSE
      )
    }
  }
  if (!analysis_mode %in% c("germline", "tumor_only", "matched_tumor_normal")) {
    stop("Unsupported analysis mode: ", analysis_mode, call. = FALSE)
  }
  invisible(TRUE)
}

sarek_nextflow_base_row <- function(record) {
  sample_key <- paste(record$patient, record$sample, sep = "::")
  data.frame(
    patient = record$patient,
    sex = record$sex,
    status = sarek_nextflow_status(record$role, sample_key),
    sample = record$sample,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

sarek_nextflow_fastq_rows <- function(record) {
  sample_key <- paste(record$patient, record$sample, sep = "::")
  file_data <- lapply(record$files, function(file) {
    path <- sarek_nextflow_file_path(file, sample_key)
    if (!grepl("\\.(fastq|fq)\\.gz$", path, ignore.case = TRUE, perl = TRUE)) {
      stop(
        "nf-core/sarek 3.9.0 requires gzip-compressed FASTQ input (.fastq.gz or .fq.gz): ",
        path,
        call. = FALSE
      )
    }
    list(
      path = path,
      lane = sarek_text(file$lane),
      read = if (is.null(file$read) || !length(file$read)) NA_integer_ else suppressWarnings(as.integer(file$read[[1]]))
    )
  })
  lanes <- vapply(file_data, function(file) file$lane, character(1))
  lanes[!nzchar(lanes)] <- "lane_1"
  reads <- vapply(file_data, function(file) file$read, integer(1))
  paths <- vapply(file_data, function(file) file$path, character(1))
  rows <- list()
  for (lane in unique(lanes)) {
    selected <- which(lanes == lane)
    lane_reads <- reads[selected]
    if (length(selected) != 2L || anyNA(lane_reads) || !identical(sort(lane_reads), c(1L, 2L))) {
      stop(
        "FASTQ lane ", sample_key, "::", lane,
        " must contain exactly one R1 and one R2.",
        call. = FALSE
      )
    }
    base <- sarek_nextflow_base_row(record)
    rows[[length(rows) + 1L]] <- data.frame(
      base,
      lane = lane,
      fastq_1 = paths[selected[lane_reads == 1L]],
      fastq_2 = paths[selected[lane_reads == 2L]],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  do.call(rbind, rows)
}

sarek_nextflow_unmapped_bam_rows <- function(record) {
  sample_key <- paste(record$patient, record$sample, sep = "::")
  paths <- vapply(record$files, sarek_nextflow_file_path, character(1), sample_key = sample_key)
  order_index <- order(paths)
  paths <- paths[order_index]
  files <- record$files[order_index]
  lanes <- vapply(files, function(file) sarek_text(file$lane), character(1))
  for (i in seq_along(lanes)) if (!nzchar(lanes[[i]])) lanes[[i]] <- paste0("lane_", i)
  if (anyDuplicated(lanes)) {
    stop("Unmapped BAM lane IDs must be unique within sample ", sample_key, ".", call. = FALSE)
  }
  base <- sarek_nextflow_base_row(record)
  do.call(rbind, lapply(seq_along(paths), function(i) {
    data.frame(base, lane = lanes[[i]], bam = paths[[i]], stringsAsFactors = FALSE, check.names = FALSE)
  }))
}

sarek_nextflow_alignment_row <- function(record) {
  sample_key <- paste(record$patient, record$sample, sep = "::")
  if (length(record$files) != 1L) {
    stop(
      "Processed BAM/CRAM samples must contain exactly one alignment file: ", sample_key,
      ". Merge lanes before using this processing state.",
      call. = FALSE
    )
  }
  file <- record$files[[1]]
  path <- sarek_nextflow_file_path(file, sample_key)
  base <- sarek_nextflow_base_row(record)
  if (identical(record$input_format, "bam")) {
    return(data.frame(
      base,
      bam = path,
      bai = sarek_nextflow_index_path(file, sample_key, "bai"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  data.frame(
    base,
    cram = path,
    crai = sarek_nextflow_index_path(file, sample_key, "crai"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

sarek_nextflow_vcf_row <- function(record) {
  sample_key <- paste(record$patient, record$sample, sep = "::")
  if (length(record$files) != 1L) {
    stop("Annotation input currently supports one VCF per manifest sample: ", sample_key, ".", call. = FALSE)
  }
  path <- sarek_nextflow_file_path(record$files[[1]], sample_key)
  if (!grepl("\\.vcf\\.gz$", path, ignore.case = TRUE, perl = TRUE)) {
    stop("Annotation input must be a sorted, bgzip-compressed .vcf.gz file: ", path, call. = FALSE)
  }
  data.frame(patient = record$patient, sample = record$sample, vcf = path, stringsAsFactors = FALSE, check.names = FALSE)
}

sarek_build_nextflow_input <- function(manifest) {
  sarek_nextflow_require_helpers()
  records <- sarek_nextflow_manifest_samples(manifest)
  for (record in records) {
    sarek_nextflow_validate_identifier(record$patient, "Patient ID")
    sarek_nextflow_validate_identifier(record$sample, "Sample ID")
  }

  steps <- vapply(records, sarek_nextflow_step_for_sample, character(1))
  if (length(unique(steps)) != 1L) {
    details <- paste(
      vapply(records, function(record) paste(record$patient, record$sample, sep = "::"), character(1)),
      steps,
      sep = " -> ",
      collapse = "; "
    )
    stop(
      "One nf-core/sarek run can use only one starting step. Split this manifest by processing state: ",
      details,
      call. = FALSE
    )
  }
  step <- unique(steps)[[1]]
  analysis_mode <- sarek_text(manifest$analysis$mode)
  sarek_nextflow_validate_mode(records, analysis_mode, step)

  formats <- unique(vapply(records, function(record) record$input_format, character(1)))
  if (length(formats) != 1L) {
    stop(
      "One generated samplesheet currently requires one input format. Split this manifest by format: ",
      paste(formats, collapse = ", "),
      call. = FALSE
    )
  }
  format <- formats[[1]]
  rows <- lapply(records, function(record) {
    if (identical(format, "fastq")) return(sarek_nextflow_fastq_rows(record))
    if (identical(format, "bam") && identical(step, "mapping")) return(sarek_nextflow_unmapped_bam_rows(record))
    if (format %in% c("bam", "cram")) return(sarek_nextflow_alignment_row(record))
    if (identical(format, "vcf")) return(sarek_nextflow_vcf_row(record))
    stop("No samplesheet writer is available for input format: ", format, call. = FALSE)
  })
  samplesheet <- do.call(rbind, rows)
  rownames(samplesheet) <- NULL

  sex_dependent_mode <- analysis_mode %in% c("tumor_only", "matched_tumor_normal")
  missing_sex_patients <- if (!sex_dependent_mode || identical(step, "annotate")) {
    character(0)
  } else {
    unique(vapply(records, function(record) {
      if (identical(record$sex, "NA")) record$patient else ""
    }, character(1)))
  }
  missing_sex_patients <- missing_sex_patients[nzchar(missing_sex_patients)]

  list(
    generator_version = SAREK_NEXTFLOW_INPUT_VERSION,
    pipeline = SAREK_NEXTFLOW_PIPELINE,
    pipeline_version = SAREK_NEXTFLOW_PIPELINE_VERSION,
    manifest_id = sarek_text(manifest$manifest_id),
    analysis_mode = analysis_mode,
    step = step,
    input_format = format,
    samplesheet = samplesheet,
    warnings = if (length(missing_sex_patients)) {
      paste0(
        "Sex chromosomes are not provided for patient(s): ",
        paste(missing_sex_patients, collapse = ", "),
        ". A caller that requires sex will be omitted before submission; the remaining callers can still run."
      )
    } else character(0)
  )
}

sarek_write_nextflow_samplesheet <- function(nextflow_input, path) {
  if (!is.list(nextflow_input) || !is.data.frame(nextflow_input$samplesheet) || !NROW(nextflow_input$samplesheet)) {
    stop("A non-empty Sarek Nextflow input is required.", call. = FALSE)
  }
  path <- sarek_text(path)
  if (!nzchar(path) || !sarek_is_absolute_path(path)) {
    stop("The Nextflow samplesheet output path must be absolute.", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".sarek_samplesheet_", tmpdir = dirname(path), fileext = ".csv")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  utils::write.table(
    nextflow_input$samplesheet,
    file = temporary,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE,
    na = "",
    qmethod = "double",
    fileEncoding = "UTF-8"
  )
  if (!file.rename(temporary, path)) {
    stop("Could not atomically save the Sarek samplesheet: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

sarek_manifest_json_to_nextflow_input <- function(manifest_path) {
  sarek_build_nextflow_input(sarek_nextflow_read_manifest(manifest_path))
}
