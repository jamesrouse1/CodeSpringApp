# Pure helpers for discovering inputs and building a user-confirmed Sarek manifest.
# This module deliberately has no Shiny or pipeline-runtime dependency.

SAREK_MANIFEST_SCHEMA_VERSION <- "1.0"
SAREK_INSPECTOR_VERSION <- "0.2.0"
SAREK_SUPPORTED_INPUT_FORMATS <- c("fastq", "ubam", "bam", "cram", "vcf", "bcf")
SAREK_SAMPLE_ROLES <- c("germline", "tumor", "normal", "unknown")
SAREK_ASSAY_TYPES <- c("WGS", "WES", "targeted", "annotation_only")
SAREK_ANALYSIS_MODES <- c("germline", "tumor_only", "matched_tumor_normal", "annotation_only")
SAREK_SEX_CHROMOSOMES <- c("NA", "XX", "XY")
SAREK_PROCESSING_STATES <- c(
  "unknown",
  "unmapped",
  "aligned",
  "duplicate_marked",
  "recalibrated",
  "analysis_ready",
  "variant_calls"
)

sarek_text <- function(value, default = "") {
  if (is.null(value) || !length(value) || is.na(value[[1]])) return(default)
  value <- trimws(as.character(value[[1]]))
  if (nzchar(value)) value else default
}

sarek_identifier <- function(value, fallback = "sample") {
  value <- gsub("[^A-Za-z0-9_.-]+", "_", sarek_text(value, fallback))
  value <- gsub("^[_.-]+|[_.-]+$", "", value)
  if (nzchar(value)) substr(value, 1L, 128L) else fallback
}

sarek_normalize_sex <- function(value, default = "NA") {
  value <- toupper(sarek_text(value, default))
  if (value %in% c("", "UNKNOWN", "NOT PROVIDED")) value <- "NA"
  value
}

sarek_detect_input_format <- function(path) {
  name <- tolower(basename(sarek_text(path)))
  if (grepl("\\.(fastq|fq)(\\.gz)?$", name, perl = TRUE)) return("fastq")
  if (grepl("\\.ubam$", name, perl = TRUE)) return("ubam")
  if (grepl("\\.bam$", name, perl = TRUE)) return("bam")
  if (grepl("\\.cram$", name, perl = TRUE)) return("cram")
  if (grepl("\\.vcf(\\.gz)?$", name, perl = TRUE)) return("vcf")
  if (grepl("\\.bcf$", name, perl = TRUE)) return("bcf")
  ""
}

sarek_normalize_existing_path <- function(path) {
  path <- sarek_text(path)
  if (!nzchar(path)) stop("An input path is empty.")
  if (!file.exists(path) && !dir.exists(path)) stop("Input path does not exist: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

sarek_is_absolute_path <- function(path) {
  path <- sarek_text(path)
  nzchar(path) && startsWith(path, "/")
}

sarek_path_is_within <- function(path, root) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- sub("/+$", "", normalizePath(root, winslash = "/", mustWork = FALSE))
  identical(path, root) || startsWith(path, paste0(root, "/"))
}

sarek_assert_allowed_paths <- function(paths, allowed_roots = NULL, must_exist = TRUE) {
  if (is.null(allowed_roots) || !length(allowed_roots)) return(invisible(TRUE))
  allowed_roots <- unique(vapply(allowed_roots, sarek_normalize_existing_path, character(1)))
  normalized <- if (isTRUE(must_exist)) {
    unique(vapply(paths, sarek_normalize_existing_path, character(1)))
  } else {
    if (any(!vapply(paths, sarek_is_absolute_path, logical(1)))) stop("Authorized paths must be absolute server paths.")
    unique(vapply(paths, normalizePath, character(1), winslash = "/", mustWork = FALSE))
  }
  outside <- normalized[!vapply(normalized, function(path) {
    any(vapply(allowed_roots, function(root) sarek_path_is_within(path, root), logical(1)))
  }, logical(1))]
  if (length(outside)) {
    stop("Input path is outside the authorized server roots: ", paste(outside, collapse = ", "))
  }
  invisible(TRUE)
}

sarek_find_index <- function(path, format = sarek_detect_input_format(path)) {
  candidates <- switch(
    format,
    bam = c(paste0(path, ".bai"), sub("\\.bam$", ".bai", path, ignore.case = TRUE)),
    cram = c(paste0(path, ".crai"), sub("\\.cram$", ".crai", path, ignore.case = TRUE)),
    vcf = c(paste0(path, ".tbi"), paste0(path, ".csi")),
    bcf = c(paste0(path, ".csi")),
    character(0)
  )
  candidates <- unique(candidates[file.exists(candidates)])
  if (length(candidates)) normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE) else ""
}

sarek_discover_input_paths <- function(paths, recursive = FALSE, max_files = 5000L, allowed_roots = NULL) {
  paths <- unique(vapply(paths, sarek_normalize_existing_path, character(1)))
  sarek_assert_allowed_paths(paths, allowed_roots)
  max_files <- suppressWarnings(as.integer(max_files))
  if (!is.finite(max_files) || max_files < 1L) stop("max_files must be a positive integer.")

  candidates <- character(0)
  for (path in paths) {
    if (dir.exists(path)) {
      candidates <- c(
        candidates,
        list.files(
          path,
          full.names = TRUE,
          recursive = isTRUE(recursive),
          all.files = FALSE,
          include.dirs = FALSE,
          no.. = TRUE
        )
      )
    } else {
      candidates <- c(candidates, path)
    }
  }

  candidates <- unique(candidates[file.exists(candidates) & !dir.exists(candidates)])
  formats <- vapply(candidates, sarek_detect_input_format, character(1))
  supported <- nzchar(formats)
  supported_paths <- unique(vapply(candidates[supported], sarek_normalize_existing_path, character(1)))
  sarek_assert_allowed_paths(supported_paths, allowed_roots)

  if (length(supported_paths) > max_files) {
    stop(
      "Input discovery found ", length(supported_paths),
      " supported files, above the safety limit of ", max_files,
      ". Select a narrower folder or increase the administrator limit."
    )
  }

  list(
    source_paths = paths,
    supported_paths = supported_paths,
    ignored_paths = candidates[!supported]
  )
}

sarek_fastq_technical_tokens <- function(stem) {
  # Some sequencing facilities append a dual index, flowcell, lane, chunk, and
  # read token to the biological sample name, for example:
  # ECO-18-BloodDNA_ACAGTTACCT-GAACTGCCGG_23W33VLT4_L006_001.R1
  #
  # Requiring the dual-index + flowcell + lane chain keeps this conservative:
  # an ordinary sample ending in _L001 is not shortened accidentally.
  pattern <- paste0(
    "(^|[_.-])",
    "([ACGTN]{6,16}-[ACGTN]{6,16})[_.-]",
    "([A-Za-z0-9]{6,20})[_.-]",
    "(L[0-9]{3})(?:[_.-]|$)"
  )
  match <- regexec(pattern, stem, ignore.case = TRUE, perl = TRUE)
  parts <- regmatches(stem, match)[[1]]
  if (length(parts) < 5L) {
    return(list(sample_prefix = "", flowcell = "", lane = ""))
  }

  start <- as.integer(match[[1]][[1]])
  sample_prefix <- if (is.finite(start) && start > 1L) {
    sub("[_.-]+$", "", substr(stem, 1L, start - 1L), perl = TRUE)
  } else {
    ""
  }
  flowcell <- toupper(parts[[4]])
  lane <- toupper(parts[[5]])
  list(
    sample_prefix = sample_prefix,
    flowcell = flowcell,
    lane = paste(flowcell, lane, sep = "_")
  )
}

sarek_filename_tokens <- function(path, format = sarek_detect_input_format(path)) {
  name <- basename(path)
  stem <- name
  stem <- sub("\\.(fastq|fq)(\\.gz)?$", "", stem, ignore.case = TRUE, perl = TRUE)
  stem <- sub("\\.(ubam|bam|cram|vcf|bcf)(\\.gz)?$", "", stem, ignore.case = TRUE, perl = TRUE)

  fastq_technical <- if (identical(format, "fastq")) {
    sarek_fastq_technical_tokens(stem)
  } else {
    list(sample_prefix = "", flowcell = "", lane = "")
  }
  lane <- if (nzchar(fastq_technical$lane)) {
    fastq_technical$lane
  } else if (grepl("(^|[_.-])L[0-9]{3}([_.-]|$)", stem, ignore.case = TRUE, perl = TRUE)) {
    sub(".*(?:^|[_.-])L([0-9]{3})(?:[_.-]|$).*", "L\\1", stem, ignore.case = TRUE, perl = TRUE)
  } else {
    ""
  }

  read <- if (grepl("(^|[_.-])R1([_.-]|$)", stem, ignore.case = TRUE, perl = TRUE)) {
    1L
  } else if (grepl("(^|[_.-])R2([_.-]|$)", stem, ignore.case = TRUE, perl = TRUE)) {
    2L
  } else if (identical(format, "fastq") &&
             grepl("(^|[_.-])1([_.-][0-9]{3})?$", stem, perl = TRUE)) {
    1L
  } else if (identical(format, "fastq") &&
             grepl("(^|[_.-])2([_.-][0-9]{3})?$", stem, perl = TRUE)) {
    2L
  } else {
    NA_integer_
  }

  sample_stem <- fastq_technical$sample_prefix
  if (!nzchar(sample_stem)) {
    sample_stem <- stem
    # Support both ..._L001_R1_001.fastq.gz and
    # ..._L001_001.R1.fastq.gz without removing numeric sample suffixes from
    # filenames that do not contain a lane token.
    for (iteration in seq_len(2L)) {
      sample_stem <- sub("([_.-])R[12]([_.-][0-9]{3})?$", "", sample_stem, ignore.case = TRUE, perl = TRUE)
      if (identical(format, "fastq")) {
        sample_stem <- sub("([_.-])[12]([_.-][0-9]{3})?$", "", sample_stem, perl = TRUE)
      }
      if (nzchar(lane)) {
        sample_stem <- sub("([_.-])[0-9]{3}$", "", sample_stem, perl = TRUE)
      }
      sample_stem <- sub("([_.-])L[0-9]{3}$", "", sample_stem, ignore.case = TRUE, perl = TRUE)
    }
  }
  sample_stem <- sub("([_.-])S[0-9]+$", "", sample_stem, ignore.case = TRUE, perl = TRUE)
  sample_stem <- sarek_identifier(sample_stem)

  role <- "unknown"
  role_confidence <- "none"
  role_basis <- ""
  if (grepl("(^|[_.-])(tumou?r|pdodna|pdo)([_.-]|$)", sample_stem, ignore.case = TRUE, perl = TRUE)) {
    role <- "tumor"
    role_confidence <- "medium"
    role_basis <- if (grepl("(^|[_.-])pdo(dna)?([_.-]|$)", sample_stem, ignore.case = TRUE, perl = TRUE)) {
      "the PDO/PDO-DNA filename token"
    } else {
      "the tumor filename token"
    }
  } else if (grepl("(^|[_.-])(normal|blooddna)([_.-]|$)", sample_stem, ignore.case = TRUE, perl = TRUE)) {
    role <- "normal"
    role_confidence <- "medium"
    role_basis <- if (grepl("(^|[_.-])blooddna([_.-]|$)", sample_stem, ignore.case = TRUE, perl = TRUE)) {
      "the BloodDNA filename token"
    } else {
      "the normal filename token"
    }
  } else if (grepl("(^|[_.-])germline([_.-]|$)", sample_stem, ignore.case = TRUE, perl = TRUE)) {
    role <- "germline"
    role_confidence <- "medium"
    role_basis <- "the germline filename token"
  } else if (grepl("([_.-])T$", sample_stem, ignore.case = TRUE, perl = TRUE)) {
    role <- "tumor"
    role_confidence <- "low"
    role_basis <- "a terminal T filename token"
  } else if (grepl("([_.-])N$", sample_stem, ignore.case = TRUE, perl = TRUE)) {
    role <- "normal"
    role_confidence <- "low"
    role_basis <- "a terminal N filename token"
  }

  patient_stem <- sample_stem
  patient_stem <- gsub(
    "(^|[_.-])(tumou?r|normal|germline|pdodna|pdo|blooddna)([_.-]|$)",
    "_",
    patient_stem,
    ignore.case = TRUE,
    perl = TRUE
  )
  if (role %in% c("tumor", "normal") && identical(role_confidence, "low")) {
    patient_stem <- sub("([_.-])[TN]$", "", patient_stem, ignore.case = TRUE, perl = TRUE)
  }
  patient_stem <- sarek_identifier(patient_stem, sample_stem)

  processing_state <- switch(
    format,
    fastq = "unmapped",
    ubam = "unmapped",
    vcf = "variant_calls",
    bcf = "variant_calls",
    "unknown"
  )

  warning <- if (identical(role, "unknown")) {
    "Sample role requires confirmation."
  } else {
    paste0("Sample role was inferred from ", role_basis, " with ", role_confidence, " confidence; confirm before submission.")
  }

  list(
    sample_id = sample_stem,
    patient_id = patient_stem,
    role = role,
    role_confidence = role_confidence,
    lane = lane,
    read = read,
    processing_state = processing_state,
    warning = warning
  )
}

sarek_build_discovery_table <- function(paths, recursive = FALSE, max_files = 5000L, allowed_roots = NULL) {
  discovery <- sarek_discover_input_paths(
    paths,
    recursive = recursive,
    max_files = max_files,
    allowed_roots = allowed_roots
  )
  if (!length(discovery$supported_paths)) {
    stop("No supported FASTQ, uBAM, BAM, CRAM, VCF, or BCF files were found.")
  }

  rows <- lapply(discovery$supported_paths, function(path) {
    format <- sarek_detect_input_format(path)
    tokens <- sarek_filename_tokens(path, format)
    info <- file.info(path)
    index <- sarek_find_index(path, format)
    data.frame(
      include = TRUE,
      patient_id = tokens$patient_id,
      sex = "NA",
      sample_id = tokens$sample_id,
      role = tokens$role,
      matched_normal_id = "",
      input_format = format,
      processing_state = tokens$processing_state,
      path = normalizePath(path, winslash = "/", mustWork = TRUE),
      index = index,
      inspection_status = if (identical(format, "bam")) "not_run" else "not_applicable",
      processing_recommendation = "",
      processing_confidence = "",
      inspection_summary = if (identical(format, "bam")) "BAM inspection has not run." else "",
      header_sample_ids = "",
      sort_order = "",
      reference_compatibility = "",
      inspection_evidence = "",
      lane = tokens$lane,
      read = tokens$read,
      size_bytes = as.numeric(info$size[[1]]),
      role_confidence = tokens$role_confidence,
      warning = tokens$warning,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  table <- do.call(rbind, rows)
  rownames(table) <- NULL
  for (patient_id in unique(table$patient_id)) {
    patient_rows <- table$patient_id == patient_id
    normal_ids <- unique(table$sample_id[patient_rows & table$role == "normal"])
    if (length(normal_ids) == 1L) {
      table$matched_normal_id[patient_rows & table$role == "tumor"] <- normal_ids[[1]]
    }
  }
  attr(table, "source_paths") <- discovery$source_paths
  attr(table, "ignored_count") <- length(discovery$ignored_paths)
  table
}

sarek_validate_confirmation_table <- function(table) {
  errors <- character(0)
  warnings <- character(0)
  required <- c(
    "include",
    "patient_id",
    "sample_id",
    "role",
    "input_format",
    "processing_state",
    "path"
  )

  missing <- setdiff(required, names(table))
  if (length(missing)) {
    return(list(
      valid = FALSE,
      errors = paste("Missing required columns:", paste(missing, collapse = ", ")),
      warnings = warnings
    ))
  }

  include <- suppressWarnings(as.logical(table$include))
  if (anyNA(include)) errors <- c(errors, "Every include value must be explicitly TRUE or FALSE.")
  selected <- table[which(!is.na(include) & include), , drop = FALSE]
  if (!NROW(selected)) errors <- c(errors, "At least one input file must be included.")
  if (!NROW(selected)) return(list(valid = FALSE, errors = errors, warnings = warnings))

  identifier_pattern <- "^[A-Za-z0-9][A-Za-z0-9_.-]*$"
  patient_ids <- trimws(as.character(selected$patient_id))
  sample_ids <- trimws(as.character(selected$sample_id))
  bad_patient <- is.na(patient_ids) | !nzchar(patient_ids) | !grepl(identifier_pattern, patient_ids, perl = TRUE)
  bad_sample <- is.na(sample_ids) | !nzchar(sample_ids) | !grepl(identifier_pattern, sample_ids, perl = TRUE)
  if (any(bad_patient)) errors <- c(errors, "Patient IDs may contain only letters, numbers, periods, underscores, and hyphens.")
  if (any(bad_sample)) errors <- c(errors, "Sample IDs may contain only letters, numbers, periods, underscores, and hyphens.")
  sex <- if ("sex" %in% names(selected)) {
    vapply(selected$sex, sarek_normalize_sex, character(1))
  } else {
    rep("NA", NROW(selected))
  }
  if (any(!sex %in% SAREK_SEX_CHROMOSOMES)) {
    errors <- c(errors, "Sex chromosomes must be XX, XY, or Not provided.")
  }
  sex_by_patient <- split(sex, patient_ids)
  inconsistent_sex <- names(sex_by_patient)[vapply(
    sex_by_patient,
    function(values) length(unique(values)) > 1L,
    logical(1)
  )]
  if (length(inconsistent_sex)) {
    errors <- c(errors, paste0(
      "Sex chromosomes are patient-level metadata and must agree across all samples for: ",
      paste(inconsistent_sex, collapse = ", ")
    ))
  }
  if (any(!selected$role %in% SAREK_SAMPLE_ROLES)) errors <- c(errors, "One or more sample roles are invalid.")
  if (any(!selected$input_format %in% SAREK_SUPPORTED_INPUT_FORMATS)) errors <- c(errors, "One or more input formats are invalid.")
  if (any(!selected$processing_state %in% SAREK_PROCESSING_STATES)) errors <- c(errors, "One or more processing states are invalid.")
  if (anyDuplicated(selected$path)) errors <- c(errors, "The same input path is included more than once.")
  if (any(!file.exists(selected$path))) errors <- c(errors, "One or more selected input files no longer exist.")
  if (any(file.access(selected$path, mode = 4) != 0)) errors <- c(errors, "One or more selected input files are not readable.")

  sample_key <- paste(selected$patient_id, selected$sample_id, sep = "::")
  formats_by_sample <- split(selected$input_format, sample_key)
  roles_by_sample <- split(selected$role, sample_key)
  states_by_sample <- split(selected$processing_state, sample_key)
  mixed_formats <- names(formats_by_sample)[vapply(formats_by_sample, function(x) length(unique(x)) > 1L, logical(1))]
  mixed_roles <- names(roles_by_sample)[vapply(roles_by_sample, function(x) length(unique(x)) > 1L, logical(1))]
  mixed_states <- names(states_by_sample)[vapply(states_by_sample, function(x) length(unique(x)) > 1L, logical(1))]
  if (length(mixed_formats)) errors <- c(errors, "A sample cannot mix different input formats.")
  if (length(mixed_roles)) errors <- c(errors, "All files assigned to one sample must use the same sample role.")
  if (length(mixed_states)) errors <- c(errors, "All files assigned to one sample must use the same processing state.")

  fastq <- selected[selected$input_format == "fastq", , drop = FALSE]
  if (NROW(fastq)) {
    if (!all(c("lane", "read") %in% names(fastq))) {
      errors <- c(errors, "FASTQ inputs require lane and read columns.")
    } else {
      lane_values <- trimws(as.character(fastq$lane))
      lane <- ifelse(is.na(lane_values) | !nzchar(lane_values), "unlabelled", lane_values)
      pair_key <- paste(fastq$patient_id, fastq$sample_id, lane, sep = "::")
      reads_by_pair <- split(fastq$read, pair_key)
      incomplete <- names(reads_by_pair)[!vapply(reads_by_pair, function(reads) {
        reads <- sort(suppressWarnings(as.integer(reads)))
        length(reads) == 2L && !anyNA(reads) && identical(reads, c(1L, 2L))
      }, logical(1))]
      if (length(incomplete)) {
        errors <- c(errors, paste0(
          "Each FASTQ sample lane must contain exactly one R1 and one R2: ",
          paste(incomplete, collapse = ", ")
        ))
      }
    }
  }

  if ("matched_normal_id" %in% names(selected)) {
    matches_by_sample <- split(trimws(as.character(selected$matched_normal_id)), sample_key)
    inconsistent_matches <- names(matches_by_sample)[vapply(matches_by_sample, function(x) {
      x <- x[!is.na(x) & nzchar(x)]
      length(unique(x)) > 1L
    }, logical(1))]
    if (length(inconsistent_matches)) {
      errors <- c(errors, "All files assigned to one tumor sample must name the same matched normal.")
    }
  }

  if (any(selected$role == "unknown", na.rm = TRUE)) warnings <- c(warnings, "Every unknown sample role must be confirmed before pipeline submission.")
  index_missing <- if ("index" %in% names(selected)) {
    is.na(selected$index) | !nzchar(trimws(as.character(selected$index)))
  } else {
    rep(TRUE, NROW(selected))
  }
  missing_index <- selected$input_format %in% c("bam", "cram", "vcf", "bcf") &
    index_missing
  missing_alignment_index <- selected$input_format %in% c("bam", "cram") & index_missing
  if (any(missing_alignment_index)) {
    errors <- c(errors, "Every included BAM or CRAM must have a detected companion index before confirmation.")
  }
  alignment_with_index <- selected$input_format %in% c("bam", "cram") & !index_missing
  unreadable_alignment_index <- rep(FALSE, NROW(selected))
  if (any(alignment_with_index)) {
    index_paths <- as.character(selected$index[alignment_with_index])
    unreadable_alignment_index[alignment_with_index] <- !file.exists(index_paths) |
      file.access(index_paths, mode = 4) != 0
  }
  if (any(unreadable_alignment_index)) {
    errors <- c(errors, "One or more detected BAM or CRAM indexes are missing or unreadable.")
  }
  if (any(missing_index & !missing_alignment_index)) {
    warnings <- c(warnings, "One or more indexed variant formats do not currently have a detected index.")
  }

  alignment_rows <- selected$input_format %in% c("bam", "cram")
  unknown_alignment_state <- alignment_rows & selected$processing_state == "unknown"
  if (any(unknown_alignment_state)) {
    errors <- c(errors, "Every included BAM or CRAM processing state must be explicitly confirmed before submission.")
  }
  if ("inspection_status" %in% names(selected)) {
    failed_bam_inspection <- selected$input_format == "bam" & selected$inspection_status == "failed"
    if (any(failed_bam_inspection, na.rm = TRUE)) {
      errors <- c(errors, "One or more BAM files failed lightweight integrity or index inspection.")
    }
    unverified_bam <- selected$input_format == "bam" & selected$inspection_status %in% c(
      "not_run", "unavailable", "deferred"
    )
    if (any(unverified_bam, na.rm = TRUE)) {
      warnings <- c(warnings, "One or more BAM files could not be automatically inspected; review their processing provenance manually.")
    }
    if ("reference_compatibility" %in% names(selected)) {
      reference_review <- selected$input_format == "bam" & grepl(
        "may not match GATK[.]GRCh38|not the GRCh38 length|could not be inferred",
        as.character(selected$reference_compatibility),
        ignore.case = TRUE,
        perl = TRUE
      )
      if (any(reference_review, na.rm = TRUE)) {
        warnings <- c(warnings, "One or more BAM headers require manual review for GATK.GRCh38 reference compatibility.")
      }
    }
  }

  list(valid = !length(errors), errors = unique(errors), warnings = unique(warnings))
}

sarek_validate_analysis_mode <- function(table, analysis_mode) {
  include <- suppressWarnings(as.logical(table$include))
  selected <- table[which(!is.na(include) & include), , drop = FALSE]
  errors <- character(0)
  roles <- unique(selected$role)
  formats <- unique(selected$input_format)

  analysis_mode <- sarek_text(analysis_mode)
  if (!analysis_mode %in% SAREK_ANALYSIS_MODES) {
    return("The analysis mode is invalid.")
  }

  if (identical(analysis_mode, "germline") && any(!roles %in% "germline")) {
    errors <- c(errors, "Germline mode requires every included sample to have the germline role.")
  }
  if (identical(analysis_mode, "tumor_only") && any(!roles %in% "tumor")) {
    errors <- c(errors, "Tumor-only mode requires every included sample to have the tumor role.")
  }
  if (identical(analysis_mode, "matched_tumor_normal")) {
    by_patient <- split(selected$role, selected$patient_id)
    incomplete <- names(by_patient)[!vapply(
      by_patient,
      function(x) any(x == "tumor") && any(x == "normal"),
      logical(1)
    )]
    if (length(incomplete)) {
      errors <- c(errors, paste0(
        "Matched tumor-normal mode requires at least one tumor and one normal for each patient: ",
        paste(incomplete, collapse = ", ")
      ))
    }

    if (!"matched_normal_id" %in% names(selected)) {
      errors <- c(errors, "Matched tumor-normal mode requires an explicit matched_normal_id column.")
    } else {
      sample_columns <- c("patient_id", "sample_id", "role", "matched_normal_id")
      sample_meta <- unique(selected[, sample_columns, drop = FALSE])
      tumors <- sample_meta[sample_meta$role == "tumor", , drop = FALSE]
      normals <- sample_meta[sample_meta$role == "normal", , drop = FALSE]
      used_normals <- character(0)

      for (i in seq_len(NROW(tumors))) {
        matched_normal <- sarek_text(tumors$matched_normal_id[[i]])
        if (!nzchar(matched_normal)) {
          errors <- c(errors, paste0("Tumor sample ", tumors$sample_id[[i]], " requires a matched normal."))
          next
        }
        valid_normal <- normals$patient_id == tumors$patient_id[[i]] & normals$sample_id == matched_normal
        if (!any(valid_normal)) {
          errors <- c(errors, paste0(
            "Tumor sample ", tumors$sample_id[[i]],
            " names a normal that is not included for the same patient: ", matched_normal
          ))
        } else {
          used_normals <- c(used_normals, paste(tumors$patient_id[[i]], matched_normal, sep = "::"))
        }
      }

      normal_keys <- paste(normals$patient_id, normals$sample_id, sep = "::")
      unused_normals <- normals$sample_id[!normal_keys %in% used_normals]
      if (length(unused_normals)) {
        errors <- c(errors, paste0(
          "Included normal samples must be explicitly paired to at least one tumor: ",
          paste(unique(unused_normals), collapse = ", ")
        ))
      }
    }
  }
  if (identical(analysis_mode, "annotation_only") && any(!formats %in% c("vcf", "bcf"))) {
    errors <- c(errors, "Annotation-only mode accepts only VCF or BCF inputs.")
  }
  unique(errors)
}

sarek_manifest_file_record <- function(row) {
  result <- list(path = as.character(row$path[[1]]))
  if ("index" %in% names(row) && nzchar(sarek_text(row$index[[1]]))) result$index <- as.character(row$index[[1]])
  if ("lane" %in% names(row) && nzchar(sarek_text(row$lane[[1]]))) result$lane <- as.character(row$lane[[1]])
  if ("read" %in% names(row) && !is.na(row$read[[1]])) result$read <- as.integer(row$read[[1]])
  result
}

sarek_build_manifest <- function(
  confirmation_table,
  manifest_id,
  created_by,
  assay_type,
  analysis_mode,
  preset,
  results_root,
  work_root,
  source_paths = attr(confirmation_table, "source_paths"),
  species = "human",
  assembly = "GRCh38",
  sarek_genome = "GATK.GRCh38",
  allowed_results_roots = NULL,
  allowed_work_roots = NULL
) {
  validation <- sarek_validate_confirmation_table(confirmation_table)
  mode_errors <- sarek_validate_analysis_mode(confirmation_table, analysis_mode)
  errors <- unique(c(validation$errors, mode_errors))
  if (!sarek_text(assay_type) %in% SAREK_ASSAY_TYPES) errors <- c(errors, "The assay type is invalid.")
  if (!nzchar(sarek_text(preset))) errors <- c(errors, "The analysis preset is required.")
  if (!nzchar(sarek_text(species)) || !nzchar(sarek_text(assembly)) || !nzchar(sarek_text(sarek_genome))) {
    errors <- c(errors, "Species, reference assembly, and Sarek genome are required.")
  }
  if (!sarek_is_absolute_path(results_root) || !sarek_is_absolute_path(work_root)) {
    errors <- c(errors, "Results and work roots must be absolute server paths.")
  }
  if (length(errors)) stop(paste(errors, collapse = "\n"))

  include <- suppressWarnings(as.logical(confirmation_table$include))
  selected <- confirmation_table[which(!is.na(include) & include), , drop = FALSE]
  results_root <- normalizePath(results_root, winslash = "/", mustWork = FALSE)
  work_root <- normalizePath(work_root, winslash = "/", mustWork = FALSE)
  sarek_assert_allowed_paths(results_root, allowed_results_roots, must_exist = FALSE)
  sarek_assert_allowed_paths(work_root, allowed_work_roots, must_exist = FALSE)

  patient_ids <- unique(selected$patient_id)
  patients <- lapply(patient_ids, function(patient_id) {
    patient_rows <- selected[selected$patient_id == patient_id, , drop = FALSE]
    patient_sex <- if ("sex" %in% names(patient_rows)) {
      sarek_normalize_sex(patient_rows$sex[[1]])
    } else {
      "NA"
    }
    sample_ids <- unique(patient_rows$sample_id)
    samples <- lapply(sample_ids, function(sample_id) {
      sample_rows <- patient_rows[patient_rows$sample_id == sample_id, , drop = FALSE]
      notes <- if ("warning" %in% names(sample_rows)) {
        note_values <- as.character(sample_rows$warning)
        unique(note_values[!is.na(note_values) & nzchar(note_values)])
      } else {
        character(0)
      }
      sample <- list(
        sample_id = sample_id,
        role = unique(sample_rows$role)[[1]],
        input_format = unique(sample_rows$input_format)[[1]],
        processing_state = unique(sample_rows$processing_state)[[1]],
        files = lapply(seq_len(NROW(sample_rows)), function(i) sarek_manifest_file_record(sample_rows[i, , drop = FALSE]))
      )
      if (length(notes)) sample$confirmation_notes <- notes
      sample
    })

    relationships <- list()
    if (identical(analysis_mode, "matched_tumor_normal")) {
      tumors <- unique(patient_rows$sample_id[patient_rows$role == "tumor"])
      for (tumor in tumors) {
        tumor_rows <- patient_rows[patient_rows$sample_id == tumor, , drop = FALSE]
        matched_normal <- unique(trimws(as.character(tumor_rows$matched_normal_id)))
        matched_normal <- matched_normal[nzchar(matched_normal)][[1]]
        relationships[[length(relationships) + 1L]] <- list(
          type = "matched_tumor_normal",
          sample_ids = c(tumor, matched_normal)
        )
      }
    }

    list(
      patient_id = patient_id,
      sex = patient_sex,
      samples = samples,
      relationships = relationships
    )
  })

  source_paths <- unique(vapply(
    if (is.null(source_paths) || !length(source_paths)) selected$path else source_paths,
    function(path) normalizePath(path, winslash = "/", mustWork = FALSE),
    character(1)
  ))

  list(
    schema_version = SAREK_MANIFEST_SCHEMA_VERSION,
    manifest_id = sarek_identifier(manifest_id, "sarek_analysis"),
    status = "confirmed",
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    created_by = sarek_text(created_by, "unknown"),
    reference = list(
      species = sarek_text(species),
      assembly = sarek_text(assembly),
      sarek_genome = sarek_text(sarek_genome)
    ),
    assay = list(
      type = sarek_text(assay_type),
      intervals = NULL
    ),
    analysis = list(
      mode = sarek_text(analysis_mode),
      preset = sarek_text(preset, "core")
    ),
    storage = list(
      results_root = results_root,
      work_root = work_root
    ),
    patients = patients,
    provenance = list(
      source_paths = source_paths,
      inspector_version = SAREK_INSPECTOR_VERSION
    )
  )
}

sarek_write_manifest <- function(manifest, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The jsonlite R package is required to write a Sarek manifest.")
  }
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".sarek_manifest_", tmpdir = dirname(path), fileext = ".json")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  jsonlite::write_json(
    manifest,
    path = temporary,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    na = "null"
  )
  if (!file.rename(temporary, path)) stop("Could not atomically save the Sarek manifest: ", path)
  path
}
