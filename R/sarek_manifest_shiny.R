# Shiny adapter for the pure Sarek manifest helpers in sarek_manifest.R.
# Keeping UI state here avoids adding pipeline-specific logic to app.R.

sarek_shiny_value <- function(value, default = "") {
  if (is.null(value) || !length(value) || is.na(value[[1]])) return(default)
  value
}

sarek_parse_path_input <- function(value) {
  lines <- unlist(strsplit(as.character(sarek_shiny_value(value)), "\n", fixed = TRUE), use.names = FALSE)
  paths <- trimws(lines)
  unique(paths[!is.na(paths) & nzchar(paths)])
}

sarek_parse_include_value <- function(value) {
  value <- tolower(trimws(as.character(sarek_shiny_value(value))))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  NA
}

sarek_editable_columns <- function() {
  c(
    "include",
    "patient_id",
    "sex",
    "sample_id",
    "role",
    "matched_normal_id",
    "input_format",
    "processing_state",
    "lane",
    "read"
  )
}

sarek_manifest_column_widths <- function() {
  c(
    include = "80px",
    patient_id = "150px",
    sex = "120px",
    sample_id = "170px",
    role = "110px",
    matched_normal_id = "180px",
    input_format = "120px",
    processing_state = "170px",
    path = "360px",
    index = "320px",
    lane = "90px",
    read = "80px",
    size_bytes = "120px",
    role_confidence = "140px",
    warning = "340px"
  )
}

sarek_manifest_column_defs <- function(table) {
  widths <- sarek_manifest_column_widths()
  widths <- widths[names(widths) %in% names(table)]
  lapply(names(widths), function(column) {
    list(targets = match(column, names(table)) - 1L, width = unname(widths[[column]]))
  })
}

sarek_manifest_status_kind <- function(message) {
  message <- trimws(as.character(sarek_shiny_value(message)))
  if (startsWith(message, "ERROR:")) return("error")
  if (startsWith(message, "ACTION REQUIRED:")) return("review")
  "info"
}

sarek_sample_key <- function(patient_id, sample_id) {
  paste(as.character(patient_id), as.character(sample_id), sep = "::")
}

sarek_fastq_pairing_status <- function(rows) {
  include <- suppressWarnings(as.logical(rows$include))
  rows <- rows[which(!is.na(include) & include), , drop = FALSE]
  if (!NROW(rows)) return("Excluded")
  rows <- rows[rows$input_format == "fastq", , drop = FALSE]
  if (!NROW(rows)) return("Not applicable")

  lane_values <- trimws(as.character(rows$lane))
  lanes <- ifelse(is.na(lane_values) | !nzchar(lane_values), "unlabelled", lane_values)
  reads_by_lane <- split(rows$read, lanes)
  problems <- unlist(lapply(names(reads_by_lane), function(lane) {
    reads <- suppressWarnings(as.integer(reads_by_lane[[lane]]))
    if (length(reads) == 2L && !anyNA(reads) && identical(sort(reads), c(1L, 2L))) {
      return(character(0))
    }
    details <- character(0)
    missing <- setdiff(c(1L, 2L), reads[!is.na(reads)])
    if (length(missing)) details <- c(details, paste0("missing R", missing))
    if (anyNA(reads)) details <- c(details, "read number not detected")
    duplicated_reads <- unique(reads[!is.na(reads) & duplicated(reads)])
    if (length(duplicated_reads)) details <- c(details, paste0("duplicate R", duplicated_reads))
    if (length(reads) > 2L && !length(duplicated_reads)) details <- c(details, "more than two files")
    paste0(lane, " ", paste(details, collapse = ", "))
  }), use.names = FALSE)

  if (length(problems)) paste0("Needs attention: ", paste(problems, collapse = "; ")) else {
    paste0("Complete: ", length(reads_by_lane), " lane", if (length(reads_by_lane) == 1L) "" else "s")
  }
}

sarek_sample_index_display <- function(rows) {
  indexed <- rows$input_format %in% c("bam", "cram", "vcf", "bcf")
  if (!any(indexed)) return("Not required")
  if (!"index" %in% names(rows)) return("Missing")
  indexes <- trimws(as.character(rows$index[indexed]))
  present <- !is.na(indexes) & nzchar(indexes)
  if (!all(present)) return("Missing")
  if (length(indexes) == 1L) return(paste0("Detected: ", basename(indexes[[1]])))
  paste0("Detected: ", length(indexes), " indexes")
}

sarek_sample_bam_inspection_display <- function(rows) {
  bam <- rows[rows$input_format == "bam", , drop = FALSE]
  if (!NROW(bam)) return("Not applicable")
  if (!"inspection_status" %in% names(bam)) return("Not run")
  statuses <- unique(trimws(as.character(bam$inspection_status)))
  statuses <- statuses[!is.na(statuses) & nzchar(statuses)]
  recommendations <- if ("processing_recommendation" %in% names(bam)) {
    unique(trimws(as.character(bam$processing_recommendation)))
  } else character(0)
  recommendations <- recommendations[!is.na(recommendations) & nzchar(recommendations) & recommendations != "unknown"]
  confidences <- if ("processing_confidence" %in% names(bam)) {
    unique(trimws(as.character(bam$processing_confidence)))
  } else character(0)
  confidences <- confidences[!is.na(confidences) & nzchar(confidences) & confidences != "none"]
  status <- if (length(statuses)) statuses[[1]] else "not_run"
  label <- gsub("_", " ", status, fixed = TRUE)
  if (identical(status, "passed") && length(recommendations)) {
    return(paste0(
      "Passed: ", gsub("_", " ", recommendations[[1]], fixed = TRUE),
      if (length(confidences)) paste0(" (", confidences[[1]], ")") else ""
    ))
  }
  paste0(toupper(substr(label, 1L, 1L)), substring(label, 2L))
}

sarek_sample_review_table <- function(table) {
  if (is.null(table) || !NROW(table)) return(data.frame())
  keys <- sarek_sample_key(table$patient_id, table$sample_id)
  groups <- split(seq_len(NROW(table)), keys)
  rows <- lapply(groups, function(index) {
    sample_rows <- table[index, , drop = FALSE]
    values <- function(column, blank = "") {
      value <- unique(trimws(as.character(sample_rows[[column]])))
      value <- value[!is.na(value) & nzchar(value)]
      if (length(value)) paste(value, collapse = ", ") else blank
    }
    data.frame(
      include = all(suppressWarnings(as.logical(sample_rows$include)), na.rm = TRUE),
      patient_id = as.character(sample_rows$patient_id[[1]]),
      sex = if ("sex" %in% names(sample_rows)) sarek_normalize_sex(sample_rows$sex[[1]]) else "NA",
      sample_id = as.character(sample_rows$sample_id[[1]]),
      role = values("role", "unknown"),
      matched_normal_id = values("matched_normal_id"),
      input_format = values("input_format"),
      processing_state = values("processing_state"),
      index = sarek_sample_index_display(sample_rows),
      BAM_inspection = sarek_sample_bam_inspection_display(sample_rows),
      file_count = NROW(sample_rows),
      FASTQ_pairing = sarek_fastq_pairing_status(sample_rows),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

sarek_sample_review_display_table <- function(table) {
  review <- sarek_sample_review_table(table)
  required <- c(
    "include", "patient_id", "sex", "sample_id", "role", "matched_normal_id",
    "input_format", "processing_state", "index", "BAM_inspection", "file_count", "FASTQ_pairing"
  )
  if (!NROW(review)) {
    return(data.frame(
      Include = character(0), `Patient ID` = character(0), `Sample ID` = character(0),
      `Sex chromosomes` = character(0), Role = character(0), `Matched normal` = character(0), `Input format` = character(0),
      `Processing state` = character(0), Index = character(0), `BAM inspection` = character(0),
      Files = integer(0), `FASTQ pairing` = character(0),
      stringsAsFactors = FALSE, check.names = FALSE
    ))
  }
  missing <- setdiff(required, names(review))
  if (length(missing)) {
    stop("The sample review table is missing required fields: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  data.frame(
    Include = ifelse(!is.na(review$include) & as.logical(review$include), "☑", "☐"),
    `Patient ID` = as.character(review$patient_id),
    `Sex chromosomes` = ifelse(review$sex == "NA", "Not provided", as.character(review$sex)),
    `Sample ID` = as.character(review$sample_id),
    Role = as.character(review$role),
    `Matched normal` = as.character(review$matched_normal_id),
    `Input format` = as.character(review$input_format),
    `Processing state` = as.character(review$processing_state),
    Index = as.character(review$index),
    `BAM inspection` = as.character(review$BAM_inspection),
    Files = as.integer(review$file_count),
    `FASTQ pairing` = as.character(review$FASTQ_pairing),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

sarek_run_history_label <- function(run) {
  run_id <- sarek_text(run[["run_id"]], "unnamed run")
  submitted <- sarek_text(run[["submitted_at"]])
  if (nzchar(submitted)) paste(run_id, submitted, sep = " — ") else run_id
}

sarek_recommend_analysis_mode <- function(table) {
  if (is.null(table) || !NROW(table)) return("")
  include <- suppressWarnings(as.logical(table$include))
  selected <- table[which(!is.na(include) & include), , drop = FALSE]
  if (!NROW(selected)) return("")
  formats <- unique(as.character(selected$input_format))
  roles <- unique(as.character(selected$role))
  if (all(formats %in% c("vcf", "bcf"))) return("annotation_only")
  if (length(roles) == 1L && identical(roles, "germline")) return("germline")
  if (length(roles) == 1L && identical(roles, "tumor")) return("tumor_only")
  if (all(c("tumor", "normal") %in% roles)) {
    by_patient <- split(selected$role, selected$patient_id)
    complete <- vapply(by_patient, function(x) any(x == "tumor") && any(x == "normal"), logical(1))
    if (all(complete)) return("matched_tumor_normal")
  }
  ""
}

sarek_analysis_mode_label <- function(mode) {
  labels <- c(
    germline = "Germline",
    tumor_only = "Tumor only",
    matched_tumor_normal = "Matched tumor-normal",
    annotation_only = "Annotation only"
  )
  if (mode %in% names(labels)) unname(labels[[mode]]) else "No clear recommendation"
}

sarek_run_activity_panel <- function(activity, state = "") {
  state <- sarek_normalize_slurm_state(state)
  if (!is.list(activity) || !isTRUE(activity$available)) {
    message <- if (state %in% c("PENDING", "CONFIGURING", "SUBMITTED")) {
      "The Nextflow controller is starting. Pipeline activity will appear shortly."
    } else if (state %in% c("RUNNING", "COMPLETING", "SUSPENDED")) {
      "The controller is active, but no readable Nextflow activity has been recorded yet."
    } else {
      "No Nextflow activity log is available for this run."
    }
    return(shiny::div(
      class = "sarek-live-activity sarek-live-activity-empty",
      shiny::h4("Live Nextflow activity"),
      shiny::tags$p(message)
    ))
  }

  current <- as.character(activity$current_labels)
  current <- current[nzchar(current)]
  current_text <- if (length(current)) {
    paste(current, collapse = ", ")
  } else if (state %in% c("RUNNING", "COMPLETING", "SUSPENDED")) {
    "Nextflow is starting or is between task submissions."
  } else {
    "No task is currently running."
  }
  estimate <- sarek_text(activity$estimate$label, "A current-task estimate is not available yet.")
  events <- as.character(activity$recent_events)
  events <- events[nzchar(events)]

  shiny::div(
    class = "sarek-live-activity",
    shiny::h4("Live Nextflow activity"),
    shiny::div(
      class = "sarek-live-task",
      shiny::tags$strong("Current step: "),
      current_text
    ),
    shiny::div(
      class = "sarek-live-counts",
      shiny::span(shiny::tags$strong(activity$running), " running"),
      shiny::span(shiny::tags$strong(activity$completed), " completed"),
      shiny::span(shiny::tags$strong(activity$failed), " failed")
    ),
    shiny::tags$p(class = "sarek-live-estimate", estimate),
    shiny::tags$p(
      class = "muted small-note",
      "Any time estimate applies only to the current task and is based on completed tasks of the same step; it is not a whole-workflow ETA."
    ),
    if (length(events)) {
      shiny::tags$details(
        class = "sarek-live-log",
        shiny::tags$summary("Recent Nextflow events"),
        shiny::tags$pre(paste(events, collapse = "\n"))
      )
    }
  )
}

sarek_apply_sample_update <- function(
  table,
  sample_key,
  include,
  patient_id,
  sample_id,
  role,
  matched_normal_id,
  processing_state,
  sex = "NA"
) {
  keys <- sarek_sample_key(table$patient_id, table$sample_id)
  rows <- which(keys == sample_key)
  if (!length(rows)) stop("The selected sample is no longer available. Select it again.")
  table$include[rows] <- isTRUE(include)
  table$patient_id[rows] <- trimws(as.character(patient_id))
  if (!"sex" %in% names(table)) table$sex <- "NA"
  sex <- sarek_normalize_sex(sex)
  if (!sex %in% SAREK_SEX_CHROMOSOMES) stop("Sex chromosomes must be XX, XY, or Not provided.")
  table$sex[table$patient_id == trimws(as.character(patient_id))] <- sex
  table$sample_id[rows] <- trimws(as.character(sample_id))
  table$role[rows] <- as.character(role)
  table$matched_normal_id[rows] <- if (identical(role, "tumor")) trimws(as.character(matched_normal_id)) else ""
  table$processing_state[rows] <- as.character(processing_state)
  table
}

sarek_apply_file_pairing_update <- function(table, path, lane, read) {
  rows <- which(as.character(table$path) == as.character(path))
  if (length(rows) != 1L) stop("The selected file is no longer available. Select it again.")
  if (!identical(as.character(table$input_format[[rows]]), "fastq")) {
    stop("Lane and read corrections apply only to FASTQ files.")
  }
  read <- suppressWarnings(as.integer(read))
  if (!is.na(read) && !read %in% c(1L, 2L)) stop("FASTQ read must be R1, R2, or not detected.")
  table$lane[[rows]] <- trimws(as.character(lane))
  table$read[[rows]] <- read
  table
}

sarek_apply_confirmation_edit <- function(table, edit) {
  if (is.null(table) || !NROW(table) || is.null(edit$row) || is.null(edit$col)) return(table)
  row <- suppressWarnings(as.integer(edit$row))
  column <- suppressWarnings(as.integer(edit$col)) + 1L
  if (
    is.na(row) || is.na(column) ||
    row < 1L || row > NROW(table) ||
    column < 1L || column > NCOL(table)
  ) {
    return(table)
  }

  column_name <- names(table)[[column]]
  if (!column_name %in% sarek_editable_columns()) return(table)
  value <- as.character(sarek_shiny_value(edit$value))

  if (identical(column_name, "include")) {
    table[[column_name]][[row]] <- sarek_parse_include_value(value)
  } else if (identical(column_name, "read")) {
    value <- trimws(value)
    table[[column_name]][[row]] <- if (nzchar(value)) suppressWarnings(as.integer(value)) else NA_integer_
  } else {
    table[[column_name]][[row]] <- trimws(value)
  }
  table
}

sarek_confirmation_summary <- function(table) {
  if (is.null(table) || !NROW(table)) {
    return(data.frame(
      Measure = c("Files", "Included files", "Patients", "Samples"),
      Value = c(0L, 0L, 0L, 0L),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  include <- suppressWarnings(as.logical(table$include))
  selected <- table[which(!is.na(include) & include), , drop = FALSE]
  sample_keys <- if (NROW(selected)) {
    unique(paste(selected$patient_id, selected$sample_id, sep = "::"))
  } else {
    character(0)
  }
  formats <- if (NROW(selected)) paste(sort(unique(selected$input_format)), collapse = ", ") else ""
  roles <- if (NROW(selected)) paste(sort(unique(selected$role)), collapse = ", ") else ""
  data.frame(
    Measure = c("Files", "Included files", "Patients", "Samples", "Formats", "Roles"),
    Value = c(
      as.character(NROW(table)),
      as.character(NROW(selected)),
      as.character(length(unique(selected$patient_id))),
      as.character(length(sample_keys)),
      formats,
      roles
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

sarek_manifest_validation_text <- function(validation, mode_errors = character(0)) {
  errors <- unique(c(sarek_shiny_value(validation$errors, character(0)), mode_errors))
  warnings <- unique(sarek_shiny_value(validation$warnings, character(0)))
  sections <- character(0)
  if (length(errors)) {
    sections <- c(sections, "ERRORS", paste0("- ", errors))
  } else {
    sections <- c(sections, "VALID", "- The current confirmation table and analysis mode are valid.")
  }
  if (length(warnings)) {
    sections <- c(sections, "", "WARNINGS", paste0("- ", warnings))
  }
  paste(sections, collapse = "\n")
}

sarek_manifest_validation_kind <- function(message) {
  message <- as.character(sarek_shiny_value(message))
  if (grepl("(^|\\n)ERRORS($|\\n)", message, perl = TRUE)) return("error")
  if (grepl("(^|\\n)WARNINGS($|\\n)", message, perl = TRUE)) return("warning")
  if (grepl("(^|\\n)(VALID|CONFIRMED)($|\\n)", message, perl = TRUE)) return("success")
  "info"
}

sarek_manifest_validation_heading <- function(kind) {
  switch(
    kind,
    error = "Validation failed — corrections required",
    warning = "Validation passed with warnings — review required",
    success = "Validation passed",
    "Validation has not run"
  )
}

sarek_confirmed_samples_table <- function(manifest) {
  if (is.null(manifest) || !is.list(manifest) || !length(manifest$patients)) return(data.frame())
  rows <- list()
  for (patient in manifest$patients) {
    patient_id <- sarek_text(patient$patient_id)
    matched_normals <- list()
    for (relationship in patient$relationships) {
      if (identical(sarek_text(relationship$type), "matched_tumor_normal") && length(relationship$sample_ids) >= 2L) {
        matched_normals[[sarek_text(relationship$sample_ids[[1]])]] <- sarek_text(relationship$sample_ids[[2]])
      }
    }
    for (sample in patient$samples) {
      sample_id <- sarek_text(sample$sample_id)
      rows[[length(rows) + 1L]] <- data.frame(
        Patient = patient_id,
        `Sex chromosomes` = if (sarek_normalize_sex(patient$sex) == "NA") "Not provided" else sarek_normalize_sex(patient$sex),
        Sample = sample_id,
        Role = sarek_text(sample$role),
        `Matched normal` = sarek_text(matched_normals[[sample_id]]),
        Input = toupper(sarek_text(sample$input_format)),
        `Processing stage` = gsub("_", " ", sarek_text(sample$processing_state), fixed = TRUE),
        Files = length(sample$files),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

sarek_manifest_field_guide <- function() {
  shiny::tags$details(
    class = "sarek-field-guide",
    open = "open",
    shiny::tags$summary("Field definitions and rules"),
    shiny::fluidRow(
      shiny::column(
        4,
        shiny::h5("Sample identity"),
        shiny::tags$dl(
          shiny::tags$dt("Include"),
          shiny::tags$dd("Checkbox controlling whether the entire sample and all its files enter the manifest."),
          shiny::tags$dt("Patient ID"),
          shiny::tags$dd("Groups related samples. A matched tumor and normal must use the same patient ID."),
          shiny::tags$dt("Sex chromosomes"),
          shiny::tags$dd("Patient-level XX or XY used by CNV callers. If not provided, ASCAT or Control-FREEC is skipped and the remaining callers still run."),
          shiny::tags$dt("Sample ID"),
          shiny::tags$dd("Identifies one biological sample within a patient. It must be unique within that patient."),
          shiny::tags$dt("Role"),
          shiny::tags$dd("Germline, tumor, normal, or unknown. Filename-based role guesses must be reviewed."),
          shiny::tags$dt("Matched normal"),
          shiny::tags$dd("For a tumor sample, select an included normal sample belonging to the same patient.")
        )
      ),
      shiny::column(
        4,
        shiny::h5("Input files"),
        shiny::tags$dl(
          shiny::tags$dt("Input format"),
          shiny::tags$dd("Detected from the file extension and read-only: FASTQ, uBAM, BAM, CRAM, VCF, or BCF."),
          shiny::tags$dt("Processing state"),
          shiny::tags$dd("How far the input has already been processed. BAM recommendations come from lightweight header inspection and still require user confirmation."),
          shiny::tags$dt("FASTQ lane and read"),
          shiny::tags$dd("Every included sample lane must contain exactly one R1 and one R2. Correct detection only when the filename was interpreted incorrectly."),
          shiny::tags$dt("Index"),
          shiny::tags$dd("Detected companion index for BAM, CRAM, VCF, or BCF. BAM/CRAM submission is blocked when its index is missing or unreadable.")
        )
      ),
      shiny::column(
        4,
        shiny::h5("Analysis and storage"),
        shiny::tags$dl(
          shiny::tags$dt("Analysis mode"),
          shiny::tags$dd("Must agree with the included roles: all germline, all tumor, or matched tumor-normal."),
          shiny::tags$dt("Results root"),
          shiny::tags$dd("Permanent destination for final pipeline outputs and reports. It must be an absolute writable server path."),
          shiny::tags$dt("Work root"),
          shiny::tags$dd("Temporary space used by Nextflow for intermediate files. The default follows the Unix home directory; users with limited home quotas should select high-capacity storage."),
          shiny::tags$dt("Manifest ID"),
          shiny::tags$dd("Short name used for the submitted run, its results folder, and its cluster job.")
        )
      )
    )
  )
}

sarek_manifest_table_output <- function(output_id) {
  if (requireNamespace("DT", quietly = TRUE)) {
    DT::DTOutput(output_id)
  } else {
    shiny::tableOutput(output_id)
  }
}

sarek_manifest_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::tags$style(shiny::HTML(
      ".sarek-confirmation-table { width:100%; max-width:100%; overflow-x:auto; padding-bottom:8px; }
       .sarek-confirmation-table .dataTables_wrapper { overflow-x:visible; }
       .sarek-confirmation-table .dataTables_scroll { width:100%; overflow-x:auto; }
       .sarek-confirmation-table .dataTables_scrollHeadInner,
       .sarek-confirmation-table table.dataTable {
         width:max-content !important;
         min-width:2550px !important;
         table-layout:auto !important;
       }
       .sarek-confirmation-table table.dataTable thead th,
       .sarek-confirmation-table table.dataTable tbody td {
         min-width:110px !important;
         max-width:none !important;
         white-space:nowrap !important;
       }
       .sarek-confirmation-table table.dataTable thead th:nth-child(8),
       .sarek-confirmation-table table.dataTable tbody td:nth-child(8),
       .sarek-confirmation-table table.dataTable thead th:nth-child(9),
       .sarek-confirmation-table table.dataTable tbody td:nth-child(9),
       .sarek-confirmation-table table.dataTable thead th:nth-child(14),
       .sarek-confirmation-table table.dataTable tbody td:nth-child(14) {
         min-width:320px !important;
       }
       .sarek-status-banner {
         display:flex;
         align-items:flex-start;
         gap:10px;
         width:100%;
         margin:8px 0 16px 0;
         padding:13px 15px;
         border:1px solid #b8d4f2;
         border-left:6px solid #1769aa;
         border-radius:8px;
         background:#eef6ff;
         color:#17324d;
         font-size:14px;
         font-weight:650;
         line-height:1.45;
         white-space:normal;
         overflow-wrap:anywhere;
         box-shadow:0 2px 7px rgba(23,50,77,.08);
       }
       .sarek-status-banner::before {
         content:'i';
         flex:0 0 22px;
         width:22px;
         height:22px;
         border-radius:50%;
         background:#1769aa;
         color:#fff;
         text-align:center;
         font-weight:800;
         line-height:22px;
       }
       .sarek-bam-inspection {
         margin:12px 0 16px 0;
         padding:14px 16px;
         border:1px solid #c9d6e4;
         border-left:6px solid #1769aa;
         border-radius:8px;
         background:#f7fbff;
       }
       .sarek-bam-inspection-pass { border-left-color:#27864a; background:#f0faf3; }
       .sarek-bam-inspection-failed { border-left-color:#b3261e; background:#fff3f2; }
       .sarek-bam-inspection-review { border-left-color:#b26a00; background:#fff8e8; }
       .sarek-bam-inspection dl { display:grid; grid-template-columns:minmax(150px,220px) 1fr; gap:6px 14px; margin:10px 0; }
       .sarek-bam-inspection dt { font-weight:700; }
       .sarek-bam-inspection dd { margin:0; overflow-wrap:anywhere; }
       .sarek-status-banner-review {
         border-color:#e7c66a;
         border-left-color:#b7791f;
         background:#fff8df;
         color:#5c3b08;
       }
       .sarek-status-banner-review::before {
         content:'!';
         background:#b7791f;
       }
       .sarek-status-banner-error {
         border-color:#efb0b0;
         border-left-color:#b42318;
         background:#fff1f0;
         color:#7a1c16;
       }
       .sarek-status-banner-error::before {
         content:'!';
         background:#b42318;
       }
       .sarek-review-checklist {
         margin:10px 0 16px 0;
         padding:12px 14px;
         border:1px solid #d7e0ea;
         border-radius:8px;
         background:#f8fafc;
       }
       .sarek-review-item { margin:5px 0; line-height:1.4; }
       .sarek-review-item-pass { color:#276749; }
       .sarek-review-item-review { color:#8a5700; font-weight:650; }
       .sarek-sample-editor {
         margin:8px 0 16px 0;
         padding:14px;
         border:1px solid #d7e0ea;
         border-radius:8px;
         background:#fff;
       }
       .sarek-sample-editor h4 { margin-top:0; }
       .sarek-sample-review-table { width:100%; max-width:100%; overflow-x:auto; }
       .sarek-sample-review-table table.dataTable {
         width:max-content !important;
         min-width:1450px !important;
         table-layout:auto !important;
       }
       .sarek-sample-review-table table.dataTable th,
       .sarek-sample-review-table table.dataTable td {
         min-width:120px !important;
         max-width:none !important;
         white-space:nowrap !important;
       }
       .sarek-edit-instruction {
         margin:10px 0;
         padding:10px 12px;
         border-left:5px solid #1769aa;
         border-radius:6px;
         background:#eef6ff;
         color:#17324d;
         font-weight:650;
       }
       .sarek-field-guide {
         margin:12px 0 18px 0;
         padding:12px 14px;
         border:1px solid #cfd9e5;
         border-radius:8px;
         background:#fbfcfe;
       }
       .sarek-field-guide summary { cursor:pointer; font-weight:700; }
       .sarek-field-guide h5 { margin-top:16px; font-weight:700; }
       .sarek-field-guide dt { margin-top:8px; color:#17324d; }
       .sarek-field-guide dd { margin-left:0; color:#536273; line-height:1.4; }
       .sarek-validation-panel {
         display:flex;
         align-items:flex-start;
         gap:12px;
         margin:10px 0 16px 0;
         padding:15px 17px;
         border:1px solid #b8d4f2;
         border-left:7px solid #1769aa;
         border-radius:8px;
         background:#eef6ff;
         color:#17324d;
         box-shadow:0 2px 8px rgba(23,50,77,.08);
       }
       .sarek-validation-panel::before {
         content:'i';
         flex:0 0 24px;
         width:24px;
         height:24px;
         border-radius:50%;
         background:#1769aa;
         color:#fff;
         text-align:center;
         font-weight:800;
         line-height:24px;
       }
       .sarek-validation-panel h4 { margin:1px 0 7px 0; font-weight:750; }
       .sarek-validation-panel pre {
         margin:0;
         padding:0;
         border:0;
         background:transparent;
         color:inherit;
         font-family:inherit;
         font-size:14px;
         line-height:1.5;
         white-space:pre-wrap;
         overflow-wrap:anywhere;
       }
       .sarek-validation-panel-success {
         border-color:#9fd4b2;
         border-left-color:#2f855a;
         background:#edf9f1;
         color:#22543d;
       }
       .sarek-validation-panel-success::before { content:'✓'; background:#2f855a; }
       .sarek-validation-panel-warning {
         border-color:#e7c66a;
         border-left-color:#b7791f;
         background:#fff8df;
         color:#5c3b08;
       }
       .sarek-validation-panel-warning::before { content:'!'; background:#b7791f; }
       .sarek-validation-panel-error {
         border-color:#efb0b0;
         border-left-color:#b42318;
         background:#fff1f0;
         color:#7a1c16;
       }
       .sarek-validation-panel-error::before { content:'!'; background:#b42318; }
       .sarek-storage-note {
         margin:2px 0 12px 0;
         color:#536273;
         font-size:13px;
         line-height:1.4;
       }
       .sarek-run-review {
         margin:10px 0 16px 0;
         padding:16px;
         border:1px solid #d7e0ea;
         border-radius:10px;
         background:#fbfcfe;
       }
       .sarek-run-review-grid {
         display:grid;
         grid-template-columns:repeat(auto-fit,minmax(190px,1fr));
         gap:10px;
       }
       .sarek-run-review-card {
         min-width:0;
         padding:11px 12px;
         border:1px solid #dce4ec;
         border-radius:8px;
         background:#fff;
       }
       .sarek-run-review-card span {
         display:block;
         margin-bottom:4px;
         color:#647386;
         font-size:12px;
         font-weight:700;
         letter-spacing:.03em;
         text-transform:uppercase;
       }
       .sarek-run-review-card strong {
         display:block;
         color:#17324d;
         line-height:1.35;
         overflow-wrap:anywhere;
       }
       .sarek-submit-panel {
         margin:10px 0 16px 0;
         padding:15px 17px;
         border:1px solid #b8d4f2;
         border-radius:9px;
         background:#f5f9ff;
       }
       .sarek-submit-panel .checkbox { margin-top:0; }
       .sarek-submit-panel .btn-primary { min-width:190px; font-weight:700; }
       .sarek-advanced-downloads {
         margin-top:12px;
         color:#536273;
       }
       .sarek-file-details {
         clear:both;
         display:block;
         margin-top:18px;
         padding-top:4px;
       }
       .sarek-file-details > summary {
         font-weight:700;
         color:#17324d;
       }
       .sarek-sample-review-table .shiny-output-error,
       .sarek-file-details .shiny-output-error {
         display:block;
         position:static;
         margin:10px 0;
         padding:10px 12px;
         white-space:normal;
         overflow-wrap:anywhere;
       }
       .sarek-run-history {
         margin:10px 0 18px 0;
         padding:14px 16px;
         border:1px solid #c9d8e8;
         border-radius:10px;
         background:#f8fbff;
       }
       .sarek-run-history > summary {
         cursor:pointer;
         color:#17324d;
         font-size:16px;
         font-weight:750;
       }
       .sarek-run-history-actions,
       .sarek-path-actions {
         display:flex;
         flex-wrap:wrap;
         gap:8px;
         margin:6px 0 10px 0;
       }
       .sarek-progress-shell {
         height:24px;
         margin:10px 0 8px 0;
         overflow:hidden;
         border-radius:12px;
         background:#e5ebf2;
       }
       .sarek-progress-bar {
         min-width:8px;
         height:100%;
         color:#fff;
         background:#2d78c4;
         font-size:12px;
         font-weight:750;
         line-height:24px;
         text-align:center;
         transition:width .35s ease;
       }
       .sarek-progress-bar.running {
         background-image:linear-gradient(45deg,rgba(255,255,255,.18) 25%,transparent 25%,transparent 50%,rgba(255,255,255,.18) 50%,rgba(255,255,255,.18) 75%,transparent 75%,transparent);
         background-size:32px 32px;
       }
       .sarek-progress-bar.success { background:#218739; }
       .sarek-progress-bar.error { background:#b42318; }
       .sarek-progress-bar.unknown { background:#687586; }
       .sarek-live-activity {
         margin:12px 0;
         padding:14px 16px;
         border:1px solid #bfd5ea;
         border-left:6px solid #2d78c4;
         border-radius:9px;
         background:#f4f9ff;
       }
       .sarek-live-activity h4 { margin:0 0 9px 0; }
       .sarek-live-activity-empty { border-left-color:#8795a5; background:#f7f8fa; }
       .sarek-live-task { margin:5px 0 10px 0; overflow-wrap:anywhere; }
       .sarek-live-counts { display:flex; flex-wrap:wrap; gap:10px; margin:7px 0; }
       .sarek-live-counts span { padding:4px 9px; border-radius:12px; background:#e3edf7; }
       .sarek-live-estimate { margin:9px 0 4px 0; font-weight:650; }
       .sarek-live-log { margin-top:10px; }
       .sarek-live-log summary { cursor:pointer; font-weight:700; }
       .sarek-live-log pre {
         max-height:260px;
         margin:8px 0 0 0;
         padding:10px;
         overflow:auto;
         border:1px solid #d7e0ea;
         border-radius:6px;
         background:#fff;
         white-space:pre-wrap;
         overflow-wrap:anywhere;
       }
       .sarek-run-paths {
         margin-top:8px;
         color:#536273;
         font-size:12px;
         line-height:1.45;
         overflow-wrap:anywhere;
       }
       .sarek-section-tabs > .nav-tabs {
         margin:0 0 20px 0;
         border-bottom:2px solid #d7e0ea;
       }
       .sarek-section-tabs > .nav-tabs > li > a {
         padding:12px 20px;
         color:#42566d;
         font-size:15px;
         font-weight:700;
       }
       .sarek-section-tabs > .nav-tabs > li.active > a,
       .sarek-section-tabs > .nav-tabs > li.active > a:hover,
       .sarek-section-tabs > .nav-tabs > li.active > a:focus {
         border:1px solid #c8d6e5;
         border-bottom-color:#fff;
         color:#1769aa;
       }
       .sarek-results-table {
         width:100%;
         max-width:100%;
         margin-top:14px;
         overflow-x:auto;
       }
       .sarek-results-controls {
         position:relative;
         z-index:3;
         overflow:visible;
       }
       .sarek-results-table table.dataTable th,
       .sarek-results-table table.dataTable td {
         white-space:nowrap;
       }
       .sarek-variant-table .dataTables_scrollBody table.dataTable {
         min-width:1900px;
       }
       .sarek-variant-table table.dataTable th,
       .sarek-variant-table table.dataTable td {
         min-width:110px;
       }
       .sarek-explorer-view {
         position:relative;
         z-index:1;
         max-height:76vh;
         overflow:auto;
         margin-top:18px;
         padding:16px;
         border:1px solid #d7e0ea;
         background:#fff;
       }"
    )),
    shiny::div(
      class = "sarek-section-tabs",
      shiny::tabsetPanel(
        id = ns("sarek_sections"),
        selected = "prepare",
        type = "tabs",
        shiny::tabPanel(
          "Prepare job",
          value = "prepare",
          shiny::div(
            class = "progress-header-row",
            shiny::div(
              shiny::h3("Prepare a Sarek analysis"),
              shiny::tags$p(
                class = "muted",
                "Discover sequencing files, review every inferred field, and submit a confirmed analysis to the cluster."
              )
            )
          ),
          sarek_manifest_field_guide(),
    shiny::tags$div(
      class = "read-source-note",
      shiny::tags$strong("Step 1: discover files"),
      shiny::tags$p(
        "Enter one readable absolute server path per line. A path may point to a file or folder."
      ),
      shiny::tags$p(
        class = "muted small-note",
        "Supported inputs: paired FASTQ, uBAM, BAM, CRAM, VCF, and BCF. Discovered BAMs receive a lightweight header and index inspection; no reads or variants are scanned."
      )
    ),
    shiny::fluidRow(
      shiny::column(
        8,
        shiny::textAreaInput(
          ns("paths"),
          "Input files or folders",
          value = "",
          rows = 8,
          placeholder = "/absolute/path/sample_T_L001_R1.fastq.gz\n/absolute/path/sample_T_L001_R2.fastq.gz"
        ),
        shiny::div(
          class = "sarek-path-actions",
          shiny::actionButton(ns("browse_input_folder"), "Browse and add server folder")
        ),
        shiny::checkboxInput(ns("recursive"), "Search inside subfolders", value = FALSE),
        shiny::actionButton(ns("discover"), "Discover inputs", class = "btn-primary")
      ),
      shiny::column(
        4,
        shiny::uiOutput(ns("status"))
      )
    ),
    shiny::tags$hr(),
    shiny::tags$div(
      class = "read-source-note",
      shiny::tags$strong("Step 2: review and correct samples"),
      shiny::tags$p(
        "The searchable table uses one row per biological sample. Select a row once to edit it below; changes apply consistently to every associated file."
      ),
      shiny::tags$p(
        class = "muted small-note",
        "File format, paths, indexes, sizes, confidence, and warnings remain read-only. Lane/read detection can be corrected inside file details."
      )
    ),
    shiny::tableOutput(ns("summary")),
    shiny::uiOutput(ns("readiness")),
    shiny::div(
      class = "sarek-edit-instruction",
      "Editable: click any sample row once, then use the checkbox, text fields, and dropdowns in the selected-sample editor below."
    ),
    shiny::div(
      class = "sarek-sample-review-table",
      sarek_manifest_table_output(ns("sample_review_table"))
    ),
    shiny::uiOutput(ns("sample_editor")),
    shiny::uiOutput(ns("bam_inspection")),
    shiny::tags$details(
      class = "sarek-file-details",
      shiny::tags$summary("Show input and index paths for the selected sample"),
      shiny::br(),
      shiny::uiOutput(ns("file_editor")),
      shiny::div(
        class = "sarek-confirmation-table",
        sarek_manifest_table_output(ns("file_detail_table"))
      )
    ),
    shiny::tags$hr(),
    shiny::tags$div(
      class = "read-source-note",
      shiny::tags$strong("Step 3: describe the intended analysis"),
      shiny::tags$p(
        "These settings define the run that will be shown for final review before submission."
      )
    ),
    shiny::fluidRow(
      shiny::column(
        4,
        shiny::textInput(ns("manifest_id"), "Manifest ID", value = "sarek_analysis"),
        shiny::selectInput(
          ns("assay_type"),
          "Assay",
          choices = c("Whole-genome sequencing" = "WGS", "Whole-exome sequencing" = "WES",
                      "Targeted panel" = "targeted", "Annotation only" = "annotation_only"),
          selected = "WGS"
        ),
        shiny::selectInput(
          ns("analysis_mode"),
          "Analysis mode",
          choices = c("Germline" = "germline", "Tumor only" = "tumor_only",
                      "Matched tumor-normal" = "matched_tumor_normal",
                      "Annotation only" = "annotation_only"),
          selected = "germline",
          selectize = FALSE
        ),
        shiny::selectInput(
          ns("preset"),
          "Analysis preset",
          choices = c(
            "Core, mode-aware variant calling" = "core",
            "Small variants + VEP (report)" = "small_variants_vep"
          ),
          selected = "core",
          selectize = FALSE
        )
      ),
      shiny::column(
        4,
        shiny::textInput(
          ns("results_root"),
          "Results root — permanent outputs",
          value = "",
          placeholder = "/absolute/path/results/sarek/user"
        ),
        shiny::actionButton(ns("browse_results_root"), "Browse server folder"),
        shiny::actionButton(ns("save_results_root"), "Save as my Sarek default"),
        shiny::uiOutput(ns("default_results_root_status")),
        shiny::tags$p(
          class = "sarek-storage-note",
          "Final results, reports, and deliverables will be stored here."
        ),
        shiny::textInput(
          ns("work_root"),
          "Work root — temporary Nextflow files",
          value = "",
          placeholder = "/high-capacity/path/work/sarek/user"
        ),
        shiny::actionButton(ns("browse_work_root"), "Browse server folder"),
        shiny::tags$p(
          class = "sarek-storage-note",
          "Intermediate files can be much larger than the final results. The default follows your Unix home; use Browse to select high-capacity storage if your home quota is limited."
        ),
        shiny::tags$p(
          class = "muted small-note",
          "Both values must be absolute server paths. Manifest creation records them but does not create folders or start Sarek."
        )
      ),
      shiny::column(
        4,
        shiny::tags$details(
          shiny::tags$summary("Reference settings"),
          shiny::br(),
          shiny::textInput(ns("species"), "Species", value = "human"),
          shiny::textInput(ns("assembly"), "Assembly", value = "GRCh38"),
          shiny::textInput(ns("sarek_genome"), "nf-core/sarek genome key", value = "GATK.GRCh38")
        ),
        shiny::br(),
        shiny::actionButton(ns("confirm"), "Validate and confirm manifest", class = "btn-primary")
      )
    ),
    shiny::h4("Validation"),
    shiny::uiOutput(ns("validation")),
    shiny::h4("Nextflow input"),
    shiny::tags$p(
      class = "muted small-note",
      "A confirmed analysis is converted into the nf-core/sarek input and starting step. Submission remains locked until the review below is acknowledged."
    ),
    shiny::uiOutput(ns("nextflow_input")),
    shiny::h4("Confirmed run review"),
    shiny::uiOutput(ns("run_review")),
    shiny::tableOutput(ns("confirmed_samples")),
    shiny::h4("Submit run"),
    shiny::uiOutput(ns("submission_ui")),
          shiny::uiOutput(ns("submission_status"))
        ),
        shiny::tabPanel(
          "Run status",
          value = "status",
          shiny::div(
            class = "progress-header-row",
            shiny::div(
              shiny::h3("Sarek run status"),
              shiny::tags$p(
                class = "muted",
                "Track active Sarek runs and review previously submitted analyses."
              )
            )
          ),
          shiny::div(
            class = "sarek-run-history",
            shiny::tags$h4("Previous and active Sarek runs"),
            shiny::fluidRow(
              shiny::column(
                7,
                shiny::textInput(ns("history_root"), "Results location to view", value = ""),
                shiny::div(
                  class = "sarek-run-history-actions",
                  shiny::actionButton(ns("browse_history_root"), "Browse server folder"),
                  shiny::actionButton(ns("refresh_runs"), "Refresh runs")
                ),
                shiny::actionButton(ns("save_history_root"), "Save this as my Sarek default"),
                shiny::uiOutput(ns("history_default_results_root_status")),
                shiny::tags$p(
                  class = "muted small-note",
                  "Runs are discovered from their private CodeSpring submission records. Active status refreshes automatically."
                )
              ),
              shiny::column(5, shiny::uiOutput(ns("run_history_selector")))
            ),
            shiny::tableOutput(ns("run_history_table")),
            shiny::uiOutput(ns("selected_run_status"))
          )
        ),
        shiny::tabPanel(
          "Files",
          value = "files",
          shiny::div(
            class = "progress-header-row",
            shiny::div(
              shiny::h3("Sarek files"),
              shiny::tags$p(
                class = "muted",
                "Choose a submitted run and find its reports, variants, alignments, and other result files."
              )
            )
          ),
          shiny::fluidRow(
            shiny::column(7, shiny::uiOutput(ns("results_run_selector"))),
            shiny::column(
              5,
              shiny::actionButton(ns("refresh_results"), "Refresh results")
            )
          ),
          shiny::uiOutput(ns("selected_results_summary")),
          shiny::div(class = "sarek-results-controls", shiny::uiOutput(ns("files_filters"))),
          shiny::div(
            class = "sarek-results-table",
            sarek_manifest_table_output(ns("results_file_table"))
          )
        ),
        shiny::tabPanel(
          "Results Explorer",
          value = "explorer",
          shiny::div(class = "sarek-results-controls", shiny::fluidRow(
            shiny::column(6, shiny::uiOutput(ns("explorer_run_selector"))),
            shiny::column(4, shiny::uiOutput(ns("explorer_sample_selector"))),
            shiny::column(2, shiny::br(), shiny::actionButton(ns("refresh_explorer"), "Refresh"))
          )),
          shiny::tabsetPanel(
            id = ns("explorer_sections"),
            shiny::tabPanel(
              "QC",
              shiny::h4("Pipeline QC"),
              shiny::div(class = "sarek-results-controls", shiny::uiOutput(ns("qc_selector"))),
              shiny::downloadButton(ns("download_qc"), "Download selected QC report"),
              shiny::div(class = "sarek-explorer-view", shiny::uiOutput(ns("qc_preview")))
            ),
            shiny::tabPanel(
              "Variants",
              shiny::h4("Variant records"),
              shiny::uiOutput(ns("variant_scope")),
              shiny::downloadButton(ns("download_variants"), "Download inspected records"),
              shiny::div(class = "sarek-results-table sarek-variant-table", sarek_manifest_table_output(ns("variants_table")))
            ),
            shiny::tabPanel(
              "ClinVar",
              shiny::h4("ClinVar annotations"),
              shiny::tags$p(class = "muted", "Optional post-run annotation. Select a Sarek VCF and an indexed ClinVar VCF for the same genome build. Matching requires identical chromosome, position, REF, and ALT; original calls are preserved."),
              shiny::div(class = "sarek-results-controls", shiny::uiOutput(ns("clinvar_controls"))),
              shiny::actionButton(ns("start_clinvar"), "Run ClinVar annotation"),
                shiny::hr(),
                shiny::h4("Annotate several VCFs"),
                shiny::uiOutput(ns("clinvar_batch_selector")),
                shiny::actionButton(ns("start_clinvar_batch"), "Submit selected VCFs"),
              shiny::uiOutput(ns("clinvar_job_message")),
              shiny::h4("Previous ClinVar jobs"),
              shiny::tableOutput(ns("clinvar_jobs")),
              shiny::uiOutput(ns("clinvar_log_selector")),
              shiny::verbatimTextOutput(ns("clinvar_log")),
              shiny::hr(),
              shiny::h4("Browse annotated output"),
              shiny::div(
                class = "sarek-results-controls",
                shiny::uiOutput(ns("clinvar_output_selector")),
                shiny::selectInput(
                  ns("clinvar_record_mode"), "Records to show",
                  choices = c(
                    "All records in preview" = "all",
                    "With ClinVar classification only" = "classified"
                  ),
                  selected = "all", selectize = FALSE
                )
              ),
              shiny::uiOutput(ns("clinvar_scope")),
              shiny::downloadButton(ns("download_clinvar"), "Export preview CSV"),
              shiny::downloadButton(ns("download_clinvar_vcf"), "Download selected VCF"),
              shiny::div(class = "sarek-results-table sarek-variant-table", sarek_manifest_table_output(ns("clinvar_table"))),
                shiny::hr(),
                shiny::h4("ClinVar classifications across files"),
                shiny::tags$p(
                  class = "muted",
                  "Choose annotated VCFs, then scan for records with CLNSIG. ",
                  "The table shows the first 1,000 matches; the CSV export includes all matches."
                ),
                shiny::uiOutput(ns("clinvar_multi_selector")),
                shiny::actionButton(ns("clinvar_scan_multi"), "Scan selected VCFs"),
                shiny::uiOutput(ns("clinvar_multi_scope")),
                shiny::downloadButton(ns("download_clinvar_multi"), "Export all matches CSV"),
                shiny::div(
                  class = "sarek-results-table sarek-variant-table",
                  shiny::div(
                    style = "max-height: 520px; overflow-y: auto; overflow-x: auto;",
                    shiny::tableOutput(ns("clinvar_multi_table"))
                  )
                )
            ),
            shiny::tabPanel(
              "Plots",
              shiny::h4("Pipeline plots"),
              shiny::div(class = "sarek-results-controls", shiny::uiOutput(ns("plot_selector"))),
              shiny::downloadButton(ns("download_plot"), "Download selected plot"),
              shiny::div(class = "sarek-explorer-view", shiny::uiOutput(ns("plot_preview")))
            ),
            shiny::tabPanel(
              "Report",
              shiny::actionButton(ns("build_oncoplot"), "Build oncoplot"),
              shiny::downloadButton(ns("download_report"), "Download HTML report"),
              shiny::div(class = "sarek-explorer-view", shiny::uiOutput(ns("report_view")))
            )
          )
        )
      )
    )
  )
}

sarek_manifest_server <- function(
  id,
  default_results_root,
  default_work_root,
  results_root_preference_path = NULL,
  created_by = "unknown",
  allowed_input_roots = NULL,
  allowed_results_roots = NULL,
  allowed_work_roots = NULL,
  max_files = 5000L,
  submit_handler = NULL,
  browse_handler = NULL,
  bam_inspector = sarek_inspect_bam,
  samtools = "samtools",
  max_auto_bam_inspections = 20L,
  run_catalog_handler = sarek_submission_catalog,
  run_status_handler = sarek_run_status,
  clinvar_submit_handler = NULL,
  run_refresh_ms = 15000L
) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    confirmation_state <- shiny::reactiveVal(NULL)
    confirmed_manifest <- shiny::reactiveVal(NULL)
    status_state <- shiny::reactiveVal(
      "Enter one or more server paths, then select Discover inputs."
    )
    validation_state <- shiny::reactiveVal(
      "No manifest has been validated yet."
    )
    selected_sample_state <- shiny::reactiveVal("")
    validation_version <- shiny::reactiveVal(0L)
    reviewed_validation_version <- shiny::reactiveVal(-1L)
    submission_state <- shiny::reactiveVal(list(
      kind = "info",
      message = "Confirm and review a valid run before submission.",
      result = NULL
    ))
    submission_in_progress <- shiny::reactiveVal(FALSE)
    default_results_root_status <- shiny::reactiveVal("")
    clinvar_job_state <- shiny::reactiveVal("")
    run_history_refresh <- shiny::reactiveVal(Sys.time())
    run_poll <- shiny::reactiveTimer(max(5000L, as.integer(run_refresh_ms)), session = session)

    nextflow_input_state <- shiny::reactive({
      manifest <- confirmed_manifest()
      if (is.null(manifest)) {
        return(list(valid = FALSE, pending = TRUE, input = NULL, error = ""))
      }
      tryCatch(
        list(
          valid = TRUE,
          pending = FALSE,
          input = sarek_build_nextflow_input(manifest),
          error = ""
        ),
        error = function(error) {
          list(
            valid = FALSE,
            pending = FALSE,
            input = NULL,
            error = conditionMessage(error)
          )
        }
      )
    })

    submission_readiness <- shiny::reactive({
      manifest <- confirmed_manifest()
      nextflow_state <- nextflow_input_state()
      if (is.null(manifest) || !isTRUE(nextflow_state$valid)) {
        return(list(valid = FALSE, error = "The confirmed Sarek input is not ready."))
      }
      tryCatch(
        {
          paths <- sarek_submission_paths(manifest)
          params <- sarek_submission_params(manifest, nextflow_state$input, paths)
          list(valid = TRUE, error = "", params = params, paths = paths)
        },
        error = function(error) list(valid = FALSE, error = conditionMessage(error))
      )
    })

    shiny::observe({
      current_results_root <- if (nzchar(sarek_text(results_root_preference_path))) {
        sarek_results_read_default_root(
          results_root_preference_path, fallback = default_results_root
        )
      } else {
        default_results_root
      }
      shiny::updateTextInput(session, "results_root", value = current_results_root)
      shiny::updateTextInput(session, "work_root", value = default_work_root)
      shiny::updateTextInput(session, "history_root", value = current_results_root)
    })

    output$default_results_root_status <- shiny::renderUI({
      message <- default_results_root_status()
      if (nzchar(message)) shiny::tags$p(class = "muted small-note", message)
    })
    shiny::observeEvent(input$save_results_root, {
      tryCatch({
        shiny::req(nzchar(sarek_text(results_root_preference_path)))
        saved <- sarek_results_save_default_root(
          input$results_root, results_root_preference_path,
          allowed_roots = allowed_results_roots
        )
        shiny::updateTextInput(session, "history_root", value = saved)
        default_results_root_status(paste("Sarek results default saved:", saved))
        run_history_refresh(Sys.time())
      }, error = function(error) {
        default_results_root_status(paste("Could not save default:", conditionMessage(error)))
      })
    }, ignoreInit = TRUE)

    output$history_default_results_root_status <- shiny::renderUI({
      message <- default_results_root_status()
      if (nzchar(message)) shiny::tags$p(class = "muted small-note", message)
    })

    shiny::observeEvent(input$save_history_root, {
      tryCatch({
        preference_path <- sarek_text(results_root_preference_path)
        if (!nzchar(preference_path)) {
          stop("Saving a Sarek default is unavailable in this session.", call. = FALSE)
        }
        saved <- sarek_results_save_default_root(
          input$history_root, preference_path,
          allowed_roots = allowed_results_roots
        )
        shiny::updateTextInput(session, "history_root", value = saved)
        shiny::updateTextInput(session, "results_root", value = saved)
        default_results_root_status(paste("Sarek results default saved:", saved))
        run_history_refresh(Sys.time())
      }, error = function(error) {
        default_results_root_status(
          paste("Could not save default:", conditionMessage(error))
        )
      })
    }, ignoreInit = TRUE)

    request_browser <- function(target, current, input_type = "text", append = FALSE) {
      if (!is.function(browse_handler)) {
        status_state("ERROR: Server folder browsing is not configured for this app session.")
        return(invisible(FALSE))
      }
      tryCatch(
        {
          browse_handler(
            target = ns(target), mode = "dir", current = current,
            input_type = input_type, append = append
          )
          invisible(TRUE)
        },
        error = function(error) {
          status_state(paste0("ERROR: Could not open the server folder browser: ", conditionMessage(error)))
          invisible(FALSE)
        }
      )
    }

    shiny::observeEvent(input$browse_input_folder, {
      paths <- sarek_parse_path_input(input$paths)
      current <- if (length(paths)) paths[[length(paths)]] else default_results_root
      request_browser("paths", current, input_type = "textarea", append = TRUE)
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$browse_results_root, {
      request_browser("results_root", sarek_shiny_value(input$results_root, default_results_root))
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$browse_work_root, {
      request_browser("work_root", sarek_shiny_value(input$work_root, default_work_root))
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$browse_history_root, {
      request_browser("history_root", sarek_shiny_value(input$history_root, default_results_root))
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$refresh_runs, {
      run_history_refresh(Sys.time())
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$refresh_results, {
      run_history_refresh(Sys.time())
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$refresh_explorer, {
      run_history_refresh(Sys.time())
    }, ignoreInit = TRUE)

    run_catalog <- shiny::reactive({
      run_poll()
      run_history_refresh()
      root <- trimws(as.character(sarek_shiny_value(input$history_root, default_results_root)))
      if (!is.function(run_catalog_handler)) return(data.frame())
      tryCatch(run_catalog_handler(root), error = function(error) data.frame())
    })

    output$run_history_selector <- shiny::renderUI({
      catalog <- run_catalog()
      if (!NROW(catalog)) {
        return(shiny::div(class = "empty-box", "No submitted Sarek runs were found in this results location."))
      }
      labels <- vapply(seq_len(NROW(catalog)), function(index) {
        sarek_run_history_label(as.list(catalog[index, , drop = FALSE]))
      }, character(1))
      choices <- stats::setNames(as.character(catalog$run_id), labels)
      selected <- as.character(sarek_shiny_value(input$selected_run))
      if (!selected %in% catalog$run_id) selected <- catalog$run_id[[1]]
      shiny::selectInput(ns("selected_run"), "Run to view", choices = choices, selected = selected, selectize = FALSE)
    })

    output$run_history_table <- shiny::renderTable({
      catalog <- run_catalog()
      if (!NROW(catalog)) return(NULL)
      shown <- utils::head(catalog, 25L)
      data.frame(
        Run = shown$run_id,
        Submitted = gsub("T", " ", sub("Z$", " UTC", shown$submitted_at)),
        `Controller PID` = if ("controller_pid" %in% names(shown)) shown$controller_pid else shown$job_id,
        `Starting step` = shown$step,
        Tools = shown$tools,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = FALSE)

    selected_run <- shiny::reactive({
      catalog <- run_catalog()
      shiny::req(NROW(catalog))
      selected <- as.character(sarek_shiny_value(input$selected_run, catalog$run_id[[1]]))
      if (!selected %in% catalog$run_id) selected <- catalog$run_id[[1]]
      catalog[catalog$run_id == selected, , drop = FALSE][1, , drop = FALSE]
    })

    output$selected_run_status <- shiny::renderUI({
      run_poll()
      run <- selected_run()
      status <- if (is.function(run_status_handler)) {
        tryCatch(run_status_handler(as.list(run)), error = function(error) {
          list(state = "", elapsed = "", source = "error", error = conditionMessage(error))
        })
      } else {
        list(state = sarek_text(run$status, "submitted"), elapsed = "", source = "record")
      }
      progress <- sarek_run_progress(status$state)
      activity <- sarek_run_activity(as.list(run), state = status$state)
      details <- c(
        if ("controller_pid" %in% names(run) && nzchar(sarek_text(run$controller_pid))) {
          paste0("Nextflow controller PID ", sarek_text(run$controller_pid))
        } else if (nzchar(sarek_text(run$job_id))) {
          paste0("Slurm job ", sarek_text(run$job_id))
        } else {
          "No active controller ID recorded"
        },
        if (nzchar(sarek_text(status$elapsed))) paste0("elapsed ", sarek_text(status$elapsed)) else NULL,
        paste0("status source: ", sarek_text(status$source, "record"))
      )
      shiny::div(
        class = "sarek-run-review",
        shiny::tags$h4(paste0(sarek_text(run$run_id), " — ", progress$label)),
        shiny::tags$p(class = "muted small-note", paste(details, collapse = " · ")),
        shiny::div(
          class = "sarek-progress-shell",
          role = "progressbar",
          `aria-valuemin` = "0",
          `aria-valuemax` = "100",
          `aria-valuenow` = as.character(progress$percent),
          shiny::div(
            class = paste("sarek-progress-bar", progress$kind),
            style = paste0("width:", progress$percent, "%;"),
            progress$label
          )
        ),
        sarek_run_activity_panel(activity, status$state),
        if (nzchar(sarek_text(status$error))) {
          shiny::tags$p(class = "sarek-validation-panel sarek-validation-panel-error", sarek_text(status$error))
        },
        shiny::div(
          class = "sarek-run-paths",
          shiny::tags$strong("Results: "), sarek_text(run$output_dir),
          shiny::tags$br(),
          shiny::tags$strong("Run record: "), sarek_text(run$run_dir),
          if (nzchar(sarek_text(run$work_dir))) shiny::tagList(
            shiny::tags$br(), shiny::tags$strong("Temporary work: "), sarek_text(run$work_dir)
          )
        ),
        shiny::actionButton(ns("open_selected_results"), "View this run's results")
      )
    })

    output$results_run_selector <- shiny::renderUI({
      catalog <- run_catalog()
      if (!NROW(catalog)) {
        return(shiny::div(class = "empty-box", "No submitted Sarek runs were found in this results location."))
      }
      labels <- vapply(seq_len(NROW(catalog)), function(index) {
        sarek_run_history_label(as.list(catalog[index, , drop = FALSE]))
      }, character(1))
      choices <- stats::setNames(as.character(catalog$run_id), labels)
      selected <- as.character(sarek_shiny_value(
        input$results_run,
        sarek_shiny_value(input$selected_run, catalog$run_id[[1]])
      ))
      if (!selected %in% catalog$run_id) selected <- catalog$run_id[[1]]
      shiny::selectInput(ns("results_run"), "Run results to view", choices = choices, selected = selected, selectize = FALSE)
    })

    selected_results_run <- shiny::reactive({
      catalog <- run_catalog()
      shiny::req(NROW(catalog))
      selected <- as.character(sarek_shiny_value(
        input$results_run,
        sarek_shiny_value(input$selected_run, catalog$run_id[[1]])
      ))
      if (!selected %in% catalog$run_id) selected <- catalog$run_id[[1]]
      catalog[catalog$run_id == selected, , drop = FALSE][1, , drop = FALSE]
    })

    results_file_catalog <- shiny::reactive({
      run_history_refresh()
      run <- selected_results_run()
      sarek_result_file_catalog(sarek_text(run$output_dir), max_files = 5000L)
    })

    files_artifacts <- shiny::reactive({
      sarek_results_artifact_catalog(as.list(selected_results_run()))
    })
    output$files_filters <- shiny::renderUI({
      catalog <- files_artifacts()
      samples <- c("All", sarek_explorer_samples(catalog))
      categories <- c("All", sort(unique(as.character(catalog$category))))
      sample <- shiny::isolate(sarek_text(input$files_sample, "All"))
      category <- shiny::isolate(sarek_text(input$files_category, "All"))
      shiny::fluidRow(
        shiny::column(6, shiny::selectInput(
          ns("files_sample"), "Sample", choices = samples,
          selected = if (sample %in% samples) sample else "All", selectize = FALSE
        )),
        shiny::column(6, shiny::selectInput(
          ns("files_category"), "Output type", choices = categories,
          selected = if (category %in% categories) category else "All", selectize = FALSE
        ))
      )
    })

    output$selected_results_summary <- shiny::renderUI({
      run <- selected_results_run()
      status <- if (is.function(run_status_handler)) {
        tryCatch(run_status_handler(as.list(run)), error = function(error) {
          list(state = "", elapsed = "", source = "error", error = conditionMessage(error))
        })
      } else {
        list(state = sarek_text(run$status, "submitted"), elapsed = "", source = "record")
      }
      progress <- sarek_run_progress(status$state)
      files <- results_file_catalog()
      output_dir <- sarek_text(run$output_dir)
      directory_ready <- nzchar(output_dir) && dir.exists(output_dir)
      truncated <- isTRUE(attr(files, "truncated"))
      panel_kind <- if (!directory_ready || !NROW(files)) "warning" else "success"
      shiny::div(
        class = paste("sarek-validation-panel", paste0("sarek-validation-panel-", panel_kind)),
        shiny::div(
          shiny::h4(paste0(sarek_text(run$run_id), " — ", progress$label)),
          if (!directory_ready) {
            shiny::tags$p("The results directory is not available yet. It may not have been created by the active run.")
          } else if (!NROW(files)) {
            shiny::tags$p("The results directory exists but does not contain readable result files yet.")
          } else {
            shiny::tags$p(
              paste0(
                NROW(files), " readable result file", if (NROW(files) == 1L) "" else "s", " found",
                if (truncated) " (display limited to the first 5,000 files)" else "", "."
              )
            )
          },
          shiny::tags$p(class = "sarek-run-paths", shiny::tags$strong("Results location: "), output_dir)
        )
      )
    })

    results_table_data <- shiny::reactive({
      files <- files_artifacts()
      if (nzchar(sarek_text(input$files_sample)) && input$files_sample != "All") {
        files <- sarek_results_filter_sample(files, input$files_sample)
      }
      if (nzchar(sarek_text(input$files_category)) && input$files_category != "All") {
        files <- files[files$category == input$files_category, , drop = FALSE]
      }
      shown <- files[, c("file", "folder", "sample", "category", "size", "modified"), drop = FALSE]
      names(shown) <- c("File", "Folder", "Sample", "Output type", "Size", "Modified")
      shown
    })

    if (requireNamespace("DT", quietly = TRUE)) {
      output$results_file_table <- DT::renderDT({
        shown <- results_table_data()
        shiny::req(NROW(shown))
        DT::datatable(
          shown,
          rownames = FALSE,
          filter = "top",
          selection = "none",
          options = list(
            pageLength = 50,
            lengthMenu = c(25, 50, 100, 250),
            deferRender = TRUE,
            autoWidth = TRUE,
            scrollX = TRUE
          ),
          class = "stripe hover compact"
        )
      })
    } else {
      output$results_file_table <- shiny::renderTable({
        shown <- results_table_data()
        if (!NROW(shown)) return(NULL)
        utils::head(shown, 100L)
      }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = FALSE)
    }

    output$explorer_run_selector <- shiny::renderUI({
      catalog <- run_catalog()
      if (!NROW(catalog)) return(shiny::div(class = "empty-box", "No submitted Sarek runs found."))
      labels <- vapply(seq_len(NROW(catalog)), function(i)
        sarek_run_history_label(as.list(catalog[i, , drop = FALSE])), character(1))
      choices <- stats::setNames(as.character(catalog$run_id), labels)
      selected <- shiny::isolate(sarek_text(
        input$explorer_run, sarek_text(input$selected_run, catalog$run_id[[1]])
      ))
      if (!selected %in% catalog$run_id) selected <- catalog$run_id[[1]]
      shiny::selectInput(ns("explorer_run"), "Run", choices, selected = selected, selectize = FALSE)
    })
    explorer_run <- shiny::reactive({
      run_history_refresh()
      input$history_root
      selected <- sarek_text(input$explorer_run)
      catalog <- shiny::isolate(run_catalog())
      shiny::req(NROW(catalog))
      if (!selected %in% catalog$run_id) selected <- catalog$run_id[[1]]
      as.list(catalog[catalog$run_id == selected, , drop = FALSE][1, , drop = FALSE])
    })
    explorer_catalog <- shiny::reactive({
      run <- explorer_run()
      sarek_results_artifact_catalog(run)
    })
    output$explorer_sample_selector <- shiny::renderUI({
      samples <- c("All", sarek_explorer_samples(explorer_catalog()))
      selected <- shiny::isolate(sarek_text(input$explorer_sample, "All"))
      shiny::selectInput(ns("explorer_sample"), "Sample",
                         samples, selected = if (selected %in% samples) selected else "All",
                         selectize = FALSE)
    })
    explorer_sample <- shiny::reactive({
      selected <- sarek_text(input$explorer_sample, "All")
      if (!selected %in% c("All", sarek_explorer_samples(explorer_catalog()))) "All" else selected
    })
    scoped_catalog <- shiny::reactive({
      catalog <- explorer_catalog()
      if (explorer_sample() != "All") catalog <- sarek_results_filter_sample(catalog, explorer_sample())
      catalog
    })
    qc_catalog <- shiny::reactive({
      sarek_explorer_qc(explorer_catalog(), explorer_sample())
    })
    plot_catalog <- shiny::reactive(sarek_explorer_plots(scoped_catalog()))
    artifact_selector <- function(catalog, id, label) {
      if (!NROW(catalog)) return(shiny::div(class = "empty-box", "No matching outputs were found for this run or sample."))
      labels <- paste0(basename(as.character(catalog$file)),
                       ifelse(!is.na(catalog$sample) & nzchar(catalog$sample),
                              paste0(" · ", catalog$sample), ""))
      labels <- make.unique(labels, sep = " · ")
      shiny::selectInput(ns(id), label,
                         stats::setNames(as.character(seq_len(NROW(catalog))), labels),
                         selected = "1", selectize = FALSE)
    }
    selected_artifact <- function(catalog, value) {
      index <- suppressWarnings(as.integer(value))
      if (!length(index) || is.na(index) || index < 1L || index > NROW(catalog)) return(NULL)
      catalog[index, , drop = FALSE]
    }
    output$qc_selector <- shiny::renderUI(artifact_selector(qc_catalog(), "qc_artifact", "QC report"))
    output$plot_selector <- shiny::renderUI(artifact_selector(plot_catalog(), "plot_artifact", "Plot"))
    output$qc_preview <- shiny::renderUI({
      row <- selected_artifact(qc_catalog(), input$qc_artifact)
      if (is.null(row)) return(shiny::div(class = "empty-box", "Choose an available QC report."))
      sarek_explorer_preview(row$path[[1]], explorer_run()$output_dir)
    })
    output$plot_preview <- shiny::renderUI({
      row <- selected_artifact(plot_catalog(), input$plot_artifact)
      if (is.null(row)) return(shiny::div(class = "empty-box", "Choose an available plot."))
      sarek_explorer_preview(row$path[[1]], explorer_run()$output_dir)
    })
    download_artifact <- function(catalog, value, target) {
      row <- selected_artifact(catalog, value)
      shiny::req(!is.null(row))
      path <- sarek_explorer_file(row$path[[1]], explorer_run()$output_dir)
      shiny::req(nzchar(path))
      shiny::req(isTRUE(file.copy(path, target, overwrite = TRUE)))
    }
    output$download_qc <- shiny::downloadHandler(
      filename = function() {
        row <- selected_artifact(qc_catalog(), input$qc_artifact)
        if (is.null(row)) "sarek-qc.html" else basename(row$file[[1]])
      },
      content = function(file) download_artifact(qc_catalog(), input$qc_artifact, file)
    )
    output$download_plot <- shiny::downloadHandler(
      filename = function() {
        row <- selected_artifact(plot_catalog(), input$plot_artifact)
        if (is.null(row)) "sarek-plot.png" else basename(row$file[[1]])
      },
      content = function(file) download_artifact(plot_catalog(), input$plot_artifact, file)
    )
    explorer_variants <- shiny::reactive({
      sarek_explorer_variants(
        sarek_clinvar_sources(explorer_catalog()),
        explorer_sample(),
        explorer_run()$output_dir
      )
    })
    report_catalog <- shiny::reactive({
      catalog <- explorer_catalog()
      annotated <- sarek_clinvar_outputs(catalog)
      rbind(annotated, catalog[!catalog$path %in% annotated$path, , drop = FALSE])
    })
    report_variants <- shiny::reactive({
      sarek_explorer_variants(report_catalog(), "All", explorer_run()$output_dir)
    })
    output$variant_scope <- shiny::renderUI({
      rows <- explorer_variants()
      shiny::tags$p(class = "muted",
        NROW(rows), " records inspected in ", attr(rows, "files_scanned"),
        " of ", attr(rows, "files_available"), " VCF files.",
        if (isTRUE(attr(rows, "limited"))) " This is a bounded preview of the run's VCFs." else "")
    })
    output$clinvar_controls <- shiny::renderUI({
      run <- explorer_run()
      manifest <- sarek_results_read_manifest(run)
      assembly <- if (is.null(manifest)) "Unknown" else sarek_text(manifest$reference$assembly, "Unknown")
      catalog <- sarek_clinvar_sources(explorer_catalog())
      if (explorer_sample() != "All") catalog <- sarek_results_filter_sample(catalog, explorer_sample())
      labels <- make.unique(paste0(basename(as.character(catalog$file)), " · ", catalog$folder), sep = " · ")
      current <- shiny::isolate(sarek_text(input$clinvar_reference))
      default_ref <- Sys.getenv(paste0("CSL_SAREK_CLINVAR_VCF_", assembly), unset = "")
      shiny::tagList(
        shiny::tags$p(class = "muted", "Run assembly: ", assembly,
                      ". The ClinVar VCF is read from the server; download and index the chosen release separately."),
        shiny::selectInput(ns("clinvar_source"), "VCF to annotate",
          choices = stats::setNames(as.character(seq_len(NROW(catalog))), labels), selectize = FALSE),
        shiny::textInput(ns("clinvar_reference"), "Indexed ClinVar VCF path on server (.vcf.gz)",
          value = if (nzchar(current)) current else default_ref),
        shiny::selectInput(ns("clinvar_assembly"), "Assembly of ClinVar VCF",
          choices = c("Confirm reference assembly" = "", "GRCh37", "GRCh38"), selected = "", selectize = FALSE)
      )
    })
    output$clinvar_batch_selector <- shiny::renderUI({
      catalog <- sarek_clinvar_sources(explorer_catalog())
      if (explorer_sample() != "All") {
        catalog <- sarek_results_filter_sample(catalog, explorer_sample())
      }
      if (!NROW(catalog)) {
        return(shiny::tags$p("No original VCFs are available for this sample."))
      }
      labels <- make.unique(
        paste0(basename(as.character(catalog$file)), " · ", catalog$folder),
        sep = " · "
      )
      paths <- as.character(catalog$path)
      selected <- shiny::isolate(input$clinvar_batch_files)
      selected <- intersect(as.character(selected), paths)
      shiny::tagList(
        shiny::selectizeInput(
          ns("clinvar_batch_files"), "Original VCFs to annotate",
          choices = stats::setNames(paths, labels),
          selected = selected, multiple = TRUE, width = "100%",
          options = list(dropdownParent = "body")
        ),
        shiny::actionButton(ns("clinvar_batch_select_all"), "Select all VCFs")
      )
    })

    shiny::observeEvent(input$clinvar_batch_select_all, {
      catalog <- sarek_clinvar_sources(explorer_catalog())
      if (explorer_sample() != "All") {
        catalog <- sarek_results_filter_sample(catalog, explorer_sample())
      }
      shiny::updateSelectizeInput(
        session, "clinvar_batch_files",
        selected = as.character(catalog$path)
      )
    }, ignoreInit = TRUE)

    output$clinvar_job_message <- shiny::renderUI({
      message <- clinvar_job_state()
      if (nzchar(message)) shiny::tags$p(class = "muted", message)
    })
    output$clinvar_jobs <- shiny::renderTable({
      run_poll()
      run_history_refresh()

      jobs <- sarek_clinvar_jobs(explorer_run())
      if (!NROW(jobs)) return(NULL)

      shown <- jobs[, c(
        "Source",
        "Output",
        "Job",
        "Status",
        "ExitCode",
        "ErrorLog"
      ), drop = FALSE]

      names(shown) <- c(
        "Source",
        "Output",
        "Job",
        "Status",
        "Exit code",
        "Error log"
      )

      shown
    }, striped = TRUE, bordered = TRUE, rownames = FALSE)

    output$clinvar_log_selector <- shiny::renderUI({
      run_poll()
      run_history_refresh()

      jobs <- sarek_clinvar_jobs(explorer_run())
      if (!NROW(jobs)) return(NULL)

      labels <- paste0(
        jobs$Source,
        " · ",
        jobs$Status,
        " · job ",
        jobs$Job
      )

      choices <- stats::setNames(
        as.character(seq_len(NROW(jobs))),
        labels
      )

      selected <- shiny::isolate(
        sarek_text(input$clinvar_log_job, "1")
      )

      if (!selected %in% unname(choices)) {
        selected <- unname(choices)[[1]]
      }

      shiny::selectInput(
        ns("clinvar_log_job"),
        "ClinVar job log",
        choices = choices,
        selected = selected,
        selectize = FALSE
      )
    })

    output$clinvar_log <- shiny::renderText({
      run_poll()
      run_history_refresh()

      jobs <- sarek_clinvar_jobs(explorer_run())

      if (!NROW(jobs)) {
        return("No ClinVar jobs have been submitted.")
      }

      index <- suppressWarnings(
        as.integer(input$clinvar_log_job)
      )

      if (!length(index) || is.na(index) ||
          index < 1L || index > NROW(jobs)) {
        return("Select a ClinVar job to view its log.")
      }

      row <- jobs[index, , drop = FALSE]

      heading <- paste(
        "Job", row$Job[[1]],
        "·", row$Status[[1]],
        if (nzchar(row$ExitCode[[1]])) {
          paste("· exit", row$ExitCode[[1]])
        } else {
          ""
        }
      )

      if (!nzchar(row$LogText[[1]])) {
        return(paste(heading, "\n\nNo log output is available yet."))
      }

      paste(heading, "\n\n", row$LogText[[1]])
    })

    shiny::observeEvent(input$start_clinvar_batch, {
      tryCatch({
        shiny::req(is.function(clinvar_submit_handler))
        run <- explorer_run()
        status <- if (is.function(run_status_handler)) {
          run_status_handler(run)
        } else {
          list(state = run$status)
        }
        if (!identical(sarek_normalize_slurm_state(status$state), "COMPLETED")) {
          stop("Wait for the Sarek run to complete before starting ClinVar annotation.", call. = FALSE)
        }
        manifest <- sarek_results_read_manifest(run)
        if (is.null(manifest)) {
          stop("The recorded Sarek manifest is unavailable.", call. = FALSE)
        }
        catalog <- sarek_clinvar_sources(explorer_catalog())
        if (explorer_sample() != "All") {
          catalog <- sarek_results_filter_sample(catalog, explorer_sample())
        }
        selected <- unique(as.character(input$clinvar_batch_files))
        selected <- selected[!is.na(selected) & nzchar(selected)]
        if (!length(selected)) {
          stop("Select at least one original VCF.", call. = FALSE)
        }
        if (!all(selected %in% as.character(catalog$path))) {
          stop("The VCF selection changed. Select the files again.", call. = FALSE)
        }
        submitted <- character()
        skipped <- character()
        failed <- character()
        previous <- sarek_clinvar_jobs(run)
        output_dir <- normalizePath(run$output_dir, winslash = "/", mustWork = TRUE)
        for (path in selected) {
          row <- catalog[which(as.character(catalog$path) == path)[1L], , drop = FALSE]
          label <- basename(row$file[[1]])
          relative <- substring(dirname(path), nchar(output_dir) + 2L)
          target <- file.path(
            output_dir, "annotation", "clinvar", relative,
            paste0(
              sub("\\.vcf(\\.gz)?$", "", basename(path), ignore.case = TRUE),
              ".clinvar.vcf.gz"
            )
          )
          completed <- NROW(previous) > 0L && any(
            !is.na(previous$Path) & previous$Path == target &
              !is.na(previous$Status) & previous$Status == "COMPLETED"
          )
          if (file.exists(target) || completed) {
            reason <- if (file.exists(target)) {
              "output already exists"
            } else {
              "completed job recorded"
            }
            skipped <- c(skipped, paste0(label, " (", reason, ")"))
            next
          }
          tryCatch({
            job <- clinvar_submit_handler(
              run, manifest, path,
              input$clinvar_reference, input$clinvar_assembly
            )
            submitted <- c(submitted, paste0(label, " (job ", job$job_id, ")"))
          }, error = function(error) {
            failed <- c(failed, paste0(label, ": ", conditionMessage(error)))
          })
        }
        clinvar_job_state(paste(
          if (length(submitted)) {
            paste("Submitted", length(submitted), "job(s):", paste(submitted, collapse = "; "))
          } else {
            "No jobs submitted."
          },
          if (length(skipped)) {
            paste("Skipped:", paste(skipped, collapse = "; "))
          } else {
            ""
          },
          if (length(failed)) {
            paste("Could not submit:", paste(failed, collapse = "; "))
          } else {
            ""
          }
        ))
        run_history_refresh(Sys.time())
      }, error = function(error) clinvar_job_state(conditionMessage(error)))
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$start_clinvar, {
      tryCatch({
        shiny::req(is.function(clinvar_submit_handler))
        run <- explorer_run()
        status <- if (is.function(run_status_handler)) run_status_handler(run) else list(state = run$status)
        if (!identical(sarek_normalize_slurm_state(status$state), "COMPLETED")) {
          stop("Wait for the Sarek run to complete before starting ClinVar annotation.", call. = FALSE)
        }
        catalog <- sarek_clinvar_sources(explorer_catalog())
        if (explorer_sample() != "All") catalog <- sarek_results_filter_sample(catalog, explorer_sample())
        row <- selected_artifact(catalog, input$clinvar_source)
        if (is.null(row)) stop("Select a source VCF first.", call. = FALSE)
        manifest <- sarek_results_read_manifest(run)
        if (is.null(manifest)) stop("The recorded Sarek manifest is unavailable.", call. = FALSE)
        submitted <- clinvar_submit_handler(run, manifest, row$path[[1]],
                                             input$clinvar_reference, input$clinvar_assembly)
        clinvar_job_state(paste("ClinVar job", submitted$job_id, "submitted for", basename(row$file[[1]]),
                                 ". Refresh after it completes to view annotated records."))
        run_history_refresh(Sys.time())
      }, error = function(error) clinvar_job_state(conditionMessage(error)))
    }, ignoreInit = TRUE)
    clinvar_view_catalog <- shiny::reactive({
      catalog <- sarek_clinvar_outputs(explorer_catalog())
      if (explorer_sample() != "All") {
        catalog <- sarek_results_filter_sample(catalog, explorer_sample())
      }
      if (!NROW(catalog)) return(catalog)
      catalog <- catalog[!duplicated(catalog$path), , drop = FALSE]
      catalog[order(catalog$folder, catalog$file), , drop = FALSE]
    })

    output$clinvar_output_selector <- shiny::renderUI({
      catalog <- clinvar_view_catalog()
      if (!NROW(catalog)) {
        return(shiny::tags$p(
          "No annotated VCF outputs were found for this run and sample. ",
          "Check the annotation jobs above and refresh after completion."
        ))
      }

      callers <- vapply(seq_len(NROW(catalog)), function(i) {
        sarek_explorer_caller(catalog$folder[[i]], catalog$file[[i]])
      }, character(1))
      labels <- paste0(
        callers, " | ", catalog$folder, "/", basename(catalog$file)
      )
      paths <- as.character(catalog$path)
      selected <- shiny::isolate(sarek_text(input$clinvar_output_file))
      if (!selected %in% paths) selected <- paths[[1]]

      shiny::selectInput(
        ns("clinvar_output_file"), "Annotated VCF to view",
        choices = stats::setNames(paths, labels),
        selected = selected, selectize = FALSE, width = "100%"
      )
    })

    output$clinvar_multi_selector <- shiny::renderUI({
      catalog <- clinvar_view_catalog()
      if (!NROW(catalog)) {
        return(shiny::tags$p("No annotated VCF outputs are available."))
      }
      callers <- vapply(seq_len(NROW(catalog)), function(i) {
        sarek_explorer_caller(catalog$folder[[i]], catalog$file[[i]])
      }, character(1))
      paths <- as.character(catalog$path)
      labels <- make.unique(
        paste0(callers, " | ", catalog$folder, "/", basename(catalog$file)),
        sep = " · "
      )
      previous <- shiny::isolate(input$clinvar_multi_files)
      selected <- if (is.null(previous)) paths else intersect(as.character(previous), paths)
      shiny::selectizeInput(
        ns("clinvar_multi_files"), "Annotated VCFs",
        choices = stats::setNames(paths, labels),
        selected = selected, multiple = TRUE, width = "100%"
      )
    })

    clinvar_multi_catalog <- shiny::reactive({
      catalog <- clinvar_view_catalog()
      selected <- unique(as.character(input$clinvar_multi_files))
      selected <- selected[!is.na(selected) & nzchar(selected)]
      catalog[catalog$path %in% selected, , drop = FALSE]
    })

    clinvar_multi_data <- shiny::eventReactive(input$clinvar_scan_multi, {
      catalog <- clinvar_multi_catalog()
      if (!NROW(catalog)) {
        return(list(
          rows = sarek_explorer_empty_variants(), matches = 0L,
          files = 0L, error = "Select at least one annotated VCF."
        ))
      }
      tryCatch(
        sarek_explorer_clinsig_files(
          catalog, explorer_run()$output_dir, preview_limit = 1000L
        ),
        error = function(error) list(
          rows = sarek_explorer_empty_variants(), matches = 0L,
          files = 0L, error = conditionMessage(error)
        )
      )
    }, ignoreInit = TRUE)

    output$clinvar_multi_scope <- shiny::renderUI({
      if (is.null(input$clinvar_scan_multi) || input$clinvar_scan_multi < 1L) {
        return(shiny::tags$p(
          class = "muted",
          "Choose annotated VCFs above, then click Scan selected VCFs."
        ))
      }
      data <- clinvar_multi_data()
      if (nzchar(data$error)) {
        return(shiny::tags$p(class = "muted", data$error))
      }
      shiny::tags$p(
        class = "muted",
        data$matches, " records with CLNSIG found across ",
        data$files, " annotated VCF(s). Showing ",
        NROW(data$rows), " in the table. CSV export includes all matches."
      )
    })

    output$clinvar_multi_table <- shiny::renderTable({
      data <- clinvar_multi_data()
      if (nzchar(data$error) || !NROW(data$rows)) return(NULL)
      data$rows
    }, striped = TRUE, bordered = TRUE, rownames = FALSE)

    output$download_clinvar_multi <- shiny::downloadHandler(
      filename = function() "sarek_clinvar_all_matches.csv",
      content = function(file) {
        catalog <- clinvar_multi_catalog()
        shiny::req(NROW(catalog) > 0L)
        sarek_explorer_clinsig_files(
          catalog, explorer_run()$output_dir,
          preview_limit = 0L, csv_path = file
        )
      }
    )

    clinvar_view_file <- shiny::reactive({
      catalog <- clinvar_view_catalog()
      selected <- sarek_text(input$clinvar_output_file)
      catalog[which(as.character(catalog$path) == selected), , drop = FALSE]
    })

    clinvar_view_data <- shiny::reactive({
      file <- clinvar_view_file()
      result <- list(
        rows = sarek_explorer_empty_variants(),
        error = "", limited = FALSE
      )
      if (NROW(file) != 1L) return(result)

      tryCatch({
        rows <- sarek_explorer_variants(
          file,
          sample = explorer_sample(),
          output_dir = explorer_run()$output_dir,
          max_files = 1L,
          max_records = 1500L,
          strict = TRUE
        )
        list(
          rows = rows, error = "",
          limited = isTRUE(attr(rows, "limited"))
        )
      }, error = function(error) {
        result$error <- conditionMessage(error)
        result
      })
    })

    clinvar_has_classification <- function(rows) {
      values <- trimws(as.character(rows$ClinVar))
      !is.na(values) & nzchar(values) & values != "."
    }

    clinvar_variants <- shiny::reactive({
      rows <- clinvar_view_data()$rows
      if (identical(input$clinvar_record_mode, "classified")) {
        rows <- rows[clinvar_has_classification(rows), , drop = FALSE]
      }
      rows
    })

    clinvar_view_job <- shiny::reactive({
      run_poll()
      run_history_refresh()
      file <- clinvar_view_file()
      if (NROW(file) != 1L) return(NULL)
      jobs <- sarek_clinvar_jobs(explorer_run())
      jobs <- jobs[which(jobs$Path == file$path[[1]]), , drop = FALSE]
      if (!NROW(jobs)) return(NULL)
      utils::tail(jobs, 1L)
    })

    output$clinvar_scope <- shiny::renderUI({
      file <- clinvar_view_file()
      if (NROW(file) != 1L) {
        return(shiny::tags$p(class = "muted", "Select an annotated VCF above."))
      }

      preview <- clinvar_view_data()
      rows <- preview$rows
      job <- clinvar_view_job()
      n_classified <- sum(clinvar_has_classification(rows))

      message <- if (nzchar(preview$error)) {
        paste("Could not read the selected VCF:", preview$error)
      } else if (!NROW(rows)) {
        "The selected VCF has a header but contains no variant records."
      } else if (!n_classified && preview$limited) {
        paste(
          "No ClinVar classifications were found in these first 1,500 records.",
          "Later records have not been inspected."
        )
      } else if (!n_classified) {
        "No ClinVar classifications were found in this file."
      } else {
        ""
      }

      shiny::tagList(
        shiny::tags$p(
          shiny::tags$strong("Viewing one file: "),
          basename(file$file[[1]])
        ),
        shiny::tags$div(
          style = "overflow-wrap:anywhere; margin-bottom:10px;",
          shiny::tags$code(file$path[[1]])
        ),
        if (!is.null(job)) shiny::tags$p(
          "Latest recorded job: ", job$Job[[1]],
          " | ", job$Status[[1]],
          " | Exit code: ", job$ExitCode[[1]]
        ),
        if (!nzchar(preview$error)) shiny::tags$p(
          NROW(rows), " records read; ",
          n_classified, " with a ClinVar classification; ",
          NROW(clinvar_variants()), " in the current preview mode. ",
          if (preview$limited) {
            "Preview capped at the first 1,500 records; these are not full-file counts."
          } else {
            "All records in this file were read."
          }
        ),
        if (nzchar(message)) shiny::tags$p(
          style = if (nzchar(preview$error)) "color:#a12622;" else "",
          message
        ),
        shiny::tags$p(
          class = "muted",
          "Table searches apply only to loaded preview records. ",
          "CSV export uses the selected preview mode before table searches. ",
          "The VCF download contains the complete selected file."
        )
      )
    })

    explorer_table <- function(id, values) {
      if (requireNamespace("DT", quietly = TRUE)) {
        output[[id]] <- DT::renderDT({
          DT::datatable(values(), rownames = FALSE, filter = "top",
                        options = list(pageLength = 25, scrollX = TRUE, deferRender = TRUE,
                          autoWidth = FALSE, columnDefs = list(
                            list(width = "220px", targets = c(10, 11, 13, 14, 15, 16)),
                            list(width = "145px", targets = "_all"))),
                        class = "stripe hover compact")
        })
      } else {
        output[[id]] <- shiny::renderTable(utils::head(values(), 100L))
      }
    }
    explorer_table("variants_table", explorer_variants)
    explorer_table("clinvar_table", clinvar_variants)
    output$download_variants <- shiny::downloadHandler(
      filename = function() "sarek-inspected-variants.csv",
      content = function(file) utils::write.csv(explorer_variants(), file, row.names = FALSE, na = "")
    )
    output$download_clinvar <- shiny::downloadHandler(
      filename = function() {
        selected <- clinvar_view_file()
        if (NROW(selected) != 1L) return("clinvar-preview.csv")
        stem <- sub("\\.vcf(\\.gz)?$", "", basename(selected$file[[1]]))
        paste0(stem, "-", sarek_text(input$clinvar_record_mode, "all"), "-preview.csv")
      },
      content = function(file) {
        shiny::req(NROW(clinvar_view_file()) == 1L)
        shiny::req(!nzchar(clinvar_view_data()$error))
        utils::write.csv(clinvar_variants(), file, row.names = FALSE, na = "")
      }
    )
    output$download_clinvar_vcf <- shiny::downloadHandler(
      filename = function() {
        selected <- clinvar_view_file()
        if (NROW(selected) != 1L) return("clinvar.vcf.gz")
        basename(selected$file[[1]])
      },
      content = function(file) {
        selected <- clinvar_view_file()
        shiny::req(NROW(selected) == 1L)
        path <- sarek_explorer_file(
          selected$path[[1]], explorer_run()$output_dir
        )
        shiny::req(nzchar(path))
        if (!file.copy(path, file, overwrite = TRUE)) {
          stop("Could not copy the selected VCF for download.")
        }
      }
    )
    report_oncoplot <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$build_oncoplot, {
      run <- explorer_run()
      data <- tryCatch(
        sarek_report_oncoplot_data(
          explorer_catalog(), sarek_text(run$output_dir),
          sarek_results_read_manifest(run)
        ),
        error = function(error) list(message = conditionMessage(error))
      )
      if (!is.null(data$message)) {
        shiny::showNotification(data$message, type = "warning", duration = 12)
      }
      report_oncoplot(list(
        run_key = sarek_text(run$run_dir, sarek_text(run$output_dir)), data = data
      ))
    })

    report_ui <- shiny::reactive({
      run <- explorer_run()
      saved <- report_oncoplot()
      summary <- if (!is.null(saved) &&
                     identical(saved$run_key, sarek_text(run$run_dir, sarek_text(run$output_dir))))
        saved$data else NULL
      sarek_explorer_report_ui(
        run, explorer_catalog(), sarek_explorer_empty_variants(),
        sarek_results_read_manifest(run), oncoplot = summary
      )
    })

    output$report_view <- shiny::renderUI(report_ui())
    output$download_report <- shiny::downloadHandler(
      filename = function() paste0(sarek_text(explorer_run()$run_id, "sarek"), "-report.html"),
      content = function(file) writeLines(sarek_explorer_report_html(report_ui()), file, useBytes = TRUE)
    )

    shiny::observeEvent(input$open_selected_results, {
      run <- selected_run()
      shiny::updateSelectInput(session, "explorer_run", selected = sarek_text(run$run_id))
      shiny::updateTabsetPanel(session, "sarek_sections", selected = "explorer")
    }, ignoreInit = TRUE)

    invalidate_confirmation <- function() {
      confirmed_manifest(NULL)
      reviewed_validation_version(-1L)
      submission_state(list(
        kind = "info",
        message = "Run settings changed. Confirm the updated analysis before submission.",
        result = NULL
      ))
      invisible(NULL)
    }

    shiny::observeEvent(input$discover, {
      paths <- sarek_parse_path_input(input$paths)
      if (!length(paths)) {
        confirmation_state(NULL)
        invalidate_confirmation()
        status_state("ERROR: Enter at least one absolute server path.")
        validation_state("No manifest has been validated yet.")
        return(invisible(NULL))
      }
      tryCatch({
        table <- sarek_build_discovery_table(
          paths,
          recursive = isTRUE(input$recursive),
          max_files = max_files,
          allowed_roots = allowed_input_roots
        )
        if (is.function(bam_inspector)) {
          table <- sarek_auto_inspect_bams(
            table,
            inspector = bam_inspector,
            samtools = samtools,
            max_bams = max_auto_bam_inspections
          )
        }
        confirmation_state(table)
        sample_review <- sarek_sample_review_table(table)
        selected_sample_state(sarek_sample_key(sample_review$patient_id[[1]], sample_review$sample_id[[1]]))
        invalidate_confirmation()
        recommended_mode <- sarek_recommend_analysis_mode(table)
        if (nzchar(recommended_mode)) {
          shiny::updateSelectInput(session, "analysis_mode", selected = recommended_mode)
        }
        ignored <- as.integer(sarek_shiny_value(attr(table, "ignored_count"), 0L))
        bam_rows <- table[table$input_format == "bam", , drop = FALSE]
        inspection_note <- ""
        if (NROW(bam_rows)) {
          statuses <- base::table(bam_rows$inspection_status)
          status_text <- paste0(names(statuses), "=", as.integer(statuses), collapse = ", ")
          inspection_note <- paste0(" BAM inspection: ", status_text, ".")
        }
        status_state(paste0(
          "ACTION REQUIRED: Discovered ", NROW(table), " supported file", if (NROW(table) == 1L) "" else "s",
          if (ignored > 0L) paste0("; ignored ", ignored, " unsupported file", if (ignored == 1L) "" else "s") else "",
          ". Review every inferred sample field and FASTQ pairing before confirmation",
          if (nzchar(recommended_mode)) paste0(". Suggested analysis mode: ", sarek_analysis_mode_label(recommended_mode)) else "",
          ".", inspection_note
        ))
        validation_state("Discovery complete. The current draft has not been confirmed.")
      }, error = function(error) {
        confirmation_state(NULL)
        invalidate_confirmation()
        status_state(paste0("ERROR: ", conditionMessage(error)))
        validation_state("No manifest has been validated yet.")
      })
      invisible(NULL)
    }, ignoreInit = TRUE)

    output$status <- shiny::renderUI({
      message <- status_state()
      kind <- sarek_manifest_status_kind(message)
      shiny::div(
        class = paste("sarek-status-banner", paste0("sarek-status-banner-", kind)),
        shiny::tags$span(message)
      )
    })
    output$summary <- shiny::renderTable(
      sarek_confirmation_summary(confirmation_state()),
      striped = TRUE,
      bordered = TRUE,
      spacing = "s"
    )

    shiny::observeEvent(input$sample_to_edit, {
      selected_sample_state(as.character(sarek_shiny_value(input$sample_to_edit)))
    }, ignoreInit = TRUE)

    selected_sample_key <- shiny::reactive({
      table <- confirmation_state()
      if (is.null(table) || !NROW(table)) return("")
      keys <- unique(sarek_sample_key(table$patient_id, table$sample_id))
      selected <- selected_sample_state()
      if (selected %in% keys) selected else keys[[1]]
    })

    output$readiness <- shiny::renderUI({
      table <- confirmation_state()
      shiny::req(!is.null(table) && NROW(table))
      review <- sarek_sample_review_table(table)
      pairing_issues <- review[startsWith(review$FASTQ_pairing, "Needs attention:"), , drop = FALSE]
      recommended <- sarek_recommend_analysis_mode(table)
      selected_mode <- as.character(sarek_shiny_value(input$analysis_mode))

      pairing_item <- if (!NROW(pairing_issues)) {
        shiny::div(class = "sarek-review-item sarek-review-item-pass", "✓ FASTQ pairing is complete for every included FASTQ sample.")
      } else {
        details <- paste0(pairing_issues$sample_id, " (", sub("^Needs attention: ", "", pairing_issues$FASTQ_pairing), ")")
        shiny::div(
          class = "sarek-review-item sarek-review-item-review",
          paste0("! FASTQ files need attention: ", paste(details, collapse = "; "), ".")
        )
      }

      mode_item <- if (nzchar(recommended) && identical(selected_mode, recommended)) {
        shiny::div(
          class = "sarek-review-item sarek-review-item-pass",
          paste0("✓ Analysis mode matches the detected roles: ", sarek_analysis_mode_label(recommended), ".")
        )
      } else if (nzchar(recommended)) {
        shiny::div(
          class = "sarek-review-item sarek-review-item-review",
          paste0(
            "! Detected roles suggest ", sarek_analysis_mode_label(recommended),
            ", but the selected mode is ", sarek_analysis_mode_label(selected_mode), "."
          )
        )
      } else {
        shiny::div(
          class = "sarek-review-item sarek-review-item-review",
          "! Roles do not support one clear analysis mode yet; review each sample role."
        )
      }

      included_samples <- !is.na(review$include) & as.logical(review$include)
      indexed_samples <- included_samples & review$input_format %in% c("bam", "cram", "vcf", "bcf")
      missing_indexes <- indexed_samples & review$index == "Missing"
      index_item <- if (!any(missing_indexes)) {
        shiny::div(class = "sarek-review-item sarek-review-item-pass", "✓ Every required input index was detected.")
      } else {
        shiny::div(
          class = "sarek-review-item sarek-review-item-review",
          paste0("! Missing index: ", paste(review$sample_id[missing_indexes], collapse = ", "), ".")
        )
      }

      bam_samples <- included_samples & review$input_format == "bam"
      bam_verified <- bam_samples & startsWith(review$BAM_inspection, "Passed:")
      bam_item <- if (!any(bam_samples)) {
        shiny::div(class = "sarek-review-item sarek-review-item-pass", "✓ No BAM inspection is required for these inputs.")
      } else if (all(bam_verified[bam_samples])) {
        shiny::div(class = "sarek-review-item sarek-review-item-pass", "✓ Lightweight BAM integrity, header, and index inspection passed.")
      } else {
        shiny::div(
          class = "sarek-review-item sarek-review-item-review",
          paste0("! BAM inspection needs review: ", paste(review$sample_id[bam_samples & !bam_verified], collapse = ", "), ".")
        )
      }

      alignment_samples <- included_samples & review$input_format %in% c("bam", "cram")
      state_item <- if (!any(alignment_samples & review$processing_state == "unknown")) {
        shiny::div(class = "sarek-review-item sarek-review-item-pass", "✓ Processing state has been explicitly selected for every alignment input.")
      } else {
        shiny::div(
          class = "sarek-review-item sarek-review-item-review",
          paste0(
            "! Confirm processing state for: ",
            paste(review$sample_id[alignment_samples & review$processing_state == "unknown"], collapse = ", "),
            "."
          )
        )
      }

      shiny::div(
        class = "sarek-review-checklist",
        shiny::tags$strong("Pre-validation checklist"),
        pairing_item,
        index_item,
        bam_item,
        state_item,
        mode_item
      )
    })

    output$sample_editor <- shiny::renderUI({
      table <- confirmation_state()
      shiny::req(!is.null(table) && NROW(table))
      review <- sarek_sample_review_table(table)
      keys <- sarek_sample_key(review$patient_id, review$sample_id)
      labels <- paste0(review$patient_id, " / ", review$sample_id, " — ", review$role)
      choices <- stats::setNames(keys, labels)
      selected <- selected_sample_key()
      row <- review[keys == selected, , drop = FALSE]
      if (!NROW(row)) row <- review[1, , drop = FALSE]
      file_keys <- sarek_sample_key(table$patient_id, table$sample_id)
      sample_files <- table[file_keys == sarek_sample_key(row$patient_id[[1]], row$sample_id[[1]]), , drop = FALSE]
      recommendations <- if ("processing_recommendation" %in% names(sample_files)) {
        unique(trimws(as.character(sample_files$processing_recommendation)))
      } else character(0)
      recommendations <- recommendations[
        !is.na(recommendations) & recommendations %in% SAREK_PROCESSING_STATES & recommendations != "unknown"
      ]
      recommended_state <- if (length(recommendations) == 1L) recommendations[[1]] else ""
      current_state <- as.character(row$processing_state[[1]])
      selected_state <- if (identical(current_state, "unknown") && nzchar(recommended_state)) {
        recommended_state
      } else {
        current_state
      }

      normal_ids <- unique(review$sample_id[
        review$patient_id == row$patient_id[[1]] & review$role == "normal"
      ])
      current_match <- as.character(row$matched_normal_id[[1]])
      normal_ids <- unique(c(normal_ids, current_match[nzchar(current_match)]))
      match_choices <- c("No matched normal" = "", stats::setNames(normal_ids, normal_ids))

      shiny::div(
        class = "sarek-sample-editor",
        shiny::h4("Edit one sample"),
        shiny::tags$p(
          class = "muted small-note",
          "Select a sample, correct its metadata once, then apply the change to all associated files."
        ),
        shiny::selectInput(ns("sample_to_edit"), "Sample to review", choices = choices, selected = selected),
        shiny::fluidRow(
          shiny::column(
            6,
            shiny::checkboxInput(ns("edit_include"), "Include this sample", value = isTRUE(row$include[[1]])),
            shiny::textInput(ns("edit_patient_id"), "Patient ID", value = row$patient_id[[1]]),
            shiny::textInput(ns("edit_sample_id"), "Sample ID", value = row$sample_id[[1]]),
            shiny::selectInput(
              ns("edit_sex"),
              "Patient sex chromosomes",
              choices = c("Not provided" = "NA", "XX" = "XX", "XY" = "XY"),
              selected = row$sex[[1]],
              selectize = FALSE
            )
          ),
          shiny::column(
            6,
            shiny::selectInput(
              ns("edit_role"),
              "Sample role",
              choices = SAREK_SAMPLE_ROLES,
              selected = row$role[[1]],
              selectize = FALSE
            ),
            shiny::selectInput(
              ns("edit_matched_normal_id"),
              "Matched normal (tumor samples only)",
              choices = match_choices,
              selected = current_match
            ),
            shiny::selectInput(
              ns("edit_processing_state"),
              "Processing state",
              choices = SAREK_PROCESSING_STATES,
              selected = selected_state
            )
          )
        ),
        shiny::tags$p(
          class = "muted small-note",
          paste0(
            "Sex chromosomes apply to every sample for this patient. Detected format: ",
            row$input_format[[1]], ". FASTQ pairing: ", row$FASTQ_pairing[[1]], "."
          )
        ),
        if (identical(current_state, "unknown") && nzchar(recommended_state)) {
          shiny::tags$p(
            class = "sarek-storage-note",
            paste0(
              "BAM inspection recommends '", gsub("_", " ", recommended_state, fixed = TRUE),
              "'. It is preselected above but is not confirmed until you select Apply sample changes."
            )
          )
        },
        shiny::actionButton(ns("apply_sample"), "Apply sample changes", class = "btn-primary")
      )
    })

    output$bam_inspection <- shiny::renderUI({
      table <- confirmation_state()
      shiny::req(!is.null(table) && NROW(table))
      keys <- sarek_sample_key(table$patient_id, table$sample_id)
      bam <- table[keys == selected_sample_key() & table$input_format == "bam", , drop = FALSE]
      if (!NROW(bam)) return(NULL)
      row <- bam[1, , drop = FALSE]
      value <- function(field, default = "") {
        if (!field %in% names(row)) return(default)
        sarek_text(row[[field]], default)
      }
      status <- value("inspection_status", "not_run")
      panel_kind <- if (identical(status, "passed")) {
        "pass"
      } else if (identical(status, "failed")) {
        "failed"
      } else {
        "review"
      }
      index_path <- value("index")
      recommendation <- value("processing_recommendation", "unknown")
      confidence <- value("processing_confidence", "none")
      summary <- value("inspection_summary", "BAM inspection has not run.")
      sample_ids <- value("header_sample_ids", "Not reported")
      sort_order <- value("sort_order", "unknown")
      reference <- value("reference_compatibility", "Not inspected")
      evidence <- value("inspection_evidence")

      shiny::div(
        class = paste("sarek-bam-inspection", paste0("sarek-bam-inspection-", panel_kind)),
        shiny::h4("Lightweight BAM inspection"),
        shiny::tags$p(shiny::tags$strong(summary)),
        shiny::tags$dl(
          shiny::tags$dt("BAM"), shiny::tags$dd(value("path")),
          shiny::tags$dt("Index"), shiny::tags$dd(if (nzchar(index_path)) index_path else "Missing"),
          shiny::tags$dt("Header sample ID"), shiny::tags$dd(sample_ids),
          shiny::tags$dt("Sort order"), shiny::tags$dd(sort_order),
          shiny::tags$dt("GRCh38 evidence"), shiny::tags$dd(reference),
          shiny::tags$dt("Suggested state"),
          shiny::tags$dd(paste0(gsub("_", " ", recommendation, fixed = TRUE), " (", confidence, " confidence)"))
        ),
        if (nzchar(evidence)) shiny::tags$p(class = "muted small-note", evidence),
        shiny::tags$p(
          class = "muted small-note",
          "This check reads the BAM header and index metadata only; it does not scan reads or call variants."
        ),
        shiny::actionButton(ns("rerun_bam_inspection"), "Re-run BAM inspection")
      )
    })

    shiny::observeEvent(input$rerun_bam_inspection, {
      table <- confirmation_state()
      shiny::req(!is.null(table) && NROW(table))
      if (!is.function(bam_inspector)) {
        status_state("ERROR: BAM inspection is not configured for this app session.")
        return(invisible(NULL))
      }
      keys <- sarek_sample_key(table$patient_id, table$sample_id)
      rows <- which(keys == selected_sample_key() & table$input_format == "bam")
      if (!length(rows)) {
        status_state("ERROR: The selected sample does not contain a BAM file.")
        return(invisible(NULL))
      }
      tryCatch({
        for (row in rows) {
          inspection <- bam_inspector(
            path = as.character(table$path[[row]]),
            index = as.character(table$index[[row]]),
            samtools = samtools
          )
          table <- sarek_apply_bam_inspection(table, table$path[[row]], inspection)
        }
        confirmation_state(table)
        invalidate_confirmation()
        statuses <- unique(table$inspection_status[rows])
        if (all(statuses == "passed")) {
          status_state("ACTION REQUIRED: BAM inspection passed. Review its evidence and confirm the recommended processing state.")
        } else {
          status_state(paste0(
            "ACTION REQUIRED: BAM inspection requires review. Status: ",
            paste(statuses, collapse = ", "), "."
          ))
        }
        validation_state("The BAM inspection changed. Confirm the processing state, then validate the manifest again.")
      }, error = function(error) {
        status_state(paste0("ERROR: BAM inspection failed unexpectedly: ", conditionMessage(error)))
      })
      invisible(NULL)
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$apply_sample, {
      table <- confirmation_state()
      shiny::req(!is.null(table) && NROW(table))
      patient_id <- trimws(as.character(sarek_shiny_value(input$edit_patient_id)))
      sample_id <- trimws(as.character(sarek_shiny_value(input$edit_sample_id)))
      if (!nzchar(patient_id) || !nzchar(sample_id)) {
        status_state("ERROR: Patient ID and sample ID cannot be blank.")
        return(invisible(NULL))
      }
      tryCatch({
        table <- sarek_apply_sample_update(
          table = table,
          sample_key = selected_sample_key(),
          include = isTRUE(input$edit_include),
          patient_id = patient_id,
          sample_id = sample_id,
          role = input$edit_role,
          matched_normal_id = input$edit_matched_normal_id,
          processing_state = input$edit_processing_state,
          sex = input$edit_sex
        )
        confirmation_state(table)
        selected_sample_state(sarek_sample_key(patient_id, sample_id))
        invalidate_confirmation()
        status_state(paste0(
          "ACTION REQUIRED: Updated ", patient_id, " / ", sample_id,
          " across all associated files. Continue reviewing the checklist before validation."
        ))
        validation_state("The draft changed. Validate it again before downloading.")
      }, error = function(error) {
        status_state(paste0("ERROR: ", conditionMessage(error)))
      })
      invisible(NULL)
    }, ignoreInit = TRUE)

    output$file_editor <- shiny::renderUI({
      table <- confirmation_state()
      shiny::req(!is.null(table) && NROW(table))
      keys <- sarek_sample_key(table$patient_id, table$sample_id)
      sample_files <- table[keys == selected_sample_key(), , drop = FALSE]
      shiny::req(NROW(sample_files))
      fastq_files <- sample_files[sample_files$input_format == "fastq", , drop = FALSE]
      if (!NROW(fastq_files)) {
        return(shiny::tags$p(class = "muted small-note", "This sample does not use FASTQ input, so lane/read correction is not needed."))
      }

      choices <- stats::setNames(fastq_files$path, basename(fastq_files$path))
      selected_path <- as.character(sarek_shiny_value(input$file_to_edit, fastq_files$path[[1]]))
      if (!selected_path %in% fastq_files$path) selected_path <- fastq_files$path[[1]]
      row <- fastq_files[fastq_files$path == selected_path, , drop = FALSE]
      selected_read <- if (is.na(row$read[[1]])) "" else as.character(row$read[[1]])

      shiny::div(
        class = "sarek-sample-editor",
        shiny::h4("Correct FASTQ lane/read detection"),
        shiny::tags$p(
          class = "muted small-note",
          "Use this only when a valid FASTQ filename was not interpreted correctly. It cannot supply a genuinely missing mate."
        ),
        shiny::selectInput(ns("file_to_edit"), "FASTQ file", choices = choices, selected = selected_path),
        shiny::fluidRow(
          shiny::column(6, shiny::textInput(ns("edit_lane"), "Lane", value = row$lane[[1]], placeholder = "L001")),
          shiny::column(
            6,
            shiny::selectInput(
              ns("edit_read"),
              "Read",
              choices = c("Not detected" = "", "R1" = "1", "R2" = "2"),
              selected = selected_read
            )
          )
        ),
        shiny::actionButton(ns("apply_file_pairing"), "Apply lane/read correction")
      )
    })

    shiny::observeEvent(input$apply_file_pairing, {
      table <- confirmation_state()
      shiny::req(!is.null(table) && NROW(table))
      tryCatch({
        read <- if (nzchar(as.character(sarek_shiny_value(input$edit_read)))) input$edit_read else NA_integer_
        table <- sarek_apply_file_pairing_update(
          table,
          path = input$file_to_edit,
          lane = input$edit_lane,
          read = read
        )
        confirmation_state(table)
        invalidate_confirmation()
        status_state("ACTION REQUIRED: FASTQ lane/read correction applied. Recheck pairing before validation.")
        validation_state("The draft changed. Validate it again before downloading.")
      }, error = function(error) {
        status_state(paste0("ERROR: ", conditionMessage(error)))
      })
      invisible(NULL)
    }, ignoreInit = TRUE)

    if (requireNamespace("DT", quietly = TRUE)) {
      output$sample_review_table <- DT::renderDT({
        table <- sarek_sample_review_display_table(confirmation_state())
        shiny::validate(shiny::need(NROW(table), "Discover inputs to populate this table."))
        DT::datatable(
          table,
          rownames = FALSE,
          class = "stripe hover compact nowrap",
          width = "100%",
          selection = list(mode = "single", target = "row"),
          options = list(
            scrollX = TRUE,
            scrollCollapse = FALSE,
            pageLength = 25,
            lengthMenu = list(c(10, 25, 50, 100), c("10", "25", "50", "100")),
            searching = TRUE,
            autoWidth = TRUE
          )
        )
      }, server = FALSE)

      shiny::observeEvent(input$sample_review_table_rows_selected, {
        selected_rows <- input$sample_review_table_rows_selected
        if (!length(selected_rows)) return(invisible(NULL))
        selected_row <- suppressWarnings(as.integer(selected_rows[[1]]))
        review <- sarek_sample_review_table(confirmation_state())
        if (!is.na(selected_row) && selected_row >= 1L && selected_row <= NROW(review)) {
          selected_sample_state(sarek_sample_key(
            review$patient_id[[selected_row]],
            review$sample_id[[selected_row]]
          ))
        }
        invisible(NULL)
      }, ignoreInit = TRUE)

      output$file_detail_table <- DT::renderDT({
        table <- confirmation_state()
        shiny::validate(shiny::need(!is.null(table) && NROW(table), "Discover inputs to populate file details."))
        keys <- sarek_sample_key(table$patient_id, table$sample_id)
        table <- table[keys == selected_sample_key(), c(
          "path", "index", "lane", "read", "size_bytes", "role_confidence", "warning"
        ), drop = FALSE]
        DT::datatable(
          table,
          rownames = FALSE,
          class = "stripe hover compact nowrap",
          width = "100%",
          selection = "none",
          options = list(
            scrollX = TRUE,
            scrollCollapse = FALSE,
            paging = FALSE,
            searching = FALSE,
            info = FALSE,
            autoWidth = TRUE,
            columnDefs = sarek_manifest_column_defs(table)
          )
        )
      }, server = FALSE)
    } else {
      output$sample_review_table <- shiny::renderTable({
        sarek_sample_review_display_table(confirmation_state())
      }, striped = TRUE, bordered = TRUE, spacing = "xs")
      output$file_detail_table <- shiny::renderTable({
        table <- confirmation_state()
        shiny::req(!is.null(table) && NROW(table))
        keys <- sarek_sample_key(table$patient_id, table$sample_id)
        table[keys == selected_sample_key(), c(
          "path", "index", "lane", "read", "size_bytes", "role_confidence", "warning"
        ), drop = FALSE]
      }, striped = TRUE, bordered = TRUE, spacing = "xs")
    }

    shiny::observeEvent(
      list(
        input$manifest_id,
        input$assay_type,
        input$analysis_mode,
        input$preset,
        input$results_root,
        input$work_root,
        input$species,
        input$assembly,
        input$sarek_genome
      ),
      {
        if (!is.null(confirmed_manifest())) {
          invalidate_confirmation()
          validation_state("The draft settings changed. Validate them again before downloading.")
        }
      },
      ignoreInit = TRUE
    )

    shiny::observeEvent(input$confirm, {
      table <- confirmation_state()
      if (is.null(table) || !NROW(table)) {
        confirmed_manifest(NULL)
        validation_state("ERRORS\n- Discover and review inputs before confirming a manifest.")
        status_state("ERROR: Manifest validation could not start. Discover and review inputs first.")
        return(invisible(NULL))
      }
      validation <- sarek_validate_confirmation_table(table)
      mode_errors <- sarek_validate_analysis_mode(table, input$analysis_mode)
      validation_state(sarek_manifest_validation_text(validation, mode_errors))
      if (!isTRUE(validation$valid) || length(mode_errors)) {
        confirmed_manifest(NULL)
        status_state("ERROR: Manifest validation failed. Review the red validation panel and correct every listed item.")
        return(invisible(NULL))
      }
      tryCatch({
        manifest <- sarek_build_manifest(
          confirmation_table = table,
          manifest_id = input$manifest_id,
          created_by = created_by,
          assay_type = input$assay_type,
          analysis_mode = input$analysis_mode,
          preset = input$preset,
          results_root = input$results_root,
          work_root = input$work_root,
          species = input$species,
          assembly = input$assembly,
          sarek_genome = input$sarek_genome,
          allowed_results_roots = allowed_results_roots,
          allowed_work_roots = allowed_work_roots
        )
        tool_resolution <- sarek_submission_tool_resolution(
          manifest$analysis$mode,
          manifest$analysis$preset,
          manifest
        )
        confirmation_warnings <- unique(c(validation$warnings, tool_resolution$warnings))
        confirmation_validation <- validation
        confirmation_validation$warnings <- confirmation_warnings
        confirmed_manifest(manifest)
        validation_version(validation_version() + 1L)
        reviewed_validation_version(-1L)
        submission_state(list(
          kind = "info",
          message = "The confirmed run is ready for final review.",
          result = NULL
        ))
        validation_state(paste0(
          sarek_manifest_validation_text(confirmation_validation, mode_errors),
          "\n\nCONFIRMED\n- Review the run summary and submission settings below."
        ))
        if (length(confirmation_warnings)) {
          status_state("ACTION REQUIRED: Validation passed with warnings. Read the amber panels before submitting this run.")
        } else {
          status_state("Manifest validation passed. Review the confirmed run before submission.")
        }
      }, error = function(error) {
        confirmed_manifest(NULL)
        validation_state(paste0("ERRORS\n- ", conditionMessage(error)))
        status_state("ERROR: Manifest validation failed. Review the red validation panel and correct every listed item.")
      })
      invisible(NULL)
    }, ignoreInit = TRUE)

    output$validation <- shiny::renderUI({
      message <- validation_state()
      kind <- sarek_manifest_validation_kind(message)
      shiny::div(
        class = paste("sarek-validation-panel", paste0("sarek-validation-panel-", kind)),
        shiny::div(
          shiny::h4(sarek_manifest_validation_heading(kind)),
          shiny::tags$pre(message)
        )
      )
    })

    output$nextflow_input <- shiny::renderUI({
      state <- nextflow_input_state()
      if (isTRUE(state$pending)) {
        return(shiny::div(
          class = "sarek-validation-panel",
          shiny::div(
            shiny::h4("Nextflow input has not been generated"),
            shiny::tags$pre("Confirm a valid manifest to generate the Sarek samplesheet and starting step.")
          )
        ))
      }
      if (!isTRUE(state$valid)) {
        return(shiny::div(
          class = "sarek-validation-panel sarek-validation-panel-error",
          shiny::div(
            shiny::h4("Nextflow input requires correction"),
            shiny::tags$pre(paste0("ERRORS\n- ", state$error))
          )
        ))
      }

      nextflow_input <- state$input
      warnings <- as.character(nextflow_input$warnings)
      warnings <- warnings[!is.na(warnings) & nzchar(warnings)]
      kind <- if (length(warnings)) "warning" else "success"
      message <- c(
        "READY",
        paste0("- Pipeline: ", nextflow_input$pipeline, " ", nextflow_input$pipeline_version),
        paste0("- Starting step: ", nextflow_input$step),
        paste0("- Input format: ", toupper(nextflow_input$input_format)),
        paste0("- Samplesheet rows: ", NROW(nextflow_input$samplesheet))
      )
      if (length(warnings)) {
        message <- c(message, "", "WARNINGS", paste0("- ", warnings))
      }
      shiny::div(
        class = paste("sarek-validation-panel", paste0("sarek-validation-panel-", kind)),
        shiny::div(
          shiny::h4(if (length(warnings)) "Nextflow input ready — review warning" else "Nextflow input ready"),
          shiny::tags$pre(paste(message, collapse = "\n"))
        )
      )
    })

    output$run_review <- shiny::renderUI({
      manifest <- confirmed_manifest()
      state <- nextflow_input_state()
      if (is.null(manifest)) {
        return(shiny::div(
          class = "sarek-run-review",
          shiny::tags$p(class = "muted small-note", "Confirm a valid analysis to display its final run summary.")
        ))
      }
      if (!isTRUE(state$valid)) {
        return(shiny::div(
          class = "sarek-validation-panel sarek-validation-panel-error",
          shiny::div(shiny::h4("Run review is unavailable"), shiny::tags$pre(state$error))
        ))
      }
      samples <- sarek_confirmed_samples_table(manifest)
      tool_resolution <- tryCatch(
        sarek_submission_tool_resolution(
          manifest$analysis$mode,
          manifest$analysis$preset,
          manifest
        ),
        error = function(error) list(
          tools = paste("Not available:", conditionMessage(error)),
          warnings = character(0)
        )
      )
      card <- function(label, value) {
        shiny::div(
          class = "sarek-run-review-card",
          shiny::tags$span(label),
          shiny::tags$strong(as.character(value))
        )
      }
      shiny::div(
        class = "sarek-run-review",
        shiny::div(
          class = "sarek-run-review-grid",
          card("Run name", manifest$manifest_id),
          card("Analysis", sarek_analysis_mode_label(manifest$analysis$mode)),
          card("Assay", manifest$assay$type),
          card("Reference", paste(manifest$reference$assembly, manifest$reference$sarek_genome, sep = " · ")),
          card("Patients", length(manifest$patients)),
          card("Samples", NROW(samples)),
          card("Starting stage", state$input$step),
          card("Pipeline tools", tool_resolution$tools)
        ),
        if (length(tool_resolution$warnings)) {
          shiny::div(
            class = "sarek-validation-panel sarek-validation-panel-warning",
            shiny::div(
              shiny::h4("Caller adjusted - run will continue"),
              shiny::tags$pre(paste(tool_resolution$warnings, collapse = "\n"))
            )
          )
        },
        shiny::tags$p(
          class = "muted small-note",
          shiny::tags$strong("Results: "), manifest$storage$results_root,
          shiny::tags$br(),
          shiny::tags$strong("Temporary work: "), manifest$storage$work_root
        )
      )
    })

    output$confirmed_samples <- shiny::renderTable({
      manifest <- confirmed_manifest()
      if (is.null(manifest)) return(NULL)
      sarek_confirmed_samples_table(manifest)
    }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = FALSE)

    output$submission_ui <- shiny::renderUI({
      manifest <- confirmed_manifest()
      if (is.null(manifest)) {
        return(shiny::div(
          class = "sarek-submit-panel",
          shiny::tags$p("Validate and confirm the analysis before submission becomes available.")
        ))
      }
      nextflow_state <- nextflow_input_state()
      if (!isTRUE(nextflow_state$valid)) {
        return(shiny::div(
          class = "sarek-submit-panel",
          shiny::tags$p("Submission is locked until the Nextflow input error is corrected.")
        ))
      }
      readiness <- submission_readiness()
      if (!isTRUE(readiness$valid)) {
        return(shiny::div(
          class = "sarek-validation-panel sarek-validation-panel-error",
          shiny::div(
            shiny::h4("Submission requires another setting"),
            shiny::tags$pre(readiness$error)
          )
        ))
      }
      state <- submission_state()
      if (identical(state$kind, "success")) {
        return(shiny::div(
          class = "sarek-submit-panel",
          shiny::tags$strong("This run has been submitted."),
          shiny::tags$p(class = "muted small-note", "Use the recorded controller PID and run activity to track this analysis.")
        ))
      }
      reviewed <- identical(reviewed_validation_version(), validation_version())
      shiny::div(
        class = "sarek-submit-panel",
        shiny::checkboxInput(
          session$ns("validation_reviewed"),
          "I reviewed the validation, warnings, run settings, samples, and storage paths above",
          value = reviewed
        ),
        if (reviewed && !isTRUE(submission_in_progress())) {
          shiny::tagList(
            shiny::actionButton(session$ns("request_submit"), "Submit Sarek run", class = "btn-primary"),
            shiny::tags$details(
              class = "sarek-advanced-downloads",
              shiny::tags$summary("Advanced: download generated samplesheet"),
              shiny::br(),
              shiny::downloadButton(session$ns("download_nextflow_input"), "Download Sarek samplesheet CSV")
            )
          )
        } else if (isTRUE(submission_in_progress())) {
          shiny::tags$p(class = "muted small-note", "Starting the detached Nextflow controller...")
        } else {
          shiny::tags$p(class = "muted small-note", "Acknowledge the review above to enable submission.")
        }
      )
    })

    output$submission_status <- shiny::renderUI({
      state <- submission_state()
      if (identical(state$kind, "info")) return(NULL)
      shiny::div(
        class = paste("sarek-validation-panel", paste0("sarek-validation-panel-", state$kind)),
        shiny::div(
          shiny::h4(if (identical(state$kind, "success")) "Run submitted" else "Submission failed"),
          shiny::tags$pre(state$message)
        )
      )
    })

    shiny::observeEvent(input$validation_reviewed, {
      if (isTRUE(input$validation_reviewed) && !is.null(confirmed_manifest())) {
        reviewed_validation_version(validation_version())
      } else {
        reviewed_validation_version(-1L)
      }
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$request_submit, {
      manifest <- confirmed_manifest()
      state <- nextflow_input_state()
      readiness <- submission_readiness()
      if (is.null(manifest) || !isTRUE(state$valid) || !isTRUE(readiness$valid) ||
          !identical(reviewed_validation_version(), validation_version())) return(invisible(NULL))
      tool_resolution <- tryCatch(
        sarek_submission_tool_resolution(
          manifest$analysis$mode,
          manifest$analysis$preset,
          manifest
        ),
        error = function(error) list(tools = conditionMessage(error), warnings = character(0))
      )
      shiny::showModal(shiny::modalDialog(
        title = "Submit this Sarek run?",
        shiny::tags$p("This creates a private run bundle and submits the Nextflow controller to Slurm."),
        shiny::tags$ul(
          shiny::tags$li(shiny::tags$strong("Run: "), manifest$manifest_id),
          shiny::tags$li(shiny::tags$strong("Starting stage: "), state$input$step),
          shiny::tags$li(shiny::tags$strong("Tools that will run: "), tool_resolution$tools),
          shiny::tags$li(shiny::tags$strong("Results root: "), manifest$storage$results_root),
          shiny::tags$li(shiny::tags$strong("Work root: "), manifest$storage$work_root)
        ),
        if (length(tool_resolution$warnings)) {
          shiny::div(
            class = "sarek-validation-panel sarek-validation-panel-warning",
            shiny::div(
              shiny::h4("Important caller change"),
              shiny::tags$pre(paste(tool_resolution$warnings, collapse = "\n"))
            )
          )
        },
        footer = shiny::tagList(
          shiny::modalButton("Cancel"),
          shiny::actionButton(session$ns("confirm_submit"), "Confirm and submit", class = "btn-primary")
        ),
        easyClose = FALSE
      ))
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$confirm_submit, {
      shiny::removeModal()
      if (isTRUE(submission_in_progress())) return(invisible(NULL))
      submission_in_progress(TRUE)
      on.exit(submission_in_progress(FALSE), add = TRUE)
      result <- tryCatch({
        manifest <- confirmed_manifest()
        state <- nextflow_input_state()
        readiness <- submission_readiness()
        if (is.null(manifest)) stop("The run is no longer confirmed.")
        if (!isTRUE(state$valid)) stop("The generated Sarek input is no longer valid: ", state$error)
        if (!isTRUE(readiness$valid)) stop(readiness$error)
        if (!identical(reviewed_validation_version(), validation_version())) {
          stop("Review acknowledgement is required before submission.")
        }
        if (!is.function(submit_handler)) stop("The Sarek submission backend is not configured.")
        submitted <- submit_handler(manifest, state$input)
        if (!is.list(submitted) || !identical(sarek_text(submitted$status), "submitted") ||
            !nzchar(sarek_text(submitted$controller_pid))) {
          stop("The Sarek submission backend did not return a detached controller PID.")
        }
        submitted
      }, error = function(error) error)
      if (inherits(result, "error")) {
        submission_state(list(
          kind = "error",
          message = paste0("ERROR\n- ", conditionMessage(result)),
          result = NULL
        ))
      } else {
        submission_state(list(
          kind = "success",
          message = paste0(
            "SUBMITTED\n- Nextflow controller PID: ", result$controller_pid,
            "\n- Run directory: ", result$run_dir,
            "\n- Results directory: ", result$output_dir
          ),
          result = result
        ))
        shiny::updateTextInput(session, "history_root", value = dirname(result$run_dir))
        run_history_refresh(Sys.time())
        shiny::updateSelectInput(session, "selected_run", selected = basename(result$run_dir))
        shiny::updateTabsetPanel(session, "sarek_sections", selected = "status")
      }
      invisible(NULL)
    }, ignoreInit = TRUE)

    output$download_nextflow_input <- shiny::downloadHandler(
      filename = function() {
        paste0(sarek_identifier(input$manifest_id, "sarek_analysis"), ".sarek.samplesheet.csv")
      },
      content = function(file) {
        if (!identical(reviewed_validation_version(), validation_version())) {
          stop("Review and acknowledge the validation results before downloading.")
        }
        state <- nextflow_input_state()
        if (!isTRUE(state$valid)) {
          stop("A valid Sarek Nextflow input is not available: ", state$error)
        }
        sarek_write_nextflow_samplesheet(state$input, file)
      },
      contentType = "text/csv"
    )

    list(
      confirmation_table = shiny::reactive(confirmation_state()),
      manifest = shiny::reactive(confirmed_manifest()),
      valid = shiny::reactive(!is.null(confirmed_manifest())),
      nextflow_input = shiny::reactive(nextflow_input_state()$input),
      nextflow_valid = shiny::reactive(isTRUE(nextflow_input_state()$valid)),
      submission_ready = shiny::reactive(isTRUE(submission_readiness()$valid)),
      submission = shiny::reactive(submission_state()$result)
    )
  })
}
