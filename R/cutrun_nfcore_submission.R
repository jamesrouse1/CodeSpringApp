CUTRUN_NFCORE_VERSION <- "3.2.2"

cutrun_nfcore_group_name <- function(row, control = FALSE) {
  fields <- c(
    trimws(as.character(row$cell_type %||% "")),
    trimws(as.character(row$mark %||% "")),
    trimws(as.character(row$condition %||% ""))
  )
  fields <- fields[nzchar(fields)]

  if (!length(fields)) {
    fields <- trimws(as.character(row$sample %||% "sample"))
  }

  value <- paste(fields, collapse = "_")
  value <- gsub("[^A-Za-z0-9_.-]+", "_", value)
  value <- gsub("^_+|_+$", "", value)

  prefix <- if (isTRUE(control)) "control" else "target"
  paste0(prefix, "__", value)
}


cutrun_nfcore_resolve_fastq <- function(project, value) {
  value <- trimws(as.character(value %||% ""))
  if (!nzchar(value)) stop("Empty FASTQ path encountered.")

  if (
    grepl(",", value, fixed = TRUE) ||
    grepl("\\n", value, fixed = TRUE) ||
    grepl("\\r", value, fixed = TRUE)
  ) {
    stop("FASTQ paths used by nf-core/cutandrun cannot contain commas or line breaks: ", value)
  }

  expanded <- path.expand(value)

  if (startsWith(expanded, "/")) {
    if (!file.exists(expanded)) {
      stop("FASTQ file does not exist: ", expanded)
    }
    path <- normalizePath(expanded, winslash = "/", mustWork = TRUE)
  } else {
    roots <- unique(c(
      dirname(project$design_matrix_path %||% ""),
      project_fastq_dirs(project)
    ))
    roots <- roots[nzchar(roots)]

    direct <- file.path(roots, value)
    direct <- direct[file.exists(direct)]

    if (length(direct) == 1L) {
      path <- normalizePath(direct[[1]], winslash = "/", mustWork = TRUE)
    } else if (length(direct) > 1L) {
      stop("Relative FASTQ path is ambiguous across configured input folders: ", value)
    } else {
      all_fastqs <- unique(unlist(
        lapply(project_fastq_dirs(project), fastq_files),
        use.names = FALSE
      ))

      basename_matches <- all_fastqs[basename(all_fastqs) == basename(value)]

      if (length(basename_matches) == 1L) {
        path <- normalizePath(
          basename_matches[[1]],
          winslash = "/",
          mustWork = TRUE
        )
      } else if (!length(basename_matches)) {
        stop("Could not resolve FASTQ file from project inputs: ", value)
      } else {
        stop(
          "FASTQ basename is ambiguous across configured input folders: ",
          basename(value)
        )
      }
    }
  }

  if (!grepl("\\.(fastq|fq)\\.gz$", path, ignore.case = TRUE)) {
    stop(
      "nf-core/cutandrun requires gzipped FASTQ files ending in .fastq.gz or .fq.gz: ",
      path
    )
  }

  if (grepl("[[:space:]]", path)) {
    stop(
      "nf-core/cutandrun does not accept FASTQ paths containing spaces: ",
      path
    )
  }

  path
}


cutrun_nfcore_lane_pairs <- function(project, filename, sample) {
  lanes <- trimws(unlist(
    strsplit(as.character(filename %||% ""), ";", fixed = TRUE),
    use.names = FALSE
  ))
  lanes <- lanes[nzchar(lanes)]

  if (!length(lanes)) {
    stop("No FASTQs are recorded for CUT&RUN sample: ", sample)
  }

  rows <- lapply(seq_along(lanes), function(i) {
    parts <- trimws(unlist(
      strsplit(lanes[[i]], ",", fixed = TRUE),
      use.names = FALSE
    ))
    parts <- parts[nzchar(parts)]

    if (length(parts) != 2L) {
      stop(
        "nf-core/cutandrun requires paired-end R1,R2 for every lane. ",
        "Sample ", sample, ", lane ", i,
        " contains ", length(parts), " FASTQ entr",
        if (length(parts) == 1L) "y." else "ies."
      )
    }

    data.frame(
      lane = i,
      fastq_1 = cutrun_nfcore_resolve_fastq(project, parts[[1]]),
      fastq_2 = cutrun_nfcore_resolve_fastq(project, parts[[2]]),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}


cutrun_nfcore_samplesheet <- function(project, design = NULL) {
  if (is.null(project) || !is_cutrun_project(project)) {
    stop("nf-core/cutandrun input generation requires a CUT&RUN project.")
  }

  if (!isTRUE(project$paired_end)) {
    stop("nf-core/cutandrun 3.2.2 supports paired-end CUT&RUN data only.")
  }

  if (is.null(design)) {
    design <- included_design_table(project)
  }

  if (!NROW(design)) {
    stop("The CUT&RUN design matrix has no included samples.")
  }

  if (!all(c("sample", "filename") %in% names(design))) {
    stop("CUT&RUN design requires sample and filename columns.")
  }

  design <- infer_cutrun_metadata(design)

  required_metadata <- c(
    "cell_type", "mark", "target_class",
    "condition", "replicate", "control_sample"
  )

  for (column in required_metadata) {
    if (!column %in% names(design)) design[[column]] <- ""
  }

  design$sample <- trimws(as.character(design$sample))

  if (any(!nzchar(design$sample))) {
    stop("Every CUT&RUN design row must have a sample name.")
  }

  if (anyDuplicated(design$sample)) {
    stop("CUT&RUN sample names must be unique before nf-core conversion.")
  }

  replicate <- suppressWarnings(as.integer(as.character(design$replicate)))

  bad_replicate <- is.na(replicate) | replicate <= 0L
  if (any(bad_replicate)) {
    stop(
      "nf-core/cutandrun requires positive integer replicate numbers. Fix: ",
      paste(design$sample[bad_replicate], collapse = ", ")
    )
  }

  design$replicate <- replicate

  control_rows <- tolower(trimws(as.character(design$target_class))) == "control"

  groups <- vapply(seq_len(NROW(design)), function(i) {
    cutrun_nfcore_group_name(
      design[i, , drop = FALSE],
      control = control_rows[[i]]
    )
  }, character(1))

  names(groups) <- design$sample

  control_groups <- character(NROW(design))

  for (i in seq_len(NROW(design))) {
    if (control_rows[[i]]) {
      control_groups[[i]] <- ""
      next
    }

    control_sample <- trimws(
      as.character(design$control_sample[[i]] %||% "")
    )

    if (!nzchar(control_sample)) {
      stop(
        "Every nf-core CUT&RUN target must identify its matched control_sample. Missing for: ",
        design$sample[[i]]
      )
    }

    control_index <- match(control_sample, design$sample)

    if (is.na(control_index)) {
      stop(
        "Control sample '", control_sample,
        "' referenced by ", design$sample[[i]],
        " is not present in the included design."
      )
    }

    if (!control_rows[[control_index]]) {
      stop(
        "control_sample '", control_sample,
        "' referenced by ", design$sample[[i]],
        " is not marked target_class=control."
      )
    }

    control_groups[[i]] <- groups[[control_index]]
  }

  output_rows <- list()
  map_rows <- list()

  for (i in seq_len(NROW(design))) {
    lanes <- cutrun_nfcore_lane_pairs(
      project,
      design$filename[[i]],
      design$sample[[i]]
    )

    for (lane_i in seq_len(NROW(lanes))) {
      output_rows[[length(output_rows) + 1L]] <- data.frame(
        group = groups[[i]],
        replicate = design$replicate[[i]],
        fastq_1 = lanes$fastq_1[[lane_i]],
        fastq_2 = lanes$fastq_2[[lane_i]],
        control = control_groups[[i]],
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      map_rows[[length(map_rows) + 1L]] <- data.frame(
        codespring_sample = design$sample[[i]],
        nfcore_group = groups[[i]],
        replicate = design$replicate[[i]],
        lane = lanes$lane[[lane_i]],
        role = if (control_rows[[i]]) "control" else "target",
        control_sample = if (control_rows[[i]]) "" else trimws(as.character(design$control_sample[[i]])),
        control_group = control_groups[[i]],
        fastq_1 = lanes$fastq_1[[lane_i]],
        fastq_2 = lanes$fastq_2[[lane_i]],
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }

  samplesheet <- do.call(rbind, output_rows)
  mapping <- do.call(rbind, map_rows)

  if (any(samplesheet$group == samplesheet$control & nzchar(samplesheet$control))) {
    stop("Internal error: nf-core target group and control group cannot be identical.")
  }

  list(
    samplesheet = samplesheet,
    mapping = mapping,
    design = design
  )
}


cutrun_nfcore_write_inputs <- function(
  project,
  run_dir = file.path(project$data_dir, "nfcore_cutandrun")
) {
  generated <- cutrun_nfcore_samplesheet(project)

  metadata_dir <- file.path(run_dir, ".codespring")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

  samplesheet_path <- file.path(metadata_dir, "samplesheet.csv")
  mapping_path <- file.path(metadata_dir, "sample_mapping.tsv")

  utils::write.table(
    generated$samplesheet,
    samplesheet_path,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = ""
  )

  utils::write.table(
    generated$mapping,
    mapping_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = ""
  )

  list(
    samplesheet = samplesheet_path,
    mapping = mapping_path,
    run_dir = normalizePath(run_dir, winslash = "/", mustWork = FALSE),
    rows = NROW(generated$samplesheet),
    samples = length(unique(generated$mapping$codespring_sample))
  )
}
