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

# -------------------------------------------------------------------------
# nf-core/cutandrun submission bundle
# -------------------------------------------------------------------------

CUTRUN_NFCORE_NEXTFLOW_VERSION <- "24.04.4"
CUTRUN_NFCORE_PIPELINE <- "nf-core/cutandrun"

cutrun_nfcore_runtime_defaults <- function() {
  backend_root <- Sys.getenv(
    "CSL_CUTRUN_NFCORE_BACKEND_ROOT",
    unset = "/grid/bsr/data/data/bsr_readable_data/CodeSpringFlow"
  )

  list(
    launcher = Sys.getenv(
      "CSL_CUTRUN_NFCORE_LAUNCHER",
      unset = file.path(backend_root, "bin", "nextflow-sarek")
    ),
    config = Sys.getenv(
      "CSL_CUTRUN_NFCORE_CONFIG",
      unset = file.path(backend_root, "conf", "cshl_slurm.config")
    ),
    nxf_home = Sys.getenv(
      "CSL_CUTRUN_NFCORE_NXF_HOME",
      unset = file.path(
        backend_root,
        "runtime",
        "nextflow",
        "cutandrun"
      )
    ),
    singularity_cache = Sys.getenv(
      "CSL_CUTRUN_NFCORE_SINGULARITY_CACHE",
      unset = file.path(
        backend_root,
        "cache",
        "singularity"
      )
    )
  )
}

cutrun_nfcore_validate_runtime <- function(runtime) {
  required <- c(
    "launcher",
    "config",
    "nxf_home",
    "singularity_cache"
  )

  missing <- required[
    !vapply(
      required,
      function(name) {
        value <- runtime[[name]]
        length(value) &&
          !is.na(value[[1]]) &&
          nzchar(trimws(as.character(value[[1]])))
      },
      logical(1)
    )
  ]

  if (length(missing)) {
    stop(
      "CUT&RUN nf-core runtime is missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  runtime <- lapply(runtime, function(value) {
    trimws(as.character(value[[1]]))
  })

  if (!grepl("^/", runtime$launcher)) {
    stop(
      "CUT&RUN nf-core launcher must use an absolute path.",
      call. = FALSE
    )
  }

  if (!file.exists(runtime$launcher) ||
      dir.exists(runtime$launcher) ||
      file.access(runtime$launcher, 1) != 0) {
    stop(
      "CUT&RUN nf-core launcher is missing or not executable: ",
      runtime$launcher,
      call. = FALSE
    )
  }

  if (!grepl("^/", runtime$config) ||
      !file.exists(runtime$config) ||
      dir.exists(runtime$config) ||
      file.access(runtime$config, 4) != 0) {
    stop(
      "CUT&RUN nf-core CSHL config is missing or unreadable: ",
      runtime$config,
      call. = FALSE
    )
  }

  for (field in c("nxf_home", "singularity_cache")) {
    if (!grepl("^/", runtime[[field]])) {
      stop(
        field,
        " must use an absolute path.",
        call. = FALSE
      )
    }

    if (!dir.create(
      runtime[[field]],
      recursive = TRUE,
      showWarnings = FALSE
    ) && !dir.exists(runtime[[field]])) {
      stop(
        "Could not create CUT&RUN nf-core runtime directory: ",
        runtime[[field]],
        call. = FALSE
      )
    }
  }

  runtime
}

cutrun_nfcore_submission_paths <- function(
  project,
  run_id = "nfcore_cutandrun"
) {
  if (is.null(project) ||
      !is.list(project) ||
      !is_cutrun_project(project)) {
    stop(
      "CUT&RUN nf-core submission requires a CUT&RUN project.",
      call. = FALSE
    )
  }

  data_dir <- trimws(as.character(project$data_dir %||% "")[1])

  if (is.na(data_dir) ||
      !nzchar(data_dir) ||
      !grepl("^/", data_dir)) {
    stop(
      "CUT&RUN project data_dir must be an absolute path.",
      call. = FALSE
    )
  }

  if (!dir.exists(data_dir)) {
    stop(
      "CUT&RUN project data directory does not exist: ",
      data_dir,
      call. = FALSE
    )
  }

  run_id <- trimws(as.character(run_id %||% "")[1])
  run_id <- gsub("[^A-Za-z0-9._-]+", "_", run_id)
  run_id <- gsub("^_+|_+$", "", run_id)

  if (!nzchar(run_id)) {
    stop(
      "CUT&RUN nf-core run ID cannot be empty.",
      call. = FALSE
    )
  }

  run_dir <- normalizePath(
    file.path(data_dir, run_id),
    winslash = "/",
    mustWork = FALSE
  )

  internal_dir <- file.path(run_dir, ".codespring")

  list(
    run_id = run_id,
    run_dir = run_dir,
    output_dir = file.path(run_dir, "results"),
    internal_dir = internal_dir,
    log_dir = file.path(internal_dir, "logs"),
    work_dir = file.path(internal_dir, "work"),
    samplesheet_path = file.path(
      internal_dir,
      "samplesheet.csv"
    ),
    mapping_path = file.path(
      internal_dir,
      "sample_mapping.tsv"
    ),
    params_path = file.path(
      internal_dir,
      "params.json"
    ),
    run_config = file.path(
      internal_dir,
      "nextflow.config"
    ),
    launch_script = file.path(
      internal_dir,
      "launch.sh"
    ),
    nextflow_log = file.path(
      internal_dir,
      "logs",
      "nextflow.log"
    ),
    trace_path = file.path(
      internal_dir,
      "logs",
      "trace.tsv"
    )
  )
}

cutrun_nfcore_reference_params <- function(project) {
  ref <- cutrun_reference_resources(project)
  species <- genome_species(project)

  if (!species %in% c("human", "mouse")) {
    stop(
      "nf-core CUT&RUN currently supports the configured ",
      "CodeSpring human or mouse reference only.",
      call. = FALSE
    )
  }

  bowtie2_prefix <- normalizePath(
    ref$bowtie2_index,
    winslash = "/",
    mustWork = FALSE
  )

  bowtie2_dir <- dirname(bowtie2_prefix)
  reference_root <- dirname(bowtie2_dir)

  fasta <- if (identical(species, "human")) {
    file.path(
      reference_root,
      "GRCh38.primary_assembly.genome.fa"
    )
  } else {
    file.path(
      reference_root,
      "GRCm39.primary_assembly.genome.fa"
    )
  }

  gtf <- trimws(as.character(ref$gtf %||% "")[1])

  required_files <- c(
    fasta = fasta,
    gtf = gtf
  )

  missing_files <- required_files[
    !file.exists(required_files)
  ]

  if (length(missing_files)) {
    stop(
      "CUT&RUN nf-core reference files are missing: ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }

  if (!dir.exists(bowtie2_dir)) {
    stop(
      "CUT&RUN Bowtie2 index directory is missing: ",
      bowtie2_dir,
      call. = FALSE
    )
  }

  index_files <- Sys.glob(
    paste0(bowtie2_prefix, "*.bt2")
  )

  if (length(index_files) < 6L) {
    index_files <- Sys.glob(
      paste0(bowtie2_prefix, "*.bt2l")
    )
  }

  if (length(index_files) < 6L) {
    stop(
      "A complete Bowtie2 index was not found for: ",
      bowtie2_prefix,
      call. = FALSE
    )
  }

  list(
    fasta = normalizePath(
      fasta,
      winslash = "/",
      mustWork = TRUE
    ),
    gtf = normalizePath(
      gtf,
      winslash = "/",
      mustWork = TRUE
    ),
    bowtie2 = normalizePath(
      bowtie2_dir,
      winslash = "/",
      mustWork = TRUE
    ),
    igenomes_ignore = TRUE
  )
}

cutrun_nfcore_params <- function(
  project,
  paths,
  generated,
  normalisation_mode = "CPM",
  include_macs2 = FALSE,
  macs2_narrow_peak = TRUE,
  seacr_stringent = "stringent"
) {
  normalisation_mode <- match.arg(
    normalisation_mode,
    c("CPM", "RPKM", "BPM", "None", "Spikein")
  )

  if (identical(normalisation_mode, "Spikein")) {
    stop(
      "Spike-in normalisation is not wired into the ",
      "CodeSpring nf-core CUT&RUN backend yet. ",
      "Use CPM for the initial integration.",
      call. = FALSE
    )
  }

  seacr_stringent <- match.arg(
    seacr_stringent,
    c("stringent", "relaxed")
  )

  if (!is.data.frame(generated$samplesheet) ||
      !NROW(generated$samplesheet)) {
    stop(
      "A generated nf-core CUT&RUN samplesheet is required.",
      call. = FALSE
    )
  }

  has_controls <- any(
    nzchar(
      trimws(
        as.character(
          generated$samplesheet$control %||% ""
        )
      )
    )
  )

  reference <- cutrun_nfcore_reference_params(project)

  params <- c(
    list(
      input = paths$samplesheet_path,
      outdir = paths$output_dir
    ),
    reference,
    list(
      normalisation_mode = normalisation_mode,
      peakcaller = if (isTRUE(include_macs2)) {
        "seacr,macs2"
      } else {
        "seacr"
      },
      use_control = has_controls,
      minimum_alignment_q_score = 20L,
      seacr_norm = "non",
      seacr_stringent = seacr_stringent,
      consensus_peak_mode = "group",
      replicate_threshold = 1L
    )
  )

  if (isTRUE(include_macs2)) {
    params$macs_gsize <- if (
      identical(genome_species(project), "human")
    ) {
      2.7e9
    } else {
      1.87e9
    }

    params$macs2_narrow_peak <- isTRUE(
      macs2_narrow_peak
    )
  }

  params
}

cutrun_nfcore_nextflow_config <- function(
  runtime,
  time_limit_hours = 48L
) {
  time_limit_hours <- suppressWarnings(
    as.integer(time_limit_hours)[1]
  )

  if (is.na(time_limit_hours) ||
      time_limit_hours < 1L ||
      time_limit_hours > 48L) {
    stop(
      "CUT&RUN nf-core process time limit must be ",
      "between 1 and 48 hours.",
      call. = FALSE
    )
  }

  c(
    paste(
      "includeConfig",
      shQuote(runtime$config)
    ),
    "",
    "process {",
    paste0(
      "  resourceLimits = [ time: ",
      time_limit_hours,
      ".h ]"
    ),
    "}"
  )
}

cutrun_nfcore_launch_script <- function(
  paths,
  runtime,
  pipeline_version = CUTRUN_NFCORE_VERSION,
  nextflow_version = CUTRUN_NFCORE_NEXTFLOW_VERSION
) {
  command <- c(
    shQuote(runtime$launcher),
    "-log",
    shQuote(paths$nextflow_log),
    "-c",
    shQuote(paths$run_config),
    "run",
    shQuote(CUTRUN_NFCORE_PIPELINE),
    "-ansi-log",
    "false",
    "-r",
    shQuote(pipeline_version),
    "-profile",
    "singularity",
    "-params-file",
    shQuote(paths$params_path),
    "-work-dir",
    shQuote(paths$work_dir),
    "-with-trace",
    shQuote(paths$trace_path)
  )

  c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "",
    paste0(
      "export NXF_HOME=",
      shQuote(runtime$nxf_home)
    ),
    paste0(
      "export NXF_SINGULARITY_CACHEDIR=",
      shQuote(runtime$singularity_cache)
    ),
    paste0(
      "export NXF_VER=",
      shQuote(nextflow_version)
    ),
    "",
    paste(
      "mkdir -p",
      shQuote(paths$output_dir),
      shQuote(paths$log_dir),
      shQuote(paths$work_dir)
    ),
    paste(
      "cd",
      shQuote(paths$run_dir)
    ),
    "",
    paste(command, collapse = " ")
  )
}

cutrun_nfcore_build_bundle <- function(
  project,
  run_id = "nfcore_cutandrun",
  runtime = cutrun_nfcore_runtime_defaults(),
  normalisation_mode = "CPM",
  include_macs2 = FALSE,
  macs2_narrow_peak = TRUE,
  seacr_stringent = "stringent",
  time_limit_hours = 48L
) {
  if (!is_cutrun_project(project)) {
    stop(
      "nf-core CUT&RUN bundle generation requires ",
      "a CUT&RUN project.",
      call. = FALSE
    )
  }

  if (!requireNamespace(
    "jsonlite",
    quietly = TRUE
  )) {
    stop(
      "The jsonlite R package is required to prepare ",
      "an nf-core CUT&RUN run.",
      call. = FALSE
    )
  }

  runtime <- cutrun_nfcore_validate_runtime(runtime)
  paths <- cutrun_nfcore_submission_paths(
    project,
    run_id = run_id
  )

  if (file.exists(paths$run_dir)) {
    stop(
      "CUT&RUN nf-core run already exists: ",
      paths$run_dir,
      call. = FALSE
    )
  }

  for (directory in c(
    paths$internal_dir,
    paths$log_dir,
    paths$output_dir,
    paths$work_dir
  )) {
    if (!dir.create(
      directory,
      recursive = TRUE,
      showWarnings = FALSE
    ) && !dir.exists(directory)) {
      stop(
        "Could not create CUT&RUN nf-core directory: ",
        directory,
        call. = FALSE
      )
    }
  }

  generated <- cutrun_nfcore_samplesheet(project)

  write.csv(
    generated$samplesheet,
    paths$samplesheet_path,
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )

  write.table(
    generated$mapping,
    paths$mapping_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )

  params <- cutrun_nfcore_params(
    project = project,
    paths = paths,
    generated = generated,
    normalisation_mode = normalisation_mode,
    include_macs2 = include_macs2,
    macs2_narrow_peak = macs2_narrow_peak,
    seacr_stringent = seacr_stringent
  )

  jsonlite::write_json(
    params,
    paths$params_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  writeLines(
    cutrun_nfcore_nextflow_config(
      runtime,
      time_limit_hours = time_limit_hours
    ),
    paths$run_config,
    useBytes = TRUE
  )

  writeLines(
    cutrun_nfcore_launch_script(
      paths,
      runtime
    ),
    paths$launch_script,
    useBytes = TRUE
  )

  Sys.chmod(
    paths$launch_script,
    mode = "0700"
  )

  list(
    paths = paths,
    params = params,
    generated = generated,
    runtime = runtime
  )
}

# -------------------------------------------------------------------------
# Detached nf-core/cutandrun lifecycle
#
# Reuses the detached Nextflow controller/status machinery already proven
# for Sarek. The nf-core bundle itself remains CUT&RUN-specific.
# -------------------------------------------------------------------------

cutrun_nfcore_controller_paths <- function(
  project,
  run_id = "nfcore_cutandrun"
) {
  paths <- cutrun_nfcore_submission_paths(project, run_id)

  project_key <- as.character(
    project$id %||%
      project$name %||%
      run_id
  )[[1]]

  project_key <- gsub(
    "[^A-Za-z0-9_.-]+",
    "_",
    project_key
  )

  project_key <- sub("^_+", "", project_key)
  project_key <- sub("_+$", "", project_key)

  if (!nzchar(project_key)) {
    project_key <- run_id
  }

  # Keep the Slurm comment reasonably short.
  project_key <- substr(project_key, 1L, 80L)

  paths$manifest_id <- project_key

  paths$stdout <- file.path(
    paths$log_dir,
    "controller.out"
  )

  paths$stderr <- file.path(
    paths$log_dir,
    "controller.err"
  )

  paths$submission_record <- file.path(
    paths$internal_dir,
    "submission.tsv"
  )

  paths$runtime_status <- file.path(
    paths$internal_dir,
    "runtime_status.tsv"
  )

  paths$active_children <- file.path(
    paths$internal_dir,
    "active_children.tsv"
  )

  paths$controller_info <- file.path(
    paths$internal_dir,
    "controller.tsv"
  )

  paths$child_tag <- paste0(
    "codespring_cutrun_",
    project_key
  )

  paths
}


cutrun_nfcore_prepare_controller_bundle <- function(
  project,
  run_id = "nfcore_cutandrun",
  runtime = cutrun_nfcore_runtime_defaults(),
  normalisation_mode = "CPM",
  include_macs2 = FALSE,
  macs2_narrow_peak = TRUE,
  seacr_stringent = "stringent",
  time_limit_hours = 48L
) {
  runtime <- cutrun_nfcore_validate_runtime(runtime)

  bundle <- cutrun_nfcore_build_bundle(
    project = project,
    run_id = run_id,
    runtime = runtime,
    normalisation_mode = normalisation_mode,
    include_macs2 = include_macs2,
    macs2_narrow_peak = macs2_narrow_peak,
    seacr_stringent = seacr_stringent,
    time_limit_hours = time_limit_hours
  )

  paths <- cutrun_nfcore_controller_paths(
    project,
    run_id
  )

  # Every nf-core/cutandrun process must execute through Slurm.
  # The lightweight Nextflow controller itself remains detached on the
  # application host, matching the established Sarek lifecycle.
  #
  # bam01 currently cannot resolve this Unix account correctly inside
  # Singularity. Allow an environment override, but keep bam01 as the
  # temporary CSHL safety default until the node is repaired.
  excluded_nodes <- trimws(
    Sys.getenv(
      "CSL_SLURM_EXCLUDE_NODES",
      unset = "bam01"
    )
  )

  child_cluster_options <- paste0(
    "--comment=",
    paths$child_tag
  )

  if (nzchar(excluded_nodes)) {
    child_cluster_options <- paste(
      child_cluster_options,
      paste0(
        "--exclude=",
        excluded_nodes
      )
    )
  }

  config_lines <- c(
    paste0(
      "includeConfig '",
      runtime$config,
      "'"
    ),
    "",
    "process {",
    "  executor = 'slurm'",
    paste0(
      "  clusterOptions = '",
      child_cluster_options,
      "'"
    ),
    paste0(
      "  resourceLimits = [ time: ",
      as.integer(time_limit_hours),
      ".h ]"
    ),
    "}"
  )

  writeLines(
    config_lines,
    paths$run_config
  )

  # Use the already-proven detached-controller wrapper. It provides:
  # - controller PID tracking
  # - runtime_status.tsv
  # - active_children.tsv
  # - initial vs -resume operation
  launch_lines <- sarek_submission_launch_script(
    paths = paths,
    runtime = runtime,
    pipeline = CUTRUN_NFCORE_PIPELINE,
    pipeline_version = CUTRUN_NFCORE_VERSION,
    nextflow_version = CUTRUN_NFCORE_NEXTFLOW_VERSION
  )

  writeLines(
    launch_lines,
    paths$launch_script
  )

  Sys.chmod(
    paths$launch_script,
    mode = "0700"
  )

  bundle$paths <- paths
  bundle
}


cutrun_nfcore_submit_run <- function(
  project,
  run_id = "nfcore_cutandrun",
  runtime = cutrun_nfcore_runtime_defaults(),
  normalisation_mode = "CPM",
  include_macs2 = FALSE,
  macs2_narrow_peak = TRUE,
  seacr_stringent = "stringent",
  time_limit_hours = 48L,
  starter = NULL
) {
  if (!is_nfcore_cutrun_project(project)) {
    stop(
      "nf-core CUT&RUN submission requires an nf-core CUT&RUN project."
    )
  }

  if (!isTRUE(project$paired_end)) {
    stop(
      "nf-core/cutandrun 3.2.2 supports paired-end CUT&RUN data only."
    )
  }

  bundle <- cutrun_nfcore_prepare_controller_bundle(
    project = project,
    run_id = run_id,
    runtime = runtime,
    normalisation_mode = normalisation_mode,
    include_macs2 = include_macs2,
    macs2_narrow_peak = macs2_narrow_peak,
    seacr_stringent = seacr_stringent,
    time_limit_hours = time_limit_hours
  )

  paths <- bundle$paths

  controller <- sarek_start_detached_controller(
    paths = paths,
    mode = "initial",
    starter = starter
  )

  submitted_at <- format(
    Sys.time(),
    "%Y-%m-%dT%H:%M:%SZ",
    tz = "UTC"
  )

  sarek_submission_write_values(
    paths$submission_record,
    c(
      status = "submitted",
      controller_pid = controller$pid,
      controller_host = controller$host,
      controller_mode = "detached",
      attempt = "1",
      mode = "initial",
      child_tag = paths$child_tag,
      submitted_at = submitted_at,
      run_dir = paths$run_dir,
      output_dir = paths$output_dir,
      work_dir = paths$work_dir,
      pipeline = CUTRUN_NFCORE_PIPELINE,
      pipeline_version = CUTRUN_NFCORE_VERSION,
      nextflow_version = CUTRUN_NFCORE_NEXTFLOW_VERSION
    )
  )

  list(
    status = "submitted",
    controller_pid = controller$pid,
    controller_host = controller$host,
    child_tag = paths$child_tag,
    run_dir = paths$run_dir,
    output_dir = paths$output_dir,
    work_dir = paths$work_dir
  )
}


cutrun_nfcore_run_status <- function(
  project,
  run_id = "nfcore_cutandrun",
  runner = NULL,
  squeue = "squeue",
  sacct = "sacct"
) {
  paths <- cutrun_nfcore_controller_paths(
    project,
    run_id
  )

  # A run has not actually started until CodeSpring has written its
  # submission record. The project/run directory existing by itself is
  # not evidence of a submitted Nextflow run.
  has_submission <- (
    file.exists(paths$submission_record) &&
    !dir.exists(paths$submission_record)
  )

  if (!has_submission) {
    return(list(
      state = "NOT_STARTED",
      has_submission = FALSE,
      controller_pid = "",
      controller_alive = FALSE,
      active_children = 0L,
      child_jobs = data.frame(),
      run_dir = paths$run_dir,
      output_dir = paths$output_dir
    ))
  }

  values <- sarek_read_key_value_file(
    paths$submission_record
  )

  run <- list(
    status = sarek_text(
      values["status"],
      "submitted"
    ),
    run_dir = paths$run_dir,
    output_dir = paths$output_dir,
    work_dir = paths$work_dir,
    runtime_status = paths$runtime_status,
    controller_pid = sarek_text(
      values["controller_pid"]
    ),
    controller_host = sarek_text(
      values["controller_host"]
    ),
    child_tag = sarek_text(
      values["child_tag"],
      paths$child_tag
    )
  )

  status <- sarek_run_status(
    run,
    runner = runner,
    squeue = squeue,
    sacct = sacct
  )

  status$has_submission <- TRUE
  status$run_dir <- paths$run_dir
  status$output_dir <- paths$output_dir

  active_children <- suppressWarnings(
    as.integer(status$active_children %||% 0L)
  )

  if (
    is.na(active_children) ||
    active_children < 0L
  ) {
    active_children <- 0L
  }

  status$active_children <- active_children

  # Never report a CUT&RUN run as terminal while its detached
  # controller or tagged Slurm children are still alive.
  if (
    isTRUE(status$controller_alive) ||
    active_children > 0L
  ) {
    status$state <- "RUNNING"
    status$source <- "live"
  }

  status
}

cutrun_nfcore_archive_trace <- function(paths) {
  trace_path <- paths$trace_path %||%
    file.path(
      paths$run_dir,
      ".codespring",
      "logs",
      "trace.tsv"
    )

  if (
    !nzchar(trace_path) ||
    !file.exists(trace_path) ||
    dir.exists(trace_path)
  ) {
    return("")
  }

  stamp <- format(
    Sys.time(),
    "%Y%m%d_%H%M%S"
  )

  archive_path <- file.path(
    dirname(trace_path),
    paste0(
      "trace.before_resume_",
      stamp,
      ".tsv"
    )
  )

  counter <- 1L

  while (file.exists(archive_path)) {
    archive_path <- file.path(
      dirname(trace_path),
      paste0(
        "trace.before_resume_",
        stamp,
        "_",
        counter,
        ".tsv"
      )
    )

    counter <- counter + 1L
  }

  if (!file.rename(
    trace_path,
    archive_path
  )) {
    stop(
      "Could not archive the existing Nextflow trace before resume: ",
      trace_path
    )
  }

  archive_path
}


cutrun_nfcore_resume_run <- function(
  project,
  run_id = "nfcore_cutandrun",
  runtime = cutrun_nfcore_runtime_defaults(),
  starter = NULL
) {
  runtime <- cutrun_nfcore_validate_runtime(
    runtime
  )

  paths <- cutrun_nfcore_controller_paths(
    project,
    run_id
  )

  if (!dir.exists(paths$run_dir)) {
    stop(
      "No nf-core CUT&RUN run exists to resume."
    )
  }

  required <- c(
    paths$launch_script,
    paths$params_path,
    paths$run_config,
    paths$work_dir
  )

  missing <- required[
    !file.exists(required) &
      !dir.exists(required)
  ]

  if (length(missing)) {
    stop(
      "Cannot resume nf-core CUT&RUN because required run state is missing: ",
      paste(missing, collapse = ", ")
    )
  }

  launch_lines <- readLines(
    paths$launch_script,
    warn = FALSE
  )

  if (
    !any(
      grepl(
        "CSL_SAREK_RESUME",
        launch_lines,
        fixed = TRUE
      )
    )
  ) {
    stop(
      "This CUT&RUN run was created before detached resume support was added."
    )
  }

  current <- cutrun_nfcore_run_status(
    project,
    run_id
  )

  active_children <- suppressWarnings(
    as.integer(current$active_children %||% 0L)
  )

  if (
    is.na(active_children) ||
    active_children < 0L
  ) {
    active_children <- 0L
  }

  if (
    isTRUE(current$controller_alive) ||
    active_children > 0L ||
    current$state %in% c(
      "RUNNING",
      "PENDING",
      "CONFIGURING",
      "COMPLETING"
    )
  ) {
    stop(
      "Cannot resume nf-core CUT&RUN while the controller or tagged Slurm tasks are still active."
    )
  }

  # Nextflow refuses to overwrite an existing -with-trace file.
  # Preserve the previous attempt before starting -resume.
  cutrun_nfcore_archive_trace(paths)

  values <- sarek_read_key_value_file(
    paths$submission_record
  )

  recorded_tag <- sarek_text(
    values["child_tag"]
  )

  if (nzchar(recorded_tag)) {
    paths$child_tag <- recorded_tag
  }

  controller <- sarek_start_detached_controller(
    paths = paths,
    mode = "resume",
    starter = starter
  )

  previous_attempt <- suppressWarnings(
    as.integer(
      sarek_text(
        values["attempt"],
        "1"
      )
    )
  )

  if (
    is.na(previous_attempt) ||
    previous_attempt < 1L
  ) {
    previous_attempt <- 1L
  }

  submitted_at <- format(
    Sys.time(),
    "%Y-%m-%dT%H:%M:%SZ",
    tz = "UTC"
  )

  sarek_submission_write_values(
    paths$submission_record,
    c(
      status = "submitted",
      controller_pid = controller$pid,
      controller_host = controller$host,
      controller_mode = "detached",
      attempt = as.character(
        previous_attempt + 1L
      ),
      mode = "resume",
      child_tag = paths$child_tag,
      submitted_at = submitted_at,
      run_dir = paths$run_dir,
      output_dir = paths$output_dir,
      work_dir = paths$work_dir,
      pipeline = CUTRUN_NFCORE_PIPELINE,
      pipeline_version = CUTRUN_NFCORE_VERSION,
      nextflow_version = CUTRUN_NFCORE_NEXTFLOW_VERSION
    )
  )

  list(
    status = "submitted",
    controller_pid = controller$pid,
    controller_host = controller$host,
    child_tag = paths$child_tag,
    run_dir = paths$run_dir,
    output_dir = paths$output_dir,
    work_dir = paths$work_dir
  )
}
