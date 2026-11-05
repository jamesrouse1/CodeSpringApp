# Result discovery and per-user preferences for the Sarek Shiny module.

SAREK_RESULTS_HELPERS_VERSION <- "0.1.0"

sarek_results_run_value <- function(run, name, default = "") {
  if (is.null(run) || !name %in% names(run)) return(default)
  sarek_text(run[[name]], default)
}

sarek_results_manifest_path <- function(run) {
  run_dir <- sarek_results_run_value(run, "run_dir")
  if (!nzchar(run_dir)) {
    output_dir <- sarek_results_run_value(run, "output_dir")
    if (nzchar(output_dir)) run_dir <- dirname(output_dir)
  }
  if (!nzchar(run_dir)) return("")
  file.path(run_dir, ".codespring", "manifest.json")
}

sarek_results_read_manifest <- function(run) {
  path <- sarek_results_manifest_path(run)
  if (
    !nzchar(path) ||
    !file.exists(path) ||
    !requireNamespace("jsonlite", quietly = TRUE)
  ) {
    return(NULL)
  }

  tryCatch(
    jsonlite::fromJSON(path, simplifyVector = FALSE),
    error = function(error) NULL
  )
}

sarek_results_manifest_samples <- function(manifest) {
  if (is.null(manifest) || !is.list(manifest) || !length(manifest$patients)) {
    return(character(0))
  }

  samples <- unlist(lapply(manifest$patients, function(patient) {
    if (is.null(patient$samples) || !length(patient$samples)) {
      return(character(0))
    }

    vapply(patient$samples, function(sample) {
      sarek_text(sample$sample_id)
    }, character(1))
  }), use.names = FALSE)

  sort(unique(samples[nzchar(samples)]))
}

sarek_results_canonical_text <- function(value) {
  value <- tolower(as.character(value))
  value <- gsub("[^a-z0-9]+", "_", value)
  gsub("^_+|_+$", "", value)
}

sarek_results_samples_for_path <- function(path, sample_ids) {
  if (!length(sample_ids)) return(character(0))

  path_key <- paste0("_", sarek_results_canonical_text(path), "_")
  sample_ids[vapply(sample_ids, function(sample_id) {
    sample_key <- sarek_results_canonical_text(sample_id)
    nzchar(sample_key) &&
      grepl(paste0("_", sample_key, "_"), path_key, fixed = TRUE)
  }, logical(1))]
}

sarek_results_category <- function(path, file_type = "") {
  key <- tolower(gsub("\\\\", "/", path))

  if (grepl("multiqc", key, fixed = TRUE)) {
    return("MultiQC")
  }

  if (grepl(
    "pipeline_info|execution_report|execution_timeline|pipeline_dag|trace[.]txt",
    key,
    perl = TRUE
  )) {
    return("Pipeline execution")
  }

  if (grepl("/(snpeff|vep|snpsift|bcfann)/", key, perl = TRUE)) {
    return("Variant annotation")
  }

  if (grepl("/(ascat|cnvkit|controlfreec)/", key, perl = TRUE)) {
    return("Copy number")
  }

  if (grepl("/(manta|tiddit|indexcov)/", key, perl = TRUE)) {
    return("Structural variants")
  }

  if (grepl(
    "/(bcftools|deepvariant|freebayes|haplotypecaller|lofreq|muse|mutect2|strelka|sentieon)/",
    key,
    perl = TRUE
  )) {
    return("Small variants")
  }

  if (grepl("/(preprocessing|mapped|markduplicates|recalibrated)/", key, perl = TRUE)) {
    return("Alignments and preprocessing")
  }

  if (grepl("/(reports|qc|fastqc|mosdepth|samtools|vcftools)/", key, perl = TRUE)) {
    return("Quality control")
  }

  if (identical(file_type, "Report")) return("Other reports")
  if (identical(file_type, "Variants")) return("Variant files")
  if (identical(file_type, "Alignment")) return("Alignment files")
  if (identical(file_type, "Regions")) return("Regions")
  if (identical(file_type, "Table / metadata")) return("Tables and metadata")

  "Other files"
}

sarek_results_preview_kind <- function(path) {
  extension <- tolower(tools::file_ext(path))

  if (extension %in% c("html", "htm")) return("html")
  if (extension == "pdf") return("pdf")
  if (extension %in% c("png", "jpg", "jpeg", "webp", "svg")) return("image")
  if (extension %in% c("csv", "tsv")) return("table")
  if (extension %in% c("txt", "json", "yaml", "yml")) return("text")

  "none"
}

sarek_results_artifact_catalog <- function(run, max_files = 5000L) {
  output_dir <- sarek_results_run_value(run, "output_dir")
  files <- sarek_result_file_catalog(output_dir, max_files = max_files)
  truncated <- isTRUE(attr(files, "truncated"))

  files$category <- character(NROW(files))
  files$sample <- character(NROW(files))
  files$sample_ids <- character(NROW(files))
  files$preview_kind <- character(NROW(files))

  if (!NROW(files)) {
    attr(files, "truncated") <- truncated
    attr(files, "manifest_samples") <- character(0)
    return(files)
  }

  manifest <- sarek_results_read_manifest(run)
  sample_ids <- sarek_results_manifest_samples(manifest)

  matched_samples <- lapply(files$path, function(path) {
    sarek_results_samples_for_path(path, sample_ids)
  })

  files$category <- vapply(seq_len(NROW(files)), function(index) {
    sarek_results_category(files$path[[index]], files$type[[index]])
  }, character(1))

  files$sample_ids <- vapply(
    matched_samples,
    paste,
    collapse = ";",
    character(1)
  )

  files$sample <- vapply(matched_samples, function(samples) {
    if (length(samples)) paste(samples, collapse = " + ") else "Run-wide"
  }, character(1))

  files$preview_kind <- vapply(
    files$path,
    sarek_results_preview_kind,
    character(1)
  )

  category_order <- match(
    files$category,
    c(
      "MultiQC",
      "Pipeline execution",
      "Variant annotation",
      "Small variants",
      "Structural variants",
      "Copy number",
      "Quality control",
      "Other reports",
      "Alignments and preprocessing",
      "Variant files",
      "Alignment files",
      "Regions",
      "Tables and metadata",
      "Other files"
    )
  )

  files <- files[
    order(
      files$preview_kind == "none",
      category_order,
      files$sample,
      files$folder,
      files$file
    ),
    ,
    drop = FALSE
  ]

  rownames(files) <- NULL
  attr(files, "truncated") <- truncated
  attr(files, "manifest_samples") <- sample_ids
  files
}

sarek_results_filter_sample <- function(catalog, selected) {
  selected <- sarek_text(selected, "__all__")

  if (!NROW(catalog) || identical(selected, "__all__")) {
    return(catalog)
  }

  if (identical(selected, "__run__")) {
    return(catalog[!nzchar(catalog$sample_ids), , drop = FALSE])
  }

  needle <- paste0(";", selected, ";")
  keep <- grepl(
    needle,
    paste0(";", catalog$sample_ids, ";"),
    fixed = TRUE
  )

  catalog[keep, , drop = FALSE]
}

sarek_results_read_default_root <- function(preference_path, fallback) {
  fallback <- normalizePath(
    path.expand(fallback),
    winslash = "/",
    mustWork = FALSE
  )

  if (
    is.null(preference_path) ||
    !length(preference_path) ||
    !nzchar(preference_path) ||
    !file.exists(preference_path)
  ) {
    return(fallback)
  }

  saved <- tryCatch(
    trimws(readLines(preference_path, warn = FALSE, n = 1L)),
    error = function(error) ""
  )

  if (!length(saved) || !nzchar(saved[[1]]) || !sarek_is_absolute_path(saved[[1]])) {
    return(fallback)
  }

  normalizePath(
    path.expand(saved[[1]]),
    winslash = "/",
    mustWork = FALSE
  )
}

sarek_results_save_default_root <- function(
    root,
    preference_path,
    allowed_roots = NULL) {
  root <- trimws(as.character(root))

  if (!nzchar(root) || !sarek_is_absolute_path(root)) {
    stop("The default Sarek results location must be an absolute server path.")
  }

  root <- normalizePath(path.expand(root), winslash = "/", mustWork = FALSE)
  sarek_assert_allowed_paths(root, allowed_roots, must_exist = FALSE)

  if (!dir.exists(root)) {
    if (!dir.create(root, recursive = TRUE, showWarnings = FALSE)) {
      stop("The selected results directory could not be created: ", root)
    }
  }

  if (file.access(root, mode = 3) != 0) {
    stop("The selected results directory is not readable and writable: ", root)
  }

  if (is.null(preference_path) || !nzchar(preference_path)) {
    stop("The Sarek results preference file is not configured.")
  }

  dir.create(dirname(preference_path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = ".sarek_results_root_",
    tmpdir = dirname(preference_path)
  )
  on.exit(unlink(temporary, force = TRUE), add = TRUE)

  writeLines(root, temporary, useBytes = TRUE)

  if (!file.rename(temporary, preference_path)) {
    if (!file.copy(temporary, preference_path, overwrite = TRUE)) {
      stop("The default Sarek results preference could not be saved.")
    }
  }

  root
}

sarek_results_preview_table <- function(path, max_rows = 500L) {
  extension <- tolower(tools::file_ext(path))

  tryCatch({
    if (extension == "csv") {
      return(utils::read.csv(
        path,
        nrows = max_rows,
        check.names = FALSE,
        stringsAsFactors = FALSE
      ))
    }

    if (extension == "tsv") {
      return(utils::read.delim(
        path,
        nrows = max_rows,
        check.names = FALSE,
        stringsAsFactors = FALSE
      ))
    }

    lines <- readLines(path, warn = FALSE, n = max_rows)
    data.frame(
      Line = seq_along(lines),
      Content = lines,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }, error = function(error) {
    data.frame(
      Message = paste("This file could not be previewed:", conditionMessage(error)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
}
