# Read-only presentation helpers for Sarek outputs. The caller must select paths
# from the run's catalog; no browser-provided path is opened directly.
sarek_explorer_file <- function(path, output_dir) {
  if (!length(path) || !length(output_dir) || is.na(path[[1]]) || is.na(output_dir[[1]])) return("")
  root <- normalizePath(as.character(output_dir[[1]]), winslash = "/", mustWork = FALSE)
  file <- normalizePath(as.character(path[[1]]), winslash = "/", mustWork = FALSE)
  if (!dir.exists(root) || !file.exists(file) || dir.exists(file) ||
      !startsWith(file, paste0(root, "/")) || file.access(file, 4) != 0) return("")
  file
}

sarek_explorer_samples <- function(catalog) {
  if (!NROW(catalog)) return(character())
  values <- as.character(catalog$sample_ids)
  values <- trimws(unlist(strsplit(values[!is.na(values)], "[,;|]", perl = TRUE)))
  sort(unique(values[nzchar(values) & !tolower(values) %in% c("all", "shared", "unknown", "na")]))
}

sarek_explorer_info <- function(value, key) {
  hit <- strsplit(value, ";", fixed = TRUE)[[1]]
  hit <- hit[startsWith(hit, paste0(key, "="))]
  if (length(hit)) sub("^[^=]+=", "", hit[[1]]) else ""
}

sarek_explorer_format <- function(format, sample, field) {
  keys <- strsplit(format, ":", fixed = TRUE)[[1]]
  values <- strsplit(sample, ":", fixed = TRUE)[[1]]
  index <- match(field, keys)
  if (is.na(index) || index > length(values)) "" else values[[index]]
}

sarek_explorer_caller <- function(folder, filename) {
  parts <- strsplit(as.character(folder), "/", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts) & parts != "."]
  if (length(parts) >= 4L && identical(parts[[1]], "annotation") &&
      identical(parts[[2]], "clinvar")) {
    call_index <- match("variant_calling", parts)
    if (!is.na(call_index) && call_index < length(parts)) return(parts[[call_index + 1L]])
  }
  if (length(parts) >= 2L && parts[[1]] %in% c("variant_calling", "annotation")) return(parts[[2]])
  if (length(parts)) return(parts[[length(parts)]])
  sub("[._].*$", "", filename)
}

sarek_explorer_empty_variants <- function() {
  data.frame(
    Sample = character(), Caller = character(), CHROM = character(),
    POS = character(), REF = character(), ALT = character(), FILTER = character(),
    Genotype = character(), Depth = character(), VAF = character(),
    Gene = character(), Consequence = character(), Impact = character(),
    ClinVar = character(), Review = character(), Condition = character(),
    Source = character(), stringsAsFactors = FALSE, check.names = FALSE
  )
}

# A bounded preview: never silently present its counts as full-run totals.
sarek_explorer_variants <- function(catalog, sample = "All", output_dir,
                                    max_files = 8L, max_records = 1500L, strict = FALSE) {
  result <- sarek_explorer_empty_variants()
  files <- catalog[grepl("\\.vcf(\\.gz)?$", catalog$file, ignore.case = TRUE) &
                     !grepl("\\.g\\.vcf", catalog$file, ignore.case = TRUE), , drop = FALSE]
  if (!identical(sample, "All")) files <- sarek_results_filter_sample(files, sample)
  total_files <- NROW(files)
  files <- utils::head(files, max_files)
  limited <- total_files > NROW(files)
  for (i in seq_len(NROW(files))) {
    path <- sarek_explorer_file(files$path[[i]], output_dir)
    if (!nzchar(path)) {
      if (strict) stop("The selected VCF is unavailable or outside the run results directory.")
      next
    }
    con <- tryCatch(if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt"),
                    error = function(error) NULL)
    if (is.null(con)) {
      if (strict) stop("The selected VCF could not be opened.")
      next
    }
    read_one <- function() {
      header <- character()
      column_line <- ""
      for (j in seq_len(10000L)) {
        line <- readLines(con, n = 1L, warn = FALSE)
        if (!length(line)) break
        if (startsWith(line, "#CHROM")) { column_line <- line; break }
        if (startsWith(line, "##INFO=<ID=CSQ,")) header <- c(header, line)
      }
      if (!nzchar(column_line)) {
        if (strict) stop("No #CHROM header found within the first 10,000 lines.")
        return(NULL)
      }
      columns <- strsplit(column_line, "\t", fixed = TRUE)[[1]]
      samples <- if (length(columns) > 9L) columns[10:length(columns)] else character()
      selected <- if (identical(sample, "All")) {
        if (length(samples) == 1L) samples[[1]] else ""
      } else if (sample %in% samples) sample else ""
      sample_index <- match(selected, columns)
      csq_fields <- if (length(header)) {
        csq <- sub('.*Format: ', '', header[[1]])
        csq <- sub('[">].*$', '', csq)
        strsplit(csq, "|", fixed = TRUE)[[1]]
      } else character()
      records <- character()
      repeat {
        chunk <- readLines(con, n = min(500L, max_records - length(records)), warn = FALSE)
        if (!length(chunk)) break
        records <- c(records, chunk[!startsWith(chunk, "#")])
        if (length(records) >= max_records) { limited <<- TRUE; break }
      }
      if (!length(records)) return(NULL)
      parsed <- lapply(records, function(line) {
        fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
        if (length(fields) < 8L) {
          if (strict) stop("A VCF record has fewer than eight columns.")
          return(NULL)
        }
        info <- fields[[8]]
        ann <- sarek_explorer_info(info, "ANN")
        csq <- sarek_explorer_info(info, "CSQ")
        gene <- consequence <- impact <- ""
        if (nzchar(ann)) {
          parts <- strsplit(strsplit(ann, ",", fixed = TRUE)[[1]][[1]], "|", fixed = TRUE)[[1]]
          if (length(parts) >= 4L) {
            consequence <- parts[[2]]; impact <- parts[[3]]; gene <- parts[[4]]
          }
        } else if (nzchar(csq) && length(csq_fields)) {
          parts <- strsplit(strsplit(csq, ",", fixed = TRUE)[[1]][[1]], "|", fixed = TRUE)[[1]]
          get <- function(name) {
            pos <- match(name, csq_fields)
            if (is.na(pos) || pos > length(parts)) "" else parts[[pos]]
          }
          consequence <- get("Consequence"); impact <- get("IMPACT"); gene <- get("SYMBOL")
        }
        format <- if (length(fields) >= 9L) fields[[9]] else ""
        genotype <- if (!is.na(sample_index) && length(fields) >= sample_index) fields[[sample_index]] else ""
        clnsig <- sarek_explorer_info(info, "CLNSIG")
        data.frame(
          Sample = if (nzchar(selected)) selected else as.character(files$sample[[i]]),
          Caller = sarek_explorer_caller(files$folder[[i]], files$file[[i]]),
          CHROM = fields[[1]], POS = fields[[2]], REF = fields[[4]], ALT = fields[[5]],
          FILTER = fields[[7]],
          Genotype = sarek_explorer_format(format, genotype, "GT"),
          Depth = sarek_explorer_format(format, genotype, "DP"),
          VAF = sarek_explorer_format(format, genotype, "AF"),
          Gene = gene, Consequence = consequence, Impact = impact,
          ClinVar = clnsig,
          Review = sarek_explorer_info(info, "CLNREVSTAT"),
          Condition = sarek_explorer_info(info, "CLNDN"),
          Source = files$file[[i]], stringsAsFactors = FALSE, check.names = FALSE
        )
      })
      parsed <- Filter(Negate(is.null), parsed)
      if (length(parsed)) do.call(rbind, parsed) else NULL
    }
    rows <- tryCatch(
      withCallingHandlers(read_one(), warning = function(w) {
        if (strict) stop(conditionMessage(w), call. = FALSE)
      }),
      error = function(error) {
        if (strict) stop(conditionMessage(error), call. = FALSE)
        NULL
      },
      finally = close(con)
    )
    if (!is.null(rows)) result <- rbind(result, rows)
  }
  attr(result, "limited") <- limited
  attr(result, "files_scanned") <- NROW(files)
  attr(result, "files_available") <- total_files
  result
}

sarek_explorer_plots <- function(catalog) {
  catalog[grepl("\\.(png|jpe?g|webp|svg|pdf)$", catalog$file, ignore.case = TRUE), , drop = FALSE]
}

# Include the same Quality control outputs shown in Files. Run-wide reports
# remain available while viewing an individual sample.
sarek_explorer_qc <- function(catalog, sample = "All") {
  if (!NROW(catalog)) return(catalog)
  names <- basename(as.character(catalog$file))
  folders <- as.character(catalog$folder)
  categories <- as.character(catalog$category)
  shared <- grepl("multiqc", paste(names, folders), ignore.case = TRUE)
  qc <- grepl("quality[[:space:]_-]*control|(^|[^a-z])qc([^a-z]|$)", categories,
              ignore.case = TRUE, perl = TRUE) |
        grepl("multiqc|fastqc|fastp|qualimap", paste(names, folders), ignore.case = TRUE)
  qc[is.na(qc)] <- FALSE
  if (!identical(sample, "All")) {
    matching <- sarek_results_filter_sample(catalog, sample)
    qc <- qc & (catalog$path %in% matching$path |
      shared | is.na(catalog$sample) | !nzchar(as.character(catalog$sample)) |
      tolower(as.character(catalog$sample)) %in% c("all", "shared"))
    qc[is.na(qc)] <- FALSE
  }
  result <- catalog[qc, , drop = FALSE]
  result[order(!grepl("multiqc.*\\.html?$", basename(as.character(result$file)),
                        ignore.case = TRUE),
               !grepl("\\.(html?|pdf|png|jpe?g|webp|svg)$", result$file,
                      ignore.case = TRUE),
               basename(as.character(result$file))), , drop = FALSE]
}

sarek_explorer_preview <- function(path, output_dir, max_bytes = 20 * 1024^2) {
  safe <- sarek_explorer_file(path, output_dir)
  if (!nzchar(safe)) return(shiny::div(class = "empty-box", "This output is unavailable."))
  bytes <- file.info(safe)$size
  if (is.na(bytes) || bytes > max_bytes) {
    return(shiny::div(class = "empty-box", "Preview is too large; use Download selected output."))
  }
  ext <- tolower(tools::file_ext(safe))
  if (ext %in% c("htm", "html")) {
    document <- paste(readLines(safe, warn = FALSE), collapse = "\n")
    document <- paste0(
      '<meta http-equiv="Content-Security-Policy" content="default-src \'none\'; ',
      'script-src \'unsafe-inline\'; style-src \'unsafe-inline\'; ',
      'img-src data: blob:; font-src data:; connect-src \'none\'">',
      document
    )
    return(shiny::tags$iframe(
      srcdoc = document, sandbox = "allow-scripts",
      style = "width:100%;height:720px;border:1px solid #d7e0ea;background:white;"
    ))
  }
  mime <- switch(ext, png = "image/png", jpg = "image/jpeg", jpeg = "image/jpeg",
                 webp = "image/webp", svg = "image/svg+xml", pdf = "application/pdf", "")
  if (!nzchar(mime) || !requireNamespace("base64enc", quietly = TRUE)) {
    return(shiny::div(class = "empty-box", "Preview unavailable; use Download selected output."))
  }
  src <- paste0("data:", mime, ";base64,", base64enc::base64encode(safe))
  if (ext == "pdf") {
    shiny::tags$object(
      data = src, type = "application/pdf",
      style = "width:100%;height:720px;border:1px solid #d7e0ea;",
      shiny::p("PDF preview unavailable in this browser. Use Download selected output.")
    )
  } else {
    shiny::tags$img(src = src, alt = basename(safe), style = "max-width:100%;max-height:720px;object-fit:contain;")
  }
}

sarek_explorer_report_table <- function(rows, limit = 50L) {
  if (!NROW(rows)) return(NULL)
  rows <- utils::head(rows, limit)
  shiny::div(class = "sarek-report-table", shiny::tags$table(
    class = "table table-striped",
    shiny::tags$thead(shiny::tags$tr(lapply(names(rows), shiny::tags$th))),
    shiny::tags$tbody(lapply(seq_len(NROW(rows)), function(i)
      shiny::tags$tr(lapply(as.character(unlist(rows[i, , drop = FALSE], use.names = FALSE)),
                            shiny::tags$td))))
  ))
}

sarek_explorer_multiqc_stats <- function(catalog, output_dir) {
  candidates <- catalog[grepl("^multiqc_general_stats\\.(txt|tsv)$",
                              basename(as.character(catalog$file)), ignore.case = TRUE), , drop = FALSE]
  for (path in candidates$path) {
    safe <- sarek_explorer_file(path, output_dir)
    if (!nzchar(safe) || file.info(safe)$size > 3 * 1024^2) next
    stats <- tryCatch(utils::read.delim(safe, nrows = 30L, check.names = FALSE,
                                        stringsAsFactors = FALSE, quote = "", comment.char = ""),
                      error = function(error) NULL)
    if (!is.null(stats) && NROW(stats) && NCOL(stats) > 1L) {
      return(stats[, seq_len(min(NCOL(stats), 9L)), drop = FALSE])
    }
  }
  NULL
}

sarek_explorer_fastp_stats <- function(catalog, output_dir) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
  candidates <- catalog[grepl("\\.json$", catalog$file, ignore.case = TRUE) &
    grepl("fastp", paste(catalog$file, catalog$folder), ignore.case = TRUE), , drop = FALSE]
  rows <- lapply(seq_len(min(NROW(candidates), 30L)), function(i) {
    safe <- sarek_explorer_file(candidates$path[[i]], output_dir)
    if (!nzchar(safe) || is.na(file.info(safe)$size) ||
        file.info(safe)$size > 3 * 1024^2) return(NULL)
    data <- tryCatch(jsonlite::fromJSON(safe), error = function(error) NULL)
    if (!is.list(data) || !is.list(data$summary)) return(NULL)
    before <- data$summary$before_filtering
    after <- data$summary$after_filtering
    if (!is.list(before) || !is.list(after)) return(NULL)
    value <- function(x) {
      if (!length(x)) return(NA_real_)
      number <- suppressWarnings(as.numeric(x[[1]]))
      if (length(number) == 1L && is.finite(number)) number else NA_real_
    }
    data.frame(Sample = as.character(candidates$sample[[i]]),
      `Reads before filtering` = value(before$total_reads),
      `Reads after filtering` = value(after$total_reads),
      `Q30 after filtering (%)` = round(100 * value(after$q30_rate), 1),
      check.names = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows)) do.call(rbind, rows) else NULL
}

sarek_explorer_report_ui <- function(run, catalog, variants, manifest = NULL, oncoplot = NULL) {
  mode <- if (!is.null(manifest)) sarek_text(manifest$analysis$mode, "Unknown") else "Unknown"
  assay <- if (!is.null(manifest)) sarek_text(manifest$assay$type, "Unknown") else "Unknown"
  reference <- if (!is.null(manifest)) sarek_text(manifest$reference$assembly, "Unknown") else "Unknown"
  files <- NROW(catalog)
  qc <- sarek_explorer_qc(catalog)
  plots <- sarek_explorer_plots(catalog)
  stats <- sarek_explorer_multiqc_stats(catalog, sarek_text(run$output_dir))
  fastp_stats <- sarek_explorer_fastp_stats(qc, sarek_text(run$output_dir))
  qc_names <- basename(as.character(qc$file))
  priority <- variants[
    toupper(variants$FILTER) == "PASS" &
      (toupper(variants$Impact) == "HIGH" |
         grepl("(^|[|,])(?:likely_)?pathogenic([|,]|$)", variants$ClinVar,
               ignore.case = TRUE, perl = TRUE)), , drop = FALSE
  ]
  inventory <- if (files) {
    counts <- sort(table(catalog$category), decreasing = TRUE)
    data.frame("Output type" = names(counts), "Files" = as.integer(counts),
               check.names = FALSE)
  } else NULL
  sample_counts <- if (files) {
    values <- as.character(catalog$sample)
    values[is.na(values) | !nzchar(values)] <- "Shared / unassigned"
    counts <- sort(table(values), decreasing = TRUE)
    data.frame("Sample" = names(counts), "Output files" = as.integer(counts),
               check.names = FALSE)
  } else NULL
  variant_counts <- if (NROW(variants)) {
    counts <- as.data.frame(table(variants$Sample, variants$Caller), stringsAsFactors = FALSE)
    counts <- counts[counts$Freq > 0L, , drop = FALSE]
    names(counts) <- c("Sample", "Caller", "Inspected records")
    counts
  } else NULL
  plot_images <- plots[grepl("\\.(png|jpe?g|webp)$", plots$file, ignore.case = TRUE), , drop = FALSE]
  plot_images <- utils::head(plot_images, 3L)
  shiny::div(
    class = "sarek-generated-report",
    shiny::h2("Sarek analysis report"),
    shiny::p(shiny::strong("Run: "), sarek_text(run$run_id),
             " · ", mode, " · ", assay, " · ", reference),
    shiny::p(shiny::strong("Generated: "), format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
    shiny::h3("Somatic alterations by tumor"),
    if (is.null(oncoplot)) {
      if (!identical(mode, "matched_tumor_normal")) {
        shiny::p("The paired somatic oncoplot applies to matched tumor/normal runs.")
      } else if (!any(grepl("mutect2_(vep|snpeff)[.]ann[.]vcf([.]gz)?$",
                             catalog$file, ignore.case = TRUE))) {
        shiny::p("No VEP or snpEff gene-annotated paired Mutect2 VCF is available for this run. ",
                 "The oncoplot needs gene consequences from a variant annotation step.")
      } else {
        shiny::p("Select Build oncoplot to summarize the annotated paired Mutect2 VCFs.")
      }
    } else sarek_report_oncoplot_ui(oncoplot),
    shiny::p("This view covers somatic small variants from paired Mutect2 only. ",
             "Structural and copy-number results remain in the output browser."),
    shiny::h3("Output inventory"),
    shiny::p(files, " readable files found in the run.",
             if (isTRUE(attr(catalog, "truncated"))) " Discovery stopped after 5,000 files." else ""),
    sarek_explorer_report_table(inventory),
    shiny::h4("Files by sample"),
    sarek_explorer_report_table(sample_counts),
    shiny::h3("Quality control"),
    shiny::p(sprintf(
      "%d QC files found: %d FastQC, %d fastp, and %d MultiQC outputs. Open individual reports in the QC tab.",
      NROW(qc),
      sum(grepl("fastqc", qc_names, ignore.case = TRUE)),
      sum(grepl("fastp", qc_names, ignore.case = TRUE)),
      sum(grepl("multiqc", qc_names, ignore.case = TRUE))
    )),
    if (NROW(qc)) shiny::p("Examples: ",
      paste(utils::head(basename(as.character(qc$file)), 8L), collapse = ", ")),
    if (!is.null(stats)) {
      shiny::div(shiny::h4("MultiQC general statistics"),
                 shiny::p("Values from the available MultiQC general statistics file; see the QC tab for the full report."),
                 sarek_explorer_report_table(stats, limit = 30L))
    },
    if (!is.null(fastp_stats)) {
      shiny::div(shiny::h4("fastp read statistics"),
                 shiny::p("Read counts before and after filtering, and Q30 rate after filtering, from available fastp JSON reports."),
                 sarek_explorer_report_table(fastp_stats, limit = 30L))
    },
    shiny::h3("Pipeline plots"),
    shiny::p(NROW(plots), " image or PDF outputs found. Browse and download each in the Plots tab."),
    if (NROW(plots)) shiny::p("Examples: ",
      paste(utils::head(basename(as.character(plots$file)), 12L), collapse = ", ")),
    shiny::h3("Interpretation scope"),
    shiny::p(switch(mode,
      germline = "Germline calls require assessment of inheritance, quality, phenotype, and evidence before interpretation.",
      tumor_only = "Tumor-only calls cannot reliably distinguish somatic from germline origin without a matched normal.",
      matched_tumor_normal = "Review tumor and normal evidence and caller filters before interpreting somatic candidates.",
      annotation_only = "Annotations describe provided variants; sequencing and variant calling were not run here.",
      "Review caller-specific outputs and input sample metadata before interpretation.")),
    shiny::p("Predicted impact is not clinical significance. Original QC reports and caller outputs remain the source for detailed review.")
  )
}

sarek_explorer_report_html <- function(report) {
  body <- htmltools::renderTags(report)$html
  paste0('<!doctype html><html><head><meta charset="utf-8"><title>Sarek output report</title>',
         '<style>body{font:16px/1.5 sans-serif;max-width:1100px;margin:2em auto;padding:0 1em;color:#243447}',
         '.sarek-report-table{overflow-x:auto;margin:1em 0}table{border-collapse:collapse;min-width:480px}',
         'td,th{padding:.5em;border:1px solid #ddd;text-align:left;white-space:nowrap}',
         'h2,h3{color:#183d62}h3{border-bottom:1px solid #d7e0ea;padding-bottom:.25em;margin-top:2em}',
         'img{max-width:100%;height:auto}</style></head><body>', body, '</body></html>')
}

# Scan complete annotated VCFs. Keep only a bounded table preview in memory;
# when csv_path is supplied, stream every matching record to that CSV.
sarek_explorer_clinsig_files <- function(catalog, output_dir,
                                         preview_limit = 1000L, csv_path = NULL) {
  preview_limit <- as.integer(preview_limit)
  if (is.na(preview_limit) || preview_limit < 0L) {
    stop("Invalid ClinVar preview limit.", call. = FALSE)
  }
  columns <- names(sarek_explorer_empty_variants())
  preview <- list()
  preview_count <- 0L
  matches <- 0L
  if (!is.null(csv_path)) {
    utils::write.table(
      sarek_explorer_empty_variants(), file = csv_path,
      sep = ",", row.names = FALSE, col.names = TRUE,
      quote = TRUE, qmethod = "double"
    )
  }

  for (i in seq_len(NROW(catalog))) {
    path <- sarek_explorer_file(catalog$path[[i]], output_dir)
    if (!nzchar(path)) {
      stop("An annotated VCF is unavailable or outside the run results directory: ",
           basename(catalog$file[[i]]), call. = FALSE)
    }
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
      gzfile(path, "rt")
    } else {
      base::file(path, "rt")
    }
    tryCatch({
      column_line <- ""
      csq_header <- ""
      for (j in seq_len(10000L)) {
        line <- readLines(con, n = 1L, warn = FALSE)
        if (!length(line)) break
        if (startsWith(line, "##INFO=<ID=CSQ,")) csq_header <- line
        if (startsWith(line, "#CHROM")) {
          column_line <- line
          break
        }
      }
      if (!nzchar(column_line)) {
        stop("No #CHROM header found in ", basename(path), call. = FALSE)
      }
      vcf_columns <- strsplit(column_line, "\t", fixed = TRUE)[[1]]
      samples <- if (length(vcf_columns) > 9L) {
        vcf_columns[10:length(vcf_columns)]
      } else {
        character()
      }
      selected_sample <- if (length(samples) == 1L) samples[[1L]] else ""
      sample_index <- match(selected_sample, vcf_columns)
      recorded_sample <- as.character(catalog$sample[[i]])
      if (is.na(recorded_sample)) recorded_sample <- ""
      csq_fields <- if (nzchar(csq_header)) {
        csq <- sub(".*Format: ", "", csq_header)
        csq <- sub('[">].*$', "", csq)
        strsplit(csq, "|", fixed = TRUE)[[1L]]
      } else {
        character()
      }

      repeat {
        chunk <- readLines(con, n = 5000L, warn = FALSE)
        if (!length(chunk)) break
        chunk <- chunk[grepl("CLNSIG=", chunk, fixed = TRUE)]
        if (!length(chunk)) next
        rows <- lapply(chunk, function(line) {
          fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
          if (length(fields) < 8L) {
            stop("Malformed VCF record in ", basename(path), call. = FALSE)
          }
          info <- fields[[8L]]
          clnsig <- sarek_explorer_info(info, "CLNSIG")
          if (!nzchar(clnsig) || identical(clnsig, ".")) return(NULL)

          ann <- sarek_explorer_info(info, "ANN")
          csq <- sarek_explorer_info(info, "CSQ")
          gene <- consequence <- impact <- ""
          if (nzchar(ann)) {
            parts <- strsplit(
              strsplit(ann, ",", fixed = TRUE)[[1L]][[1L]],
              "|", fixed = TRUE
            )[[1L]]
            if (length(parts) >= 4L) {
              consequence <- parts[[2L]]
              impact <- parts[[3L]]
              gene <- parts[[4L]]
            }
          } else if (nzchar(csq) && length(csq_fields)) {
            parts <- strsplit(
              strsplit(csq, ",", fixed = TRUE)[[1L]][[1L]],
              "|", fixed = TRUE
            )[[1L]]
            get_csq <- function(name) {
              index <- match(name, csq_fields)
              if (is.na(index) || index > length(parts)) "" else parts[[index]]
            }
            consequence <- get_csq("Consequence")
            impact <- get_csq("IMPACT")
            gene <- get_csq("SYMBOL")
          }

          format <- if (length(fields) >= 9L) fields[[9L]] else ""
          genotype <- if (!is.na(sample_index) &&
                          length(fields) >= sample_index) {
            fields[[sample_index]]
          } else {
            ""
          }
          data.frame(
            Sample = if (nzchar(selected_sample)) selected_sample else recorded_sample,
            Caller = sarek_explorer_caller(catalog$folder[[i]], catalog$file[[i]]),
            CHROM = fields[[1L]], POS = fields[[2L]],
            REF = fields[[4L]], ALT = fields[[5L]], FILTER = fields[[7L]],
            Genotype = sarek_explorer_format(format, genotype, "GT"),
            Depth = sarek_explorer_format(format, genotype, "DP"),
            VAF = sarek_explorer_format(format, genotype, "AF"),
            Gene = gene, Consequence = consequence, Impact = impact,
            ClinVar = clnsig,
            Review = sarek_explorer_info(info, "CLNREVSTAT"),
            Condition = sarek_explorer_info(info, "CLNDN"),
            Source = catalog$file[[i]],
            stringsAsFactors = FALSE, check.names = FALSE
          )
        })
        rows <- Filter(Negate(is.null), rows)
        if (!length(rows)) next
        matches <- matches + length(rows)
        if (!is.null(csv_path)) {
          utils::write.table(
            do.call(rbind, rows)[, columns, drop = FALSE],
            file = csv_path, append = TRUE, sep = ",",
            row.names = FALSE, col.names = FALSE,
            quote = TRUE, qmethod = "double"
          )
        }
        remaining <- preview_limit - preview_count
        if (remaining > 0L) {
          kept <- utils::head(rows, remaining)
          preview <- c(preview, kept)
          preview_count <- preview_count + length(kept)
        }
      }
    }, finally = close(con))
  }

  list(
    rows = if (length(preview)) {
      do.call(rbind, preview)[, columns, drop = FALSE]
    } else {
      sarek_explorer_empty_variants()
    },
    matches = matches, files = NROW(catalog), error = ""
  )
}

# Build a gene-by-tumor summary from complete annotated VCFs on demand.
sarek_report_consequence <- function(value) {
  terms <- strsplit(value, "&", fixed = TRUE)[[1]]
  if (any(terms %in% c("frameshift_variant", "stop_gained", "stop_lost", "start_lost"))) return("Truncating")
  if (any(terms %in% c("splice_acceptor_variant", "splice_donor_variant"))) return("Splice")
  if ("missense_variant" %in% terms) return("Missense")
  if (any(terms %in% c("inframe_insertion", "inframe_deletion", "protein_altering_variant"))) return("In-frame / other coding")
  ""
}

sarek_report_vcf_header <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  csq_fields <- character()
  has_ann <- FALSE
  has_columns <- FALSE
  for (i in seq_len(10000L)) {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    if (startsWith(line, "#CHROM")) { has_columns <- TRUE; break }
    if (startsWith(line, "##INFO=<ID=CSQ,") && grepl("Format: ", line, fixed = TRUE)) {
      format <- sub('".*$', "", sub("^.*Format: ", "", line))
      csq_fields <- strsplit(format, "|", fixed = TRUE)[[1]]
    }
    if (startsWith(line, "##INFO=<ID=ANN,")) has_ann <- TRUE
  }
  if (!has_columns) return(NULL)
  if (all(c("SYMBOL", "Consequence") %in% csq_fields)) {
    list(type = "CSQ", fields = csq_fields)
  } else if (has_ann) {
    list(type = "ANN", fields = character())
  } else NULL
}

sarek_report_oncoplot_data <- function(catalog, output_dir, manifest, top_n = 15L) {
  unavailable <- function(message) list(message = message)
  if (is.null(manifest) || !length(manifest$patients)) {
    return(unavailable("Sample roles are missing from the manifest; cannot assign calls to tumors."))
  }
  tumors <- unique(unlist(lapply(manifest$patients, function(patient) {
    if (!length(patient$samples)) return(character())
    vapply(Filter(function(sample) identical(tolower(sarek_text(sample$role)), "tumor"),
                  patient$samples), function(sample) sarek_text(sample$sample_id), character(1))
  }), use.names = FALSE))
  tumors <- tumors[nzchar(tumors)]
  if (!length(tumors)) return(unavailable("This plot requires tumor samples identified in the manifest."))

  name <- tolower(as.character(catalog$file))
  path <- tolower(as.character(catalog$path))
  paired <- grepl("_vs_", path, fixed = TRUE)
  annotated <- grepl("_vep\\.ann\\.vcf(\\.gz)?$|_snpeff\\.ann\\.vcf(\\.gz)?$", name)
  filtered <- grepl("filtered.*\\.vcf(\\.gz)?$", name)
  files <- catalog[grepl("mutect2", name, fixed = TRUE) & paired &
                   (annotated | filtered), , drop = FALSE]
  if (!NROW(files)) return(unavailable(
    "No paired Mutect2 VCFs are available. Gene annotation with VEP or snpEff is required."
  ))
  files$tumor <- vapply(files$sample_ids, function(ids) {
    if (is.na(ids)) return("")
    matched <- intersect(tumors, strsplit(ids, ";", fixed = TRUE)[[1]])
    if (length(matched) == 1L) matched[[1]] else ""
  }, character(1))
  files <- files[nzchar(files$tumor), , drop = FALSE]

  hits <- new.env(hash = TRUE, parent = emptyenv())
  severity <- c("Truncating", "Splice", "Missense", "In-frame / other coding")
  included <- character()
  sources <- character()
  for (tumor in tumors) {
    group <- files[files$tumor == tumor, , drop = FALSE]
    if (!NROW(group)) next
    priority <- ifelse(grepl("_vep\\.ann\\.vcf", tolower(group$file)), 0L,
                ifelse(grepl("_snpeff\\.ann\\.vcf", tolower(group$file)), 1L, 2L))
    group <- group[order(priority, group$file), , drop = FALSE]
    selected <- NULL
    header <- NULL
    for (i in seq_len(NROW(group))) {
      safe <- sarek_explorer_file(group$path[[i]], output_dir)
      if (!nzchar(safe)) next
      candidate <- sarek_report_vcf_header(safe)
      if (!is.null(candidate)) {
        selected <- safe
        header <- candidate
        break
      }
    }
    if (is.null(selected)) next
    included <- c(included, tumor)
    sources <- c(sources, basename(selected))
    con <- if (grepl("\\.gz$", selected, ignore.case = TRUE)) gzfile(selected, "rt") else file(selected, "rt")
    tryCatch({
      gene_index <- match("SYMBOL", header$fields)
      effect_index <- match("Consequence", header$fields)
      repeat {
        chunk <- readLines(con, n = 5000L, warn = FALSE)
        if (!length(chunk)) break
        chunk <- chunk[!startsWith(chunk, "#")]
        for (line in chunk) {
          columns <- strsplit(line, "\t", fixed = TRUE)[[1]]
          if (length(columns) < 8L) stop("Malformed VCF: ", basename(selected))
          if (!identical(columns[[7]], "PASS")) next
          value <- sarek_explorer_info(columns[[8]], header$type)
          if (!nzchar(value)) next
          for (entry in strsplit(value, ",", fixed = TRUE)[[1]]) {
            fields <- strsplit(entry, "|", fixed = TRUE)[[1]]
            if (identical(header$type, "CSQ")) {
              if (length(fields) < max(gene_index, effect_index)) next
              gene <- fields[[gene_index]]
              effect <- fields[[effect_index]]
            } else {
              if (length(fields) < 4L) next
              gene <- fields[[4]]
              effect <- fields[[2]]
            }
            kind <- sarek_report_consequence(effect)
            if (!nzchar(gene) || gene == "." || !nzchar(kind)) next
            key <- paste(gene, tumor, sep = "\t")
            previous <- hits[[key]]
            if (is.null(previous) || match(kind, severity) < match(previous, severity)) {
              hits[[key]] <- kind
            }
          }
        }
      }
    }, finally = close(con))
  }
  if (!length(included)) return(unavailable(
    "No gene-annotated paired Mutect2 VCFs were found. Run VEP or snpEff annotation on the called variants."
  ))
  values <- as.list(hits, all.names = TRUE)
  if (!length(values)) return(unavailable(
    "Annotated paired VCFs contain no PASS protein-altering or essential-splice calls."
  ))
  coordinates <- strsplit(names(values), "\t", fixed = TRUE)
  genes <- vapply(coordinates, `[[`, character(1), 1L)
  samples <- vapply(coordinates, `[[`, character(1), 2L)
  counts <- sort(table(genes), decreasing = TRUE)
  top <- names(counts)[seq_len(min(length(counts), as.integer(top_n)))]
  matrix <- matrix("", nrow = length(top), ncol = length(included),
                   dimnames = list(top, included))
  for (i in seq_along(values)) {
    if (genes[[i]] %in% top) matrix[genes[[i]], samples[[i]]] <- values[[i]]
  }
  sample_counts <- as.integer(table(factor(samples, levels = included)))
  list(matrix = matrix, gene_counts = as.integer(counts[top]),
       sample_counts = sample_counts,
       included = included, expected = length(tumors), sources = sources)
}

sarek_report_oncoplot_ui <- function(plot) {
  if (!is.null(plot$message)) return(shiny::p(plot$message))
  palette <- c("Truncating" = "#9f1239", "Splice" = "#9a3412",
               "Missense" = "#1d4ed8", "In-frame / other coding" = "#047857")
  values <- plot$matrix
  labels <- colnames(values)
  cell <- function(gene, sample) {
    kind <- values[gene, sample]
    shiny::tags$td(
      title = if (nzchar(kind)) paste(gene, "in", sample, "-", kind) else "No displayed alteration",
      style = paste0("height:30px;min-width:110px;text-align:center;border:1px solid #dce4ee;",
                     if (nzchar(kind)) paste0("background:", palette[[kind]], ";color:white;") else "background:#f1f5f9;"),
      if (nzchar(kind)) kind else ""
    )
  }
  shiny::div(
    shiny::p(
      length(plot$included), " of ", plot$expected,
      " manifest tumor samples have annotated paired Mutect2 VCFs. ",
      if (length(plot$included) == 1L)
        "This is a single-tumor alteration profile; cohort patterns require more tumors."
      else paste0("Showing the top ", nrow(values),
                  " genes by number of altered tumors.")
    ),
    shiny::div(style = "overflow-x:auto;max-width:100%;",
      shiny::tags$table(style = "border-spacing:3px;border-collapse:separate;",
        shiny::tags$thead(
          shiny::tags$tr(
            shiny::tags$th(style = "min-width:210px;text-align:left;", "Gene"),
            lapply(labels, function(label)
              shiny::tags$th(style = "min-width:110px;text-align:center;", label))
          ),
          shiny::tags$tr(
            shiny::tags$th(style = "white-space:nowrap;text-align:left;",
                           "Altered coding genes"),
            lapply(seq_along(labels), function(j)
              shiny::tags$th(style = "text-align:center;", plot$sample_counts[[j]]))
          )
        ),
        shiny::tags$tbody(lapply(seq_len(nrow(values)), function(i) {
          gene <- rownames(values)[[i]]
          shiny::tags$tr(shiny::tags$th(style = "white-space:nowrap;",
            paste0(gene, " (", plot$gene_counts[[i]], "/", length(labels), ")")),
            lapply(labels, function(sample) cell(gene, sample)))
        }))
      )
    ),
    shiny::div(style = "display:flex;gap:16px;flex-wrap:wrap;margin:12px 0;",
      lapply(names(palette), function(kind) shiny::tags$span(
        shiny::tags$span(style = paste0("display:inline-block;width:14px;height:14px;vertical-align:middle;margin-right:5px;background:", palette[[kind]], ";")),
        kind
      )),
      shiny::tags$span(
        shiny::tags$span(style = "display:inline-block;width:14px;height:14px;vertical-align:middle;margin-right:5px;background:#f1f5f9;border:1px solid #dce4ee;"),
        "No displayed alteration"
      )
    ),
    shiny::p("Cells show the highest-severity consequence for each gene and tumor. ",
             "Only PASS calls with coding or essential-splice consequences are shown; ",
             "the colors do not indicate clinical actionability."),
    shiny::tags$details(shiny::tags$summary("Source VCFs"),
      shiny::tags$ul(lapply(plot$sources, shiny::tags$li))),
    if (length(plot$included) < plot$expected) shiny::p(
      "Some tumor samples lack an annotated paired VCF and are omitted; blank cells only refer to included samples."
    )
  )
}
