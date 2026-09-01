#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "bin/sarek-manifest-to-input.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    paste(
      "Usage:",
      "Rscript --vanilla bin/sarek-manifest-to-input.R",
      "<confirmed-manifest.json> <output-samplesheet.csv>"
    ),
    call. = FALSE
  )
}

source(file.path(repo_root, "R", "sarek_manifest.R"), local = FALSE)
source(file.path(repo_root, "R", "sarek_nextflow_input.R"), local = FALSE)

nextflow_input <- sarek_manifest_json_to_nextflow_input(args[[1]])
samplesheet_path <- sarek_write_nextflow_samplesheet(nextflow_input, args[[2]])

cat("Sarek Nextflow input created\n")
cat("  Manifest ID: ", nextflow_input$manifest_id, "\n", sep = "")
cat("  Pipeline: ", nextflow_input$pipeline, " ", nextflow_input$pipeline_version, "\n", sep = "")
cat("  Step: ", nextflow_input$step, "\n", sep = "")
cat("  Input format: ", nextflow_input$input_format, "\n", sep = "")
cat("  Rows: ", NROW(nextflow_input$samplesheet), "\n", sep = "")
cat("  Samplesheet: ", samplesheet_path, "\n", sep = "")
if (length(nextflow_input$warnings)) {
  cat("Warnings:\n")
  cat(paste0("  - ", nextflow_input$warnings, collapse = "\n"), "\n", sep = "")
}
