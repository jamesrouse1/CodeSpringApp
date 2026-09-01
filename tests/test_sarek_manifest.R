#!/usr/bin/env Rscript

source(file.path("R", "sarek_manifest.R"))
jsonlite::fromJSON(
  file.path("schemas", "sarek_manifest_v1.schema.json"),
  simplifyVector = FALSE
)

expect_error <- function(expr, pattern) {
  error <- tryCatch({
    force(expr)
    NULL
  }, error = function(e) e)
  stopifnot(inherits(error, "error"))
  stopifnot(grepl(pattern, conditionMessage(error), fixed = TRUE))
  invisible(error)
}

expect_invalid <- function(result, pattern) {
  stopifnot(!isTRUE(result$valid))
  stopifnot(any(grepl(pattern, result$errors, fixed = TRUE)))
  invisible(result)
}

test_root <- tempfile("sarek_manifest_test_")
dir.create(test_root)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

input_dir <- file.path(test_root, "inputs")
dir.create(input_dir)
files <- file.path(input_dir, c(
  "P001_T_L001_R1_001.fastq.gz",
  "P001_T_L001_R2_001.fastq.gz",
  "P001_N.bam",
  "P002.vcf.gz"
))
stopifnot(all(file.create(files)))
stopifnot(file.create(paste0(files[[3]], ".bai")))
stopifnot(file.create(paste0(files[[4]], ".tbi")))

discovery <- sarek_build_discovery_table(input_dir, allowed_roots = test_root)
stopifnot(NROW(discovery) == 4L)
stopifnot(sum(discovery$input_format == "fastq") == 2L)
stopifnot(sum(discovery$role == "tumor") == 2L)
stopifnot(sum(discovery$role == "normal") == 1L)
stopifnot(sum(discovery$role == "unknown") == 1L)
stopifnot(all(discovery$matched_normal_id[discovery$role == "tumor"] == "P001_N"))

# Facility-style FASTQs may put the read token after the chunk and include a
# dual index plus flowcell. These technical tokens must not create one sample
# per lane or make repeated L001 lanes from different flowcells collide.
facility_dir <- file.path(test_root, "facility_fastqs")
dir.create(facility_dir)
facility_units <- c(
  "23NMNMLT4_L001",
  "23W33VLT4_L006",
  "23W33VLT4_L007",
  "23W33VLT4_L008",
  "23WFVLLT4_L001"
)
facility_files <- unlist(lapply(c("BloodDNA", "PDODNA"), function(specimen) {
  dual_index <- if (identical(specimen, "BloodDNA")) {
    "ACAGTTACCT-GAACTGCCGG"
  } else {
    "AAACAGTGCA-AGATCTGTGG"
  }
  unlist(lapply(facility_units, function(unit) {
    file.path(
      facility_dir,
      paste0("ECO-18-", specimen, "_", dual_index, "_", unit, "_001.R", 1:2, ".fastq.gz")
    )
  }), use.names = FALSE)
}), use.names = FALSE)
stopifnot(all(file.create(facility_files)))

facility_discovery <- sarek_build_discovery_table(facility_dir, allowed_roots = test_root)
stopifnot(NROW(facility_discovery) == 20L)
stopifnot(identical(sort(unique(facility_discovery$patient_id)), "ECO-18"))
stopifnot(identical(
  sort(unique(facility_discovery$sample_id)),
  c("ECO-18-BloodDNA", "ECO-18-PDODNA")
))
stopifnot(identical(
  sort(unique(facility_discovery$lane)),
  sort(facility_units)
))
stopifnot(all(facility_discovery$role[facility_discovery$sample_id == "ECO-18-BloodDNA"] == "normal"))
stopifnot(all(facility_discovery$role[facility_discovery$sample_id == "ECO-18-PDODNA"] == "tumor"))
stopifnot(all(
  facility_discovery$matched_normal_id[facility_discovery$role == "tumor"] == "ECO-18-BloodDNA"
))
stopifnot(sarek_validate_confirmation_table(facility_discovery)$valid)

matched <- discovery[discovery$patient_id == "P001", , drop = FALSE]
matched$processing_state[matched$input_format == "bam"] <- "recalibrated"
validation <- sarek_validate_confirmation_table(matched)
stopifnot(validation$valid)
stopifnot(!length(sarek_validate_analysis_mode(matched, "matched_tumor_normal")))

manifest <- sarek_build_manifest(
  confirmation_table = matched,
  manifest_id = "verification_run",
  created_by = Sys.info()[["user"]],
  assay_type = "WGS",
  analysis_mode = "matched_tumor_normal",
  preset = "core",
  results_root = file.path(test_root, "results"),
  work_root = file.path(test_root, "work"),
  allowed_results_roots = test_root,
  allowed_work_roots = test_root
)
stopifnot(length(manifest$patients) == 1L)
stopifnot(length(manifest$patients[[1]]$samples) == 2L)
stopifnot(length(manifest$patients[[1]]$relationships) == 1L)
stopifnot(identical(
  manifest$patients[[1]]$relationships[[1]]$sample_ids,
  c("P001_T", "P001_N")
))

manifest_path <- sarek_write_manifest(manifest, file.path(test_root, "manifest.json"))
stopifnot(file.exists(manifest_path))
jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)

without_warning <- matched[, setdiff(names(matched), "warning"), drop = FALSE]
manifest_without_warning <- sarek_build_manifest(
  confirmation_table = without_warning,
  manifest_id = "without_warning",
  created_by = "tester",
  assay_type = "WGS",
  analysis_mode = "matched_tumor_normal",
  preset = "core",
  results_root = file.path(test_root, "results_without_warning"),
  work_root = file.path(test_root, "work_without_warning")
)
stopifnot(length(manifest_without_warning$patients) == 1L)

bad_include <- matched
bad_include$include[[1]] <- NA
expect_invalid(sarek_validate_confirmation_table(bad_include), "explicitly TRUE or FALSE")

mixed_role <- matched
mixed_role$role[[2]] <- "normal"
expect_invalid(sarek_validate_confirmation_table(mixed_role), "same sample role")

mixed_state <- matched
mixed_state$processing_state[[2]] <- "aligned"
expect_invalid(sarek_validate_confirmation_table(mixed_state), "same processing state")

incomplete_fastq <- matched
incomplete_fastq$include[incomplete_fastq$read == 2L & !is.na(incomplete_fastq$read)] <- FALSE
expect_invalid(sarek_validate_confirmation_table(incomplete_fastq), "exactly one R1 and one R2")

unpaired <- matched
unpaired$matched_normal_id[unpaired$role == "tumor"] <- ""
stopifnot(any(grepl(
  "requires a matched normal",
  sarek_validate_analysis_mode(unpaired, "matched_tumor_normal"),
  fixed = TRUE
)))

expect_error(
  sarek_build_manifest(
    matched,
    "relative_storage",
    "tester",
    "WGS",
    "matched_tumor_normal",
    "core",
    "relative/results",
    file.path(test_root, "work")
  ),
  "absolute server paths"
)

expect_error(
  sarek_build_manifest(
    matched,
    "bad_assay",
    "tester",
    "RNA-seq",
    "matched_tumor_normal",
    "core",
    file.path(test_root, "results_bad_assay"),
    file.path(test_root, "work_bad_assay")
  ),
  "assay type is invalid"
)

outside <- tempfile("sarek_outside_")
stopifnot(file.create(outside))
on.exit(unlink(outside, force = TRUE), add = TRUE)
expect_error(
  sarek_build_discovery_table(outside, allowed_roots = test_root),
  "outside the authorized server roots"
)

missing_index_file <- file.path(input_dir, "P003_N.bam")
stopifnot(file.create(missing_index_file))
missing_index <- sarek_build_discovery_table(missing_index_file)
missing_index$role <- "germline"
missing_index$patient_id <- "P003"
missing_index$sample_id <- "P003"
index_validation <- sarek_validate_confirmation_table(missing_index)
expect_invalid(index_validation, "must have a detected companion index")

cat("Sarek manifest tests: PASS\n")
