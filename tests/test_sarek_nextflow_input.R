args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tests/test_sarek_nextflow_input.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

if (!requireNamespace("jsonlite", quietly = TRUE)) stop("The jsonlite package is required for this test.")

source(file.path(repo_root, "R", "sarek_manifest.R"), local = FALSE)
source(file.path(repo_root, "R", "sarek_nextflow_input.R"), local = FALSE)

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

expect_error <- function(expression, pattern, message) {
  observed <- tryCatch(
    {
      force(expression)
      ""
    },
    error = function(error) conditionMessage(error)
  )
  if (!nzchar(observed) || !grepl(pattern, observed, ignore.case = TRUE, perl = TRUE)) {
    stop(message, " Observed: ", observed, call. = FALSE)
  }
}

make_file <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.create(path)) stop("Could not create test file: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

make_sample <- function(sample_id, role, input_format, processing_state, files) {
  list(
    sample_id = sample_id,
    role = role,
    input_format = input_format,
    processing_state = processing_state,
    files = files
  )
}

make_patient <- function(patient_id, samples) {
  list(patient_id = patient_id, samples = samples, relationships = list())
}

make_manifest <- function(patients, mode = "germline", manifest_id = "nextflow_test") {
  list(
    schema_version = SAREK_MANIFEST_SCHEMA_VERSION,
    manifest_id = manifest_id,
    status = "confirmed",
    analysis = list(mode = mode, preset = "core"),
    patients = patients
  )
}

test_root <- tempfile("sarek-nextflow-test-")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

# Paired FASTQ manifests become one mapping row per lane.
fastq_r1 <- make_file(file.path(test_root, "germline_L001_R1.fastq.gz"))
fastq_r2 <- make_file(file.path(test_root, "germline_L001_R2.fastq.gz"))
fastq_manifest <- make_manifest(list(make_patient("patient1", list(make_sample(
  "sample1",
  "germline",
  "fastq",
  "unmapped",
  list(
    list(path = fastq_r1, lane = "L001", read = 1L),
    list(path = fastq_r2, lane = "L001", read = 2L)
  )
)))))

fastq_input <- sarek_build_nextflow_input(fastq_manifest)
assert_true(fastq_input$step == "mapping", "FASTQ input did not select the mapping step.")
assert_true(fastq_input$input_format == "fastq", "FASTQ input format was not preserved.")
assert_true(
  identical(names(fastq_input$samplesheet), c("patient", "sex", "status", "sample", "lane", "fastq_1", "fastq_2")),
  "FASTQ samplesheet columns do not match nf-core/sarek 3.9.0."
)
assert_true(NROW(fastq_input$samplesheet) == 1L, "One FASTQ pair should produce one samplesheet row.")
assert_true(fastq_input$samplesheet$sex[[1]] == "NA", "Missing sex should be represented as NA.")
assert_true(fastq_input$samplesheet$status[[1]] == 0L, "Germline status should map to 0.")
assert_true(fastq_input$samplesheet$fastq_1[[1]] == fastq_r1, "R1 path was not assigned to fastq_1.")
assert_true(fastq_input$samplesheet$fastq_2[[1]] == fastq_r2, "R2 path was not assigned to fastq_2.")
samplesheet_path <- file.path(test_root, "generated", "samplesheet.csv")
sarek_write_nextflow_samplesheet(fastq_input, samplesheet_path)
written <- utils::read.csv(samplesheet_path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = character(0))
assert_true(identical(names(written), names(fastq_input$samplesheet)), "Written samplesheet headers changed.")
assert_true(written$fastq_1[[1]] == fastq_r1 && written$fastq_2[[1]] == fastq_r2, "Written FASTQ paths changed.")

# A confirmed JSON manifest can be read and converted without Shiny.
manifest_json <- file.path(test_root, "confirmed.manifest.json")
jsonlite::write_json(fastq_manifest, manifest_json, auto_unbox = TRUE, pretty = TRUE, null = "null")
json_input <- sarek_manifest_json_to_nextflow_input(manifest_json)
assert_true(identical(json_input$samplesheet, fastq_input$samplesheet), "JSON manifest conversion changed the samplesheet.")

# Analysis-ready BAMs start at variant calling and require their BAI indexes.
normal_bam <- make_file(file.path(test_root, "patient_N.bam"))
normal_bai <- make_file(paste0(normal_bam, ".bai"))
tumor_bam <- make_file(file.path(test_root, "patient_T.bam"))
tumor_bai <- make_file(paste0(tumor_bam, ".bai"))
matched_manifest <- make_manifest(
  list(make_patient("patient2", list(
    make_sample("normal", "normal", "bam", "analysis_ready", list(list(path = normal_bam, index = normal_bai))),
    make_sample("tumor", "tumor", "bam", "analysis_ready", list(list(path = tumor_bam, index = tumor_bai)))
  ))),
  mode = "matched_tumor_normal",
  manifest_id = "matched_test"
)
matched_input <- sarek_build_nextflow_input(matched_manifest)
assert_true(matched_input$step == "variant_calling", "Analysis-ready BAM did not select variant_calling.")
assert_true(
  identical(names(matched_input$samplesheet), c("patient", "sex", "status", "sample", "bam", "bai")),
  "Processed BAM columns do not match nf-core/sarek 3.9.0."
)
assert_true(identical(matched_input$samplesheet$status, c(0L, 1L)), "Normal/tumor status mapping is incorrect.")

# Processing states choose the earliest safe Sarek restart point.
aligned_manifest <- matched_manifest
aligned_manifest$analysis$mode <- "germline"
aligned_manifest$patients <- list(make_patient("patient3", list(
  make_sample("aligned", "germline", "bam", "aligned", list(list(path = normal_bam, index = normal_bai)))
)))
assert_true(
  sarek_build_nextflow_input(aligned_manifest)$step == "markduplicates",
  "Aligned BAM did not select markduplicates."
)
aligned_manifest$patients[[1]]$samples[[1]]$processing_state <- "duplicate_marked"
assert_true(
  sarek_build_nextflow_input(aligned_manifest)$step == "prepare_recalibration",
  "Duplicate-marked BAM did not select prepare_recalibration."
)

# Annotation-only conversion uses the compact VCF samplesheet.
vcf <- make_file(file.path(test_root, "calls.vcf.gz"))
annotation_manifest <- make_manifest(
  list(make_patient("patient4", list(make_sample(
    "calls", "unknown", "vcf", "variant_calls", list(list(path = vcf))
  )))),
  mode = "annotation_only"
)
annotation_input <- sarek_build_nextflow_input(annotation_manifest)
assert_true(annotation_input$step == "annotate", "VCF input did not select annotate.")
assert_true(
  identical(names(annotation_input$samplesheet), c("patient", "sample", "vcf")),
  "Annotation samplesheet columns are incorrect."
)

# Known incompatibilities and unsafe mixed starts are stopped explicitly.
plain_fastq <- make_file(file.path(test_root, "plain_R1.fastq"))
plain_fastq_2 <- make_file(file.path(test_root, "plain_R2.fastq"))
plain_manifest <- make_manifest(list(make_patient("patient5", list(make_sample(
  "plain", "germline", "fastq", "unmapped",
  list(list(path = plain_fastq, lane = "L001", read = 1L), list(path = plain_fastq_2, lane = "L001", read = 2L))
)))))
expect_error(
  sarek_build_nextflow_input(plain_manifest),
  "gzip-compressed FASTQ",
  "Uncompressed FASTQ input should be rejected before launch."
)

missing_index_manifest <- aligned_manifest
missing_index_manifest$patients[[1]]$samples[[1]]$files[[1]]$index <- NULL
expect_error(
  sarek_build_nextflow_input(missing_index_manifest),
  "index is required",
  "Processed BAM without an index should be rejected."
)

mixed_manifest <- fastq_manifest
mixed_manifest$patients[[1]]$samples[[2]] <- make_sample(
  "ready", "germline", "bam", "analysis_ready", list(list(path = normal_bam, index = normal_bai))
)
expect_error(
  sarek_build_nextflow_input(mixed_manifest),
  "only one starting step",
  "A manifest with incompatible starting steps should be rejected."
)

ubam <- make_file(file.path(test_root, "reads.ubam"))
ubam_manifest <- make_manifest(list(make_patient("patient6", list(make_sample(
  "ubam", "germline", "ubam", "unmapped", list(list(path = ubam))
)))))
expect_error(
  sarek_build_nextflow_input(ubam_manifest),
  "requires a .bam filename",
  "The nf-core/sarek .ubam schema incompatibility should be explicit."
)

cat("Sarek manifest-to-Nextflow input tests: PASS\n")
