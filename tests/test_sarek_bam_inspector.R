args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tests/test_sarek_bam_inspector.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

source(file.path(repo_root, "R", "sarek_manifest.R"), local = FALSE)
source(file.path(repo_root, "R", "sarek_bam_inspector.R"), local = FALSE)

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

test_root <- tempfile("sarek-bam-inspector-")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

fake_samtools <- file.path(test_root, "samtools")
writeLines(
  c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "case \"$1\" in",
    "  quickcheck)",
    "    exit 0",
    "    ;;",
    "  view)",
    "    printf '@HD\\tVN:1.6\\tSO:coordinate\\n'",
    "    printf '@SQ\\tSN:chr1\\tLN:248956422\\n'",
    "    printf '@SQ\\tSN:chrM\\tLN:16569\\n'",
    "    printf '@RG\\tID:rg1\\tSM:TUMOR_01\\n'",
    "    printf '@PG\\tID:markdup\\tPN:MarkDuplicates\\n'",
    "    printf '@PG\\tID:bqsr\\tPN:ApplyBQSR\\tCL:gatk ApplyBQSR\\n'",
    "    ;;",
    "  idxstats)",
    "    printf 'chr1\\t248956422\\t100\\t0\\n'",
    "    ;;",
    "  *)",
    "    exit 2",
    "    ;;",
    "esac"
  ),
  fake_samtools
)
Sys.chmod(fake_samtools, mode = "0755")

bam <- file.path(test_root, "patient_T.bam")
bai <- paste0(bam, ".bai")
invisible(file.create(bam, bai))

inspection <- sarek_inspect_bam(bam, bai, samtools = fake_samtools)
assert_true(inspection$status == "passed", "A valid lightweight BAM inspection did not pass.")
assert_true(isTRUE(inspection$quickcheck_ok), "samtools quickcheck was not recorded as successful.")
assert_true(isTRUE(inspection$index_ok), "The usable BAM index was not recorded.")
assert_true(identical(inspection$sample_ids, "TUMOR_01"), "The BAM @RG sample ID was not extracted.")
assert_true(inspection$sort_order == "coordinate", "The BAM sort order was not extracted.")
assert_true(isTRUE(inspection$grch38_compatible), "GATK.GRCh38-compatible reference evidence was not recognized.")
assert_true(
  inspection$recommendation == "recalibrated" && inspection$confidence == "high",
  "Applied BQSR provenance did not produce the expected processing recommendation."
)

markdup_only <- c(
  "@HD\tVN:1.6\tSO:coordinate",
  "@SQ\tSN:chr1\tLN:248956422",
  "@PG\tID:markdup\tPN:MarkDuplicates",
  "@PG\tID:recal_table\tPN:BaseRecalibrator"
)
markdup_state <- sarek_infer_bam_processing_state(markdup_only)
assert_true(
  markdup_state$recommendation == "duplicate_marked",
  "BaseRecalibrator alone was incorrectly treated as proof that BQSR was applied."
)

plain_contigs <- sarek_bam_reference_evidence(c("@SQ\tSN:1\tLN:248956422"))
assert_true(
  is.na(plain_contigs$grch38_compatible) && plain_contigs$contig_style == "non-chr-prefixed",
  "A non-chr GRCh38-length BAM did not surface its possible GATK resource naming mismatch."
)

missing_index <- sarek_inspect_bam(bam, "", samtools = fake_samtools)
assert_true(missing_index$status == "failed", "A missing BAM index was not rejected by inspection.")

unavailable <- sarek_inspect_bam(bam, bai, samtools = file.path(test_root, "missing-samtools"))
assert_true(unavailable$status == "unavailable", "An unavailable samtools executable was not reported cleanly.")

table <- sarek_build_discovery_table(bam, allowed_roots = test_root)
table <- sarek_auto_inspect_bams(table, samtools = fake_samtools, max_bams = 1L)
assert_true(table$inspection_status[[1]] == "passed", "Automatic BAM inspection was not attached to discovery.")
assert_true(
  table$processing_recommendation[[1]] == "recalibrated",
  "The automatic processing recommendation was not attached to discovery."
)
validation_before_confirmation <- sarek_validate_confirmation_table(table)
assert_true(
  !isTRUE(validation_before_confirmation$valid),
  "An inspected BAM was accepted before the recommended processing state was explicitly confirmed."
)

table$role[[1]] <- "tumor"
table$processing_state[[1]] <- table$processing_recommendation[[1]]
validation_after_confirmation <- sarek_validate_confirmation_table(table)
assert_true(
  isTRUE(validation_after_confirmation$valid),
  paste("A valid inspected and corroborated BAM did not pass validation:", paste(validation_after_confirmation$errors, collapse = "; "))
)

cat("Sarek BAM inspector tests: PASS\n")
