args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tests/test_sarek_clinvar.R"
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
for (name in c("sarek_manifest.R", "sarek_nextflow_input.R", "sarek_submission.R",
               "sarek_explorer.R", "sarek_clinvar.R")) {
  source(file.path(root, "R", name), local = FALSE)
}

location <- tempfile("sarek-clinvar-")
dir.create(location)
on.exit(unlink(location, recursive = TRUE), add = TRUE)
run_dir <- file.path(location, "run")
out <- file.path(run_dir, "results")
source_dir <- file.path(out, "variant_calling", "mutect2")
dir.create(source_dir, recursive = TRUE)
dir.create(file.path(run_dir, ".codespring"))
writeLines(c("field\tvalue", "status\tsubmitted"), file.path(run_dir, ".codespring", "submission.tsv"))
source_vcf <- file.path(source_dir, "S1.vcf")
writeLines(c("##fileformat=VCFv4.2", "##reference=GRCh38",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
  "chr1\t42\t.\tC\tT\t.\tPASS\t."), source_vcf)
reference <- file.path(location, "clinvar.vcf.gz")
con <- gzfile(reference, "wt")
writeLines(c("##fileformat=VCFv4.2", "##reference=GRCh38",
  '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Classification">',
  '##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="Review">',
  '##INFO=<ID=CLNDN,Number=.,Type=String,Description="Condition">',
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
  "chr1\t42\t.\tC\tT\t.\t.\tCLNSIG=Pathogenic"), con)
close(con)
file.create(paste0(reference, ".tbi"))
run <- list(run_id = "run", run_dir = run_dir, output_dir = out)
manifest <- list(reference = list(assembly = "GRCh38"))
rejected <- tryCatch({
  sarek_clinvar_submit(run, manifest, source_vcf, reference, "GRCh37",
    bcftools = "/bin/true", runner = function(command, args) "43210")
  FALSE
}, error = function(error) grepl("matches this run", conditionMessage(error)))
stopifnot(rejected)
submitted <- sarek_clinvar_submit(run, manifest, source_vcf, reference, "GRCh38",
  bcftools = "bcftools", runner = function(command, args) {
    stopifnot("--parsable" %in% args, any(grepl("--partition=", args, fixed = TRUE)))
    "43210"
  })
stopifnot(identical(submitted$job_id, "43210"),
  grepl("/annotation/clinvar/variant_calling/mutect2/", submitted$output, fixed = TRUE))
script <- list.files(file.path(run_dir, ".codespring", "clinvar"), pattern = "\\.sh$", full.names = TRUE)
stopifnot(length(script) == 1L,
  identical(readLines(script, n = 1L), "#!/bin/bash -l"),
  any(grepl("annotate --pair-logic exact", readLines(script), fixed = TRUE)),
  any(grepl("module load EB5", readLines(script), fixed = TRUE)),
  any(grepl("module load BCFtools/1.23.1-GCC-15.2.0", readLines(script), fixed = TRUE)))
jobs <- sarek_clinvar_jobs(run)
stopifnot(NROW(jobs) == 1L, jobs$Job[[1]] == "43210")
cat("Sarek ClinVar annotation tests: PASS\n")
