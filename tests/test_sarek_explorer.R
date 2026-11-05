args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tests/test_sarek_explorer.R"
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(root, "R", "sarek_manifest.R"), local = FALSE)
source(file.path(root, "R", "sarek_nextflow_input.R"), local = FALSE)
source(file.path(root, "R", "sarek_submission.R"), local = FALSE)
source(file.path(root, "R", "sarek_results.R"), local = FALSE)
source(file.path(root, "R", "sarek_explorer.R"), local = FALSE)
stopifnot(requireNamespace("shiny", quietly = TRUE))

location <- tempfile("sarek-explorer-")
dir.create(location)
on.exit(unlink(location, recursive = TRUE), add = TRUE)
output_dir <- file.path(location, "results")
dir.create(file.path(output_dir, "variant_calling", "mutect2"), recursive = TRUE)
dir.create(file.path(output_dir, "multiqc", "multiqc_data"), recursive = TRUE)
dir.create(file.path(output_dir, "fastp"), recursive = TRUE)
dir.create(file.path(output_dir, "fastqc"), recursive = TRUE)
vcf <- file.path(output_dir, "variant_calling", "mutect2", "S1.vcf")
writeLines(c(
  "##fileformat=VCFv4.2",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
  "chr1\t42\t.\tC\tT\t100\tPASS\tANN=T|missense_variant|HIGH|GENE1;CLNSIG=Pathogenic;CLNDN=Disease\tGT:DP:AF\t0/1:35:0.3",
  "chr1\t43\t.\tC\tA\t20\tLowQual\tANN=A|synonymous_variant|LOW|GENE2\tGT:DP:AF\t0/1:12:0.1"
), vcf)
plot <- file.path(output_dir, "variant_calling", "mutect2", "S1.cnv.png")
file.create(plot)
qc_html <- file.path(output_dir, "multiqc", "multiqc_report.html")
writeLines("<html><body>MultiQC</body></html>", qc_html)
qc_data <- file.path(output_dir, "multiqc", "multiqc_data", "multiqc_general_stats.txt")
writeLines(c("Sample\tCoverage\tMapped reads", "S1\t33.5\t1000"), qc_data)
fastp_json <- file.path(output_dir, "fastp", "S1_fastp.json")
writeLines('{"summary":{"before_filtering":{"total_reads":1000},"after_filtering":{"total_reads":900,"q30_rate":0.89}}}', fastp_json)
fastqc_html <- file.path(output_dir, "fastqc", "S1_fastqc.html")
writeLines("<html><body>FastQC</body></html>", fastqc_html)
catalog <- data.frame(
  file = c("S1.vcf", "S1.cnv.png", "multiqc_report.html", "multiqc_general_stats.txt",
           "S1_fastp.json", "S1_fastqc.html"),
  folder = c(rep("variant_calling/mutect2", 2), "multiqc", "multiqc/multiqc_data", "fastp", "fastqc"),
  sample = c("S1", "S1", "Shared", "Shared", "S1", "S1"),
  sample_ids = c("S1", "S1", "", "", "S1", "S1"),
  category = c("Variants", "Plots", rep("Quality control", 4)),
  path = c(vcf, plot, qc_html, qc_data, fastp_json, fastqc_html),
  stringsAsFactors = FALSE
)
rows <- sarek_explorer_variants(catalog, output_dir = output_dir)
stopifnot(NROW(rows) == 2L, rows$Gene[[1]] == "GENE1",
          rows$ClinVar[[1]] == "Pathogenic", rows$Genotype[[1]] == "0/1",
          identical(attr(rows, "limited"), FALSE))
partial <- sarek_explorer_variants(catalog, output_dir = output_dir, max_records = 1L)
stopifnot(NROW(partial) == 1L, isTRUE(attr(partial, "limited")))
stopifnot(NROW(sarek_explorer_plots(catalog)) == 1L)
stopifnot(NROW(sarek_explorer_qc(catalog, "S1")) == 4L)
stopifnot(NROW(sarek_explorer_qc(catalog, "another_sample")) == 2L)
outside <- file.path(location, "outside.png")
file.create(outside)
link <- file.path(output_dir, "external.png")
if (isTRUE(file.symlink(outside, link))) {
  stopifnot(identical(sarek_explorer_file(link, output_dir), ""))
}
report <- sarek_explorer_report_ui(
  list(run_id = "sample_run", output_dir = output_dir), catalog, rows,
  list(analysis = list(mode = "tumor_only"),
       assay = list(type = "WGS"), reference = list(assembly = "GRCh38"))
)
html <- sarek_explorer_report_html(report)
stopifnot(grepl("GENE1", html, fixed = TRUE),
          grepl("Tumor-only calls cannot reliably distinguish", html, fixed = TRUE),
          grepl("S1.cnv.png", html, fixed = TRUE),
          grepl("MultiQC general statistics", html, fixed = TRUE),
          grepl("33.5", html, fixed = TRUE),
          grepl("fastp read statistics", html, fixed = TRUE),
          grepl("89", html, fixed = TRUE),
          grepl("Variants by sample and caller", html, fixed = TRUE),
          !grepl(output_dir, html, fixed = TRUE))
without_multiqc <- catalog[!grepl("multiqc", catalog$folder, ignore.case = TRUE), , drop = FALSE]
without_multiqc_html <- sarek_explorer_report_html(sarek_explorer_report_ui(
  list(run_id = "fastp_fastqc_run", output_dir = output_dir), without_multiqc, rows))
stopifnot(grepl("fastp read statistics", without_multiqc_html, fixed = TRUE),
          grepl("FastQC", without_multiqc_html, fixed = TRUE),
          !grepl("MultiQC general statistics", without_multiqc_html, fixed = TRUE))
cat("Sarek results explorer tests: PASS\n")
