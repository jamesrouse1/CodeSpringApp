args <- commandArgs(trailingOnly = TRUE)

repo_root <- normalizePath(
  file.path(
    dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
    ".."
  ),
  mustWork = TRUE
)

if (length(args) && dir.exists(args[[1]])) {
  Sys.setenv(CSL_CODESPRINGLAB_ROOT = normalizePath(args[[1]], mustWork = TRUE))
}

app_env <- new.env(parent = globalenv())
sys.source(file.path(repo_root, "app.R"), envir = app_env)

# app.R sources runtime helper files globally for the Shiny application.
# Source this helper explicitly into the isolated test environment so its
# functions resolve the app helpers defined in app_env.
sys.source(
  file.path(repo_root, "R", "cutrun_nfcore_submission.R"),
  envir = app_env
)

assert <- function(value, message) {
  if (!isTRUE(value)) stop("ASSERTION FAILED: ", message, call. = FALSE)
}

root <- tempfile("cutrun_nfcore_")
fastq_dir <- file.path(root, "fastq")
data_dir <- file.path(root, "results")
manifest_dir <- file.path(root, "manifest")

dir.create(fastq_dir, recursive = TRUE)
dir.create(data_dir, recursive = TRUE)
dir.create(manifest_dir, recursive = TRUE)

files <- c(
  "H3K27ac_rep1_L001_R1.fastq.gz",
  "H3K27ac_rep1_L001_R2.fastq.gz",
  "H3K27ac_rep1_L002_R1.fastq.gz",
  "H3K27ac_rep1_L002_R2.fastq.gz",
  "H3K27ac_rep2_R1.fastq.gz",
  "H3K27ac_rep2_R2.fastq.gz",
  "IgG_rep1_R1.fastq.gz",
  "IgG_rep1_R2.fastq.gz",
  "IgG_rep2_R1.fastq.gz",
  "IgG_rep2_R2.fastq.gz"
)

invisible(file.create(file.path(fastq_dir, files)))

design <- data.frame(
  sample = c("H3K27ac_rep1", "H3K27ac_rep2", "IgG_rep1", "IgG_rep2"),
  cell_type = c("Epithelial", "Epithelial", "Epithelial", "Epithelial"),
  mark = c("H3K27ac", "H3K27ac", "IgG", "IgG"),
  target_class = c("histone_narrow", "histone_narrow", "control", "control"),
  seacr_stringency = c("auto", "auto", "auto", "auto"),
  condition = c("Vehicle", "Vehicle", "Vehicle", "Vehicle"),
  replicate = c(1, 2, 1, 2),
  control_sample = c("IgG_rep1", "IgG_rep2", "", ""),
  filename = c(
    paste(
      "H3K27ac_rep1_L001_R1.fastq.gz,H3K27ac_rep1_L001_R2.fastq.gz",
      "H3K27ac_rep1_L002_R1.fastq.gz,H3K27ac_rep1_L002_R2.fastq.gz",
      sep = ";"
    ),
    "H3K27ac_rep2_R1.fastq.gz,H3K27ac_rep2_R2.fastq.gz",
    "IgG_rep1_R1.fastq.gz,IgG_rep1_R2.fastq.gz",
    "IgG_rep2_R1.fastq.gz,IgG_rep2_R2.fastq.gz"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

design_path <- file.path(manifest_dir, "design_matrix.txt")
utils::write.table(
  design,
  design_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

project <- list(
  id = "cutrun/nfcore_test",
  name = "nfcore_test",
  analysis = "CUT&RUN",
  analysis_key = "cutrun",
  cutrun_backend = "nfcore",
  paired_end = TRUE,
  fastq_dir = fastq_dir,
  fastq_dirs = fastq_dir,
  data_dir = data_dir,
  design_matrix_path = design_path
)

generated <- app_env$cutrun_nfcore_samplesheet(project, design)

sheet <- generated$samplesheet
mapping <- generated$mapping

assert(
  identical(
    names(sheet),
    c("group", "replicate", "fastq_1", "fastq_2", "control")
  ),
  "nf-core samplesheet header is incorrect"
)

assert(
  NROW(sheet) == 5L,
  "two-lane target plus three single-lane samples should create five samplesheet rows"
)

assert(
  sum(mapping$codespring_sample == "H3K27ac_rep1") == 2L,
  "pooled lanes were not preserved as nf-core technical-replicate rows"
)

assert(
  length(unique(sheet$group[grepl("^target__", sheet$group)])) == 1L,
  "biological target replicates did not share one nf-core group"
)

assert(
  length(unique(sheet$group[grepl("^control__", sheet$group)])) == 1L,
  "IgG biological replicates did not share one nf-core control group"
)

target_rows <- grepl("^target__", sheet$group)

assert(
  all(nzchar(sheet$control[target_rows])),
  "target rows are missing their nf-core control group"
)

assert(
  all(sheet$control[target_rows] %in% sheet$group),
  "target control group does not correspond to a samplesheet group"
)

control_rows <- grepl("^control__", sheet$group)

assert(
  all(!nzchar(sheet$control[control_rows])),
  "control rows must not themselves specify a control"
)

assert(
  all(file.exists(sheet$fastq_1)) &&
    all(file.exists(sheet$fastq_2)),
  "generated samplesheet does not contain real FASTQ paths"
)

assert(
  all(startsWith(sheet$fastq_1, "/")) &&
    all(startsWith(sheet$fastq_2, "/")),
  "nf-core samplesheet FASTQs must be absolute paths"
)

written <- app_env$cutrun_nfcore_write_inputs(
  project,
  run_dir = file.path(data_dir, "nfcore_cutandrun")
)

assert(file.exists(written$samplesheet), "samplesheet.csv was not written")
assert(file.exists(written$mapping), "sample_mapping.tsv was not written")

disk_sheet <- read.csv(
  written$samplesheet,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

assert(
  identical(
    names(disk_sheet),
    c("group", "replicate", "fastq_1", "fastq_2", "control")
  ),
  "written nf-core CSV header is incorrect"
)

assert(
  NROW(disk_sheet) == 5L,
  "written nf-core samplesheet row count changed"
)

unlink(root, recursive = TRUE, force = TRUE)

cat("CUT&RUN nf-core samplesheet: PASS\n")
