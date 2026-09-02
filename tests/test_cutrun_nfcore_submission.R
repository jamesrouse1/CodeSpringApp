args <- commandArgs(trailingOnly = TRUE)

repo_root <- normalizePath(
  ".",
  winslash = "/",
  mustWork = TRUE
)

app_env <- new.env(parent = globalenv())

sys.source(
  file.path(repo_root, "app.R"),
  envir = app_env
)

# app.R sources runtime helpers globally. Re-source this one into
# the isolated test environment so dependencies resolve locally.
sys.source(
  file.path(
    repo_root,
    "R",
    "cutrun_nfcore_submission.R"
  ),
  envir = app_env
)

assert <- function(value, message) {
  if (!isTRUE(value)) {
    stop(message, call. = FALSE)
  }
}

root <- tempfile("cutrun_nfcore_submission_")
dir.create(root, recursive = TRUE)

on.exit(
  unlink(root, recursive = TRUE, force = TRUE),
  add = TRUE
)

# -------------------------------------------------------------------------
# Fake CodeSpring project
# -------------------------------------------------------------------------

project <- list(
  id = "nfcore-cutrun-test",
  name = "nfcore-cutrun-test",
  analysis_key = "cutrun",
  analysis = "CUT&RUN",
  cutrun_backend = "nfcore",
  data_dir = root
)

# -------------------------------------------------------------------------
# Fake but structurally complete mouse reference
# -------------------------------------------------------------------------

reference_root <- file.path(root, "reference")
bowtie2_dir <- file.path(
  reference_root,
  "bowtie2_index"
)

dir.create(
  bowtie2_dir,
  recursive = TRUE
)

fasta <- file.path(
  reference_root,
  "GRCm39.primary_assembly.genome.fa"
)

gtf <- file.path(
  reference_root,
  "gencode.vM39.primary_assembly.annotation.gtf"
)

writeLines(">chr1\nACGT", fasta)
writeLines(
  'chr1\ttest\tgene\t1\t4\t.\t+\t.\tgene_id "x";',
  gtf
)

bowtie2_prefix <- file.path(
  bowtie2_dir,
  "GRCm39_gencodeM39"
)

index_suffixes <- c(
  ".1.bt2",
  ".2.bt2",
  ".3.bt2",
  ".4.bt2",
  ".rev.1.bt2",
  ".rev.2.bt2"
)

invisible(
  file.create(
    paste0(
      bowtie2_prefix,
      index_suffixes
    )
  )
)

app_env$genome_species <- function(project) {
  "mouse"
}

app_env$cutrun_reference_resources <- function(project) {
  list(
    bowtie2_index = bowtie2_prefix,
    gtf = gtf,
    chrom_sizes = "",
    macs2_genome = "mm",
    blacklist = ""
  )
}

# -------------------------------------------------------------------------
# Mock the separately tested samplesheet layer
# -------------------------------------------------------------------------

fastq_1 <- file.path(root, "target_R1.fastq.gz")
fastq_2 <- file.path(root, "target_R2.fastq.gz")
control_1 <- file.path(root, "igg_R1.fastq.gz")
control_2 <- file.path(root, "igg_R2.fastq.gz")

invisible(
  file.create(
    fastq_1,
    fastq_2,
    control_1,
    control_2
  )
)

generated <- list(
  samplesheet = data.frame(
    group = c(
      "target__epithelial_h3k27ac_vehicle",
      "control__epithelial_igg_vehicle"
    ),
    replicate = c(1L, 1L),
    fastq_1 = c(fastq_1, control_1),
    fastq_2 = c(fastq_2, control_2),
    control = c(
      "control__epithelial_igg_vehicle",
      ""
    ),
    stringsAsFactors = FALSE
  ),
  mapping = data.frame(
    codespring_sample = c(
      "H3K27ac_rep1",
      "IgG_rep1"
    ),
    nfcore_group = c(
      "target__epithelial_h3k27ac_vehicle",
      "control__epithelial_igg_vehicle"
    ),
    replicate = c(1L, 1L),
    lane = c(1L, 1L),
    role = c("target", "control"),
    control_sample = c("IgG_rep1", ""),
    control_group = c(
      "control__epithelial_igg_vehicle",
      ""
    ),
    fastq_1 = c(fastq_1, control_1),
    fastq_2 = c(fastq_2, control_2),
    stringsAsFactors = FALSE
  )
)

app_env$cutrun_nfcore_samplesheet <- function(
  project,
  design = NULL
) {
  generated
}

# -------------------------------------------------------------------------
# Fake runtime. We test generated commands, not Nextflow itself here.
# -------------------------------------------------------------------------

launcher <- file.path(
  root,
  "nextflow-sarek"
)

writeLines(
  c(
    "#!/usr/bin/env bash",
    "exit 0"
  ),
  launcher
)

Sys.chmod(
  launcher,
  mode = "0700"
)

cluster_config <- file.path(
  root,
  "cshl_slurm.config"
)

writeLines(
  c(
    "process {",
    "  executor = 'slurm'",
    "  queue = 'cpuq'",
    "}",
    "executor {",
    "  queueSize = 20",
    "}"
  ),
  cluster_config
)

runtime <- list(
  launcher = launcher,
  config = cluster_config,
  nxf_home = file.path(
    root,
    "nextflow_home"
  ),
  singularity_cache = file.path(
    root,
    "singularity_cache"
  )
)

# -------------------------------------------------------------------------
# Build bundle
# -------------------------------------------------------------------------

bundle <- app_env$cutrun_nfcore_build_bundle(
  project = project,
  run_id = "bundle_test",
  runtime = runtime
)

paths <- bundle$paths

required_files <- c(
  paths$samplesheet_path,
  paths$mapping_path,
  paths$params_path,
  paths$run_config,
  paths$launch_script
)

assert(
  all(file.exists(required_files)),
  "CUT&RUN nf-core bundle did not create all required files."
)

assert(
  file.access(paths$launch_script, 1) == 0,
  "CUT&RUN nf-core launch script is not executable."
)

params <- jsonlite::read_json(
  paths$params_path,
  simplifyVector = TRUE
)

assert(
  identical(
    params$normalisation_mode,
    "CPM"
  ),
  "Initial nf-core CUT&RUN normalisation is not CPM."
)

assert(
  identical(
    params$peakcaller,
    "seacr"
  ),
  "SEACR is not the primary/default nf-core CUT&RUN caller."
)

assert(
  identical(
    params$use_control,
    TRUE
  ),
  "Matched CUT&RUN controls were not enabled."
)

assert(
  identical(
    params$fasta,
    normalizePath(
      fasta,
      winslash = "/",
      mustWork = TRUE
    )
  ),
  "nf-core CUT&RUN FASTA reference is incorrect."
)

assert(
  identical(
    params$gtf,
    normalizePath(
      gtf,
      winslash = "/",
      mustWork = TRUE
    )
  ),
  "nf-core CUT&RUN GTF reference is incorrect."
)

assert(
  identical(
    params$bowtie2,
    normalizePath(
      bowtie2_dir,
      winslash = "/",
      mustWork = TRUE
    )
  ),
  "nf-core CUT&RUN Bowtie2 directory is incorrect."
)

config_text <- paste(
  readLines(
    paths$run_config,
    warn = FALSE
  ),
  collapse = "\n"
)

assert(
  grepl(
    cluster_config,
    config_text,
    fixed = TRUE
  ),
  "Generated config does not include the CSHL Slurm config."
)

assert(
  grepl(
    "resourceLimits = [ time: 48.h ]",
    config_text,
    fixed = TRUE
  ),
  "Generated config does not enforce the 48-hour ceiling."
)

launch_text <- paste(
  readLines(
    paths$launch_script,
    warn = FALSE
  ),
  collapse = "\n"
)

assert(
  grepl(
    "NXF_VER='25.10.2'",
    launch_text,
    fixed = TRUE
  ),
  "Launch script does not pin Nextflow 25.10.2."
)

assert(
  grepl(
    "nf-core/cutandrun",
    launch_text,
    fixed = TRUE
  ),
  "Launch script does not run nf-core/cutandrun."
)

assert(
  grepl(
    "-r '3.2.2'",
    launch_text,
    fixed = TRUE
  ),
  "Launch script does not pin nf-core/cutandrun 3.2.2."
)

assert(
  grepl(
    "-params-file",
    launch_text,
    fixed = TRUE
  ) &&
  grepl(
    paths$params_path,
    launch_text,
    fixed = TRUE
  ),
  "Launch script does not use the generated params.json."
)

assert(
  grepl(
    "-profile singularity",
    launch_text,
    fixed = TRUE
  ),
  "Launch script does not use the Singularity profile."
)

cat("CUT&RUN nf-core submission bundle: PASS\n")
