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

assert <- function(value, message) {
  if (!isTRUE(value)) stop("ASSERTION FAILED: ", message, call. = FALSE)
}

native_legacy <- list(
  analysis_key = "cutrun",
  analysis = "CUT&RUN"
)

native_explicit <- list(
  analysis_key = "cutrun",
  analysis = "CUT&RUN",
  cutrun_backend = "native"
)

nfcore_explicit <- list(
  analysis_key = "cutrun",
  analysis = "CUT&RUN",
  cutrun_backend = "nfcore"
)

rna_project <- list(
  analysis_key = "rna",
  analysis = "RNA-seq"
)

assert(
  identical(app_env$cutrun_backend(native_legacy), "native"),
  "legacy CUT&RUN projects without backend metadata must default to native"
)

assert(
  identical(app_env$cutrun_backend(native_explicit), "native"),
  "explicit native CUT&RUN backend was not preserved"
)

assert(
  identical(app_env$cutrun_backend(nfcore_explicit), "nfcore"),
  "explicit nf-core CUT&RUN backend was not preserved"
)

assert(
  app_env$is_native_cutrun_project(native_legacy),
  "legacy CUT&RUN project was not recognized as native"
)

assert(
  app_env$is_nfcore_cutrun_project(nfcore_explicit),
  "nf-core CUT&RUN project was not recognized as nf-core"
)

assert(
  identical(app_env$cutrun_backend(rna_project), ""),
  "non-CUT&RUN projects must not receive a CUT&RUN backend"
)

assert(
  identical(app_env$normalize_cutrun_backend("garbage"), "native"),
  "invalid backend values must safely fall back to native"
)

writer_source <- paste(
  deparse(body(app_env$write_project_config), width.cutoff = 500L),
  collapse = "\n"
)

loader_source <- paste(
  deparse(body(app_env$legacy_project_from_config), width.cutoff = 500L),
  collapse = "\n"
)

creator_source <- paste(
  deparse(body(app_env$new_project_from_inputs), width.cutoff = 500L),
  collapse = "\n"
)

assert(
  grepl("cutrun_backend", writer_source, fixed = TRUE),
  "project config writer does not persist CUT&RUN backend"
)

assert(
  grepl("cutrun_backend", loader_source, fixed = TRUE),
  "project config loader does not restore CUT&RUN backend"
)

assert(
  grepl("new_cutrun_backend", creator_source, fixed = TRUE),
  "new-project constructor does not support a CUT&RUN backend"
)

cat("CUT&RUN backend abstraction: PASS\n")
