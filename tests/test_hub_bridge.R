# Tests for R/hub_bridge.R (codespring-run), deep links and one-time open tickets.
# Usage: Rscript tests/test_hub_bridge.R [/path/to/CodeSpringLab]
# Nothing is submitted to Slurm: start() is exercised with a fake controller starter.

args <- commandArgs(trailingOnly = TRUE)
repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), ".."), mustWork = TRUE)
lab_root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(file.path(repo_root, "..", "CodeSpringLab"), mustWork = TRUE)
test_home <- tempfile("hub_bridge_home_")
dir.create(test_home)
Sys.setenv(CSL_CODESPRINGLAB_ROOT = lab_root, CSL_WEB_APP_ROOT = repo_root, CSL_WEB_HOME = file.path(path.expand("~"), ".codespringweb_hubtest"),
           CSL_RESULTS_ROOT = file.path(test_home, "results"), CSL_HUB_POLL_SECONDS = "1")
on.exit(unlink(c(test_home, Sys.getenv("CSL_WEB_HOME")), recursive = TRUE), add = TRUE)

app <- new.env(parent = globalenv())
old <- setwd(repo_root)
invisible(utils::capture.output(sys.source("app.R", envir = app)))
sys.source(file.path("R", "hub_bridge.R"), envir = app)
setwd(old)

assert <- function(value, message) if (!isTRUE(value)) stop("ASSERTION FAILED: ", message, call. = FALSE)
pass <- function(message) cat("ok -", message, "\n")

fastq <- file.path(lab_root, "scripts_DoNotTouch", "test", "fastq")
design <- file.path(lab_root, "scripts_DoNotTouch", "test", "manifest", "design_matrix.txt")
base_spec <- function(...) {
  spec <- list(pipeline = "rnaseq", project_name = "hubtest_rna", fastq_dirs = list(fastq), paired_end = TRUE,
               species = "mouse", design_matrix_path = design,
               deseq2 = list(column = "treatment", reference = "control", comparison = "treated"))
  changes <- list(...)
  for (key in names(changes)) spec[key] <- list(changes[[key]])  # replace whole fields (NULL removes)
  Filter(Negate(is.null), spec)
}

# ---- catalog ------------------------------------------------------------------------------------
cat_ <- app$hub_catalog()
assert(identical(cat_$pipelines[[1]]$id, "rnaseq") && isTRUE(cat_$pipelines[[1]]$runnable_from_hub), "RNA-seq is the runnable pipeline")
assert(!any(vapply(cat_$pipelines[-1], function(p) isTRUE(p$runnable_from_hub), logical(1))), "other pipelines are app-only for now")
assert("mouse_gencodeM39" %in% vapply(cat_$rnaseq$references, `[[`, "", "genome_version"), "references come from the app's catalog")
assert(!"__custom__" %in% unlist(cat_$rnaseq$adapters_r1), "custom adapter sequences are not offered to the Hub")
pass("catalog lists BSR pipelines and the app's own options")

# ---- check ----------------------------------------------------------------------------------------
ok <- app$hub_check_rnaseq(base_spec())
assert(isTRUE(ok$ok), paste("the bundled example passes:", paste(ok$errors, collapse = "; ")))
plan <- app$hub_plan(ok)
assert(identical(vapply(plan$steps, `[[`, "", "step"), c("Cutadapt", "FastQC", "STAR", "featureCounts", "DESeq2")), "plan has the BSR RNA-seq steps in order")
star <- plan$steps[[3]]
assert(identical(star$resources_per_job[["mem-per-cpu"]], "50G") && identical(star$resources_per_job[["cpus-per-task"]], "4"), "resources are read from CodeSpringLab's #SBATCH headers")
assert(startsWith(plan$results_dir, Sys.getenv("CSL_RESULTS_ROOT")), "results go to the user's results root")
pass("check accepts the example and plans CodeSpringLab's steps and resources")

refused <- function(spec, pattern, label) {
  res <- app$hub_check_rnaseq(spec)
  assert(!isTRUE(res$ok) && any(grepl(pattern, res$errors)), paste(label, "->", paste(res$errors, collapse = " | ")))
}
refused(base_spec(results_root = "/tmp"), "results_root can't be chosen", "caller can't choose the results folder")
refused(base_spec(project_name = "../x"), "project_name", "unsafe project names refused")
refused(base_spec(genome_version = "mm10"), "genome_version", "only catalog references")
refused(base_spec(species = "zebrafish"), "species", "only supported species")
refused(base_spec(trim = list(adapter_r1 = "AGATCGGAAGAGC")), "adapter_r1", "only adapter presets")
refused(base_spec(gsea = list(geneset = "my_sets.gmt")), "gsea.geneset", "only listed gene sets")
refused(base_spec(deseq2 = list(column = "treatment", reference = "control", comparison = "control")), "differ", "comparison must differ")
refused(base_spec(fastq_dirs = list("/nonexistent")), "FASTQ folders", "FASTQ folders must exist")
refused(base_spec(design_matrix_path = NULL, design = list(columns = list("sample", "treatment", "filename"),
  rows = list(list("a; rm -rf ~", "x", "control1_R1_001.fastq.gz,control1_R2_001.fastq.gz"), list("b", "y", "../../etc/passwd")))),
  "Sample names|not paths", "shell-unsafe sample names and paths refused")
refused(base_spec(design_matrix_path = NULL, design = list(columns = list("sample", "treatment", "filename"),
  rows = list(list("a", "x", "missing_R1.fq.gz,missing_R2.fq.gz"), list("b", "x", "control1_R1_001.fastq.gz,control1_R2_001.fastq.gz"))), deseq2 = NULL),
  "aren't in the FASTQ folders", "every read file must exist")
pass("check refuses anything outside the BSR pipeline's options")

# ---- start (fake controller), status, resume ---------------------------------------------------
started_commands <- character(0)
starter <- function(command) { started_commands <<- c(started_commands, command); "424242" }
res <- app$hub_start(base_spec(), starter = starter)
assert(isTRUE(res$ok) && identical(res$controller_pid, "424242"), "start launches one controller")
assert(grepl("codespring-run' drive", started_commands[[1]], fixed = TRUE) && grepl("setsid", started_commands[[1]]), "controller is the detached runner")
p <- app$discover_projects()[[res$project_id]]
assert(!is.null(p), "project appears in the user's CodeSpringApp project list")
assert(file.exists(file.path(p$data_dir, "manifest", "design_matrix.txt")), "design copied into the project")
again <- app$hub_check_rnaseq(base_spec())
assert(!isTRUE(again$ok) && grepl("already exists", again$errors), "an existing project is never overwritten")
st <- app$hub_status(res$project_id)
assert(identical(st$state, "stopped"), "a controller that isn't running is reported as stopped")
resumed <- app$hub_resume(res$project_id, starter = starter)
assert(isTRUE(resumed$ok) && length(started_commands) == 2, "a stopped run can be resumed")
failing <- app$hub_start(base_spec(project_name = "hubtest_fail"), starter = function(command) "no pid")
assert(!isTRUE(failing$ok) && !dir.exists(file.path(Sys.getenv("CSL_RESULTS_ROOT"), "hubtest_fail")), "a failed start leaves nothing behind")
pass("start, status and resume")

# ---- deep links and one-time tickets ----------------------------------------------------------
target <- app$deep_link_target(paste0("?token=x&project=", res$project_id, "&view=results"))
assert(!is.null(target) && identical(target$tab, "Results Explorer"), "deep link opens the Results Explorer")
assert(is.null(app$deep_link_target("?project=rna/not_mine")), "deep links only open the user's own projects")
assert(identical(app$deep_link_target(paste0("?project=", res$project_id, "&view=../../x"))$tab, "Results Explorer"), "unknown views fall back safely")

tickets <- file.path(Sys.getenv("CSL_WEB_HOME"), "open_tickets")
dir.create(tickets, recursive = TRUE, showWarnings = FALSE)
ticket <- paste(sample(c(letters, 0:9), 40, TRUE), collapse = "")
writeLines(as.character(as.numeric(Sys.time()) + 60), file.path(tickets, ticket))
assert(app$consume_open_ticket(ticket, "page"), "a fresh ticket opens the page")
assert(app$consume_open_ticket(ticket, "session"), "...and its live session")
assert(!app$consume_open_ticket(ticket, "page") && !app$consume_open_ticket(ticket, "session"), "a ticket works only once")
old_ticket <- paste(sample(c(letters, 0:9), 40, TRUE), collapse = "")
writeLines(as.character(as.numeric(Sys.time()) - 1), file.path(tickets, old_ticket))
assert(!app$consume_open_ticket(old_ticket, "page"), "an expired ticket is refused")
assert(!app$consume_open_ticket("../../etc/passwd", "page"), "ticket names are validated")
pass("deep links and one-time tickets")

# ---- managed shared install ------------------------------------------------------------------
suite <- tempfile("suite_"); dir.create(file.path(suite, "CodeSpringLab", "scripts_DoNotTouch"), recursive = TRUE)
assert(!app$is_managed_install(file.path(suite, "CodeSpringLab")), "unmarked folders are not managed installs")
file.create(file.path(suite, "CODESPRING_MANAGED_INSTALL"))
assert(app$is_managed_install(file.path(suite, "CodeSpringLab")), "the marker file identifies a managed install")
unlink(suite, recursive = TRUE)
pass("managed install marker")

unlink(file.path(Sys.getenv("CSL_RESULTS_ROOT")), recursive = TRUE)
cat("All hub bridge tests passed.\n")
