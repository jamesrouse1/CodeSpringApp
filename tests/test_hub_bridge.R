# Tests for R/hub_bridge.R (codespring-run), deep links and one-time open tickets.
# Usage: Rscript tests/test_hub_bridge.R [/path/to/CodeSpringLab]
# Nothing is submitted to Slurm: start() is exercised with a fake controller starter.

args <- commandArgs(trailingOnly = TRUE)
repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), ".."), mustWork = TRUE)
lab_root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(file.path(repo_root, "..", "CodeSpringLab"), mustWork = TRUE)
test_home <- tempfile("hub_bridge_home_")
dir.create(test_home)
mail_log <- file.path(test_home, "mail.log")
fake_mailer <- file.path(test_home, "fake_mail.sh")
writeLines(c("#!/bin/sh", paste0("{ echo \"ARGS: $*\"; cat; echo; } >> '", mail_log, "'")), fake_mailer)
Sys.chmod(fake_mailer, "0755")
# A queue that is always empty, so the wave tests don't depend on the tester's own Slurm jobs.
fake_bin <- file.path(test_home, "bin")
dir.create(fake_bin)
writeLines(c("#!/bin/sh", "exit 0"), file.path(fake_bin, "squeue"))
Sys.chmod(file.path(fake_bin, "squeue"), "0755")
Sys.setenv(PATH = paste(fake_bin, Sys.getenv("PATH"), sep = ":"))
Sys.setenv(CSL_CODESPRINGLAB_ROOT = lab_root, CSL_WEB_APP_ROOT = repo_root, CSL_WEB_HOME = file.path(path.expand("~"), ".codespringweb_hubtest"),
           CSL_RESULTS_ROOT = file.path(test_home, "results"), CSL_HUB_POLL_SECONDS = "1", CSL_HUB_MAX_SUBMIT = "5",
           CSL_HUB_MAILER = fake_mailer, CSL_HUB_REPORT_FILE_MB = "0.01")
on.exit(unlink(c(test_home, Sys.getenv("CSL_WEB_HOME")), recursive = TRUE), add = TRUE)

app <- new.env(parent = globalenv())
old <- setwd(repo_root)
invisible(utils::capture.output(sys.source("app.R", envir = app)))
sys.source(file.path("R", "hub_bridge.R"), envir = app)
sys.source(file.path("R", "hub_report.R"), envir = app)
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
refused(base_spec(gsea = list(genesets = list("MSigDB_Hallmark_2020", "my_sets.gmt"))), "gsea.genesets", "every gene set must be listed")
refused(base_spec(gsea = list(genesets = as.list(head(app$GSEAPY_GENESET_OPTIONS, 7)))), "At most 6 gene sets", "gene sets are capped")
refused(base_spec(quantifiers = list("salmon")), "quantifiers", "only RSEM and Kallisto")
refused(base_spec(run_mode = "fast"), "run_mode", "only all or step")
pass("check refuses anything outside the BSR pipeline's options")

# ---- spec version 2: several gene sets, RSEM/Kallisto, one step at a time ------------------------
v2 <- app$hub_check_rnaseq(base_spec(project_name = "hubtest_v2", gsea = list(genesets = list("MSigDB_Hallmark_2020", "Reactome_Pathways_2024")),
                                     quantifiers = list("rsem", "kallisto"), run_mode = "step"))
assert(isTRUE(v2$ok), paste("v2 options accepted:", paste(v2$errors, collapse = "; ")))
plan2 <- app$hub_plan(v2)
assert(identical(vapply(plan2$steps, `[[`, "", "step"), c("Cutadapt", "FastQC", "STAR", "featureCounts", "DESeq2", "GSEA", "RSEM", "Kallisto")), "extra quantification runs at the end")
assert(plan2$steps[[6]]$jobs == 2L && grepl("Reactome", plan2$steps[[6]]$what), "one GSEA job per gene set")
assert(identical(plan2$run_mode, "step"), "plan shows the run mode")
none <- app$hub_check_rnaseq(base_spec(project_name = "hubtest_none", quantifiers = list("none"), gsea = list(geneset = "MSigDB_Hallmark_2020")))
assert(isTRUE(none$ok) && !length(none$quantifiers) && identical(app$hub_gsea_sets(none$gsea), "MSigDB_Hallmark_2020"), "'none' and the v1 single gene set still work")
pass("several gene sets, RSEM/Kallisto and run mode are checked and planned")

# ---- start (fake controller), status, resume ---------------------------------------------------
started_commands <- character(0)
starter <- function(command) { started_commands <<- c(started_commands, command); "4194399" }
res <- app$hub_start(base_spec(), starter = starter)
assert(isTRUE(res$ok) && identical(res$controller_pid, "4194399"), "start launches one controller")
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

# ---- one step at a time: the controller pauses after each step (fake submitters, no Slurm) -------
submitted <- character(0)
for (fn in c("submit_cutadapt_jobs", "submit_fastqc_jobs", "submit_star_jobs", "submit_featurecounts_jobs", "submit_deseq2_job",
             "submit_gseapy_job", "submit_rsem_jobs", "submit_kallisto_jobs")) {
  local({ name <- fn; assign(name, function(...) { submitted <<- c(submitted, name); "Submitted batch job 1" }, envir = app) })
}
real_run_step <- app$hub_run_step  # no waiting for Slurm here: each step just submits
assign("hub_run_step", function(p, name, entry, spec) { entry$run(app$project_samples(p)); list(ok = TRUE) }, envir = app)
stepped <- app$hub_start(base_spec(project_name = "hubtest_step", run_mode = "step", gsea = list(genesets = list("MSigDB_Hallmark_2020", "KEGG_2021_Human"))), starter = starter)
assert(isTRUE(stepped$ok), paste("step-mode run starts:", paste(stepped$errors, collapse = "; ")))
spec_path <- file.path(stepped$results_dir, "log", "codespring_hub_spec.json")
app$hub_drive(spec_path)
st <- app$hub_status(stepped$project_id)
assert(identical(st$state, "paused") && identical(st$paused_after, "Cutadapt") && identical(submitted, "submit_cutadapt_jobs"), "pauses after the first step")
refused_resume <- app$hub_resume(stepped$project_id, starter = starter)
assert(!isTRUE(refused_resume$ok), "a paused run is continued, not resumed")
cont <- app$hub_continue(stepped$project_id, starter = starter)
assert(isTRUE(cont$ok) && identical(cont$next_step, "FastQC"), "continue names the next step")
app$hub_drive(spec_path)
st <- app$hub_status(stepped$project_id)
assert(identical(st$paused_after, "FastQC") && identical(tail(submitted, 1), "submit_fastqc_jobs"), "...and runs only that step")
for (i in 1:4) { app$hub_continue(stepped$project_id, starter = starter); app$hub_drive(spec_path) }
st <- app$hub_status(stepped$project_id)
assert(identical(st$state, "complete"), "the last step finishes the run instead of pausing")
assert(sum(submitted == "submit_gseapy_job") == 2L, "GSEA submits one job per gene set")
assert(!isTRUE(app$hub_continue(stepped$project_id, starter = starter)$ok), "a finished run can't be continued")
assign("hub_run_step", real_run_step, envir = app)
pass("one step at a time: pause, continue, finish")

# ---- deep links and one-time tickets ----------------------------------------------------------
target <- app$deep_link_target(paste0("?token=x&project=", res$project_id, "&view=results"))
assert(!is.null(target) && identical(target$tab, "Results Explorer"), "deep link opens the Results Explorer")
assert(is.null(app$deep_link_target("?project=rna/not_mine")), "deep links only open the user's own projects")
assert(identical(app$deep_link_target(paste0("?project=", res$project_id, "&view=../../x"))$tab, "Results Explorer"), "unknown views fall back safely")

tickets <- file.path(Sys.getenv("CSL_WEB_HOME"), "open_tickets")
dir.create(tickets, recursive = TRUE, showWarnings = FALSE)
new_ticket <- function(ttl = 60) {
  t <- paste(sample(c(letters, LETTERS, 0:9, "_", "-"), 43, TRUE), collapse = "")
  writeLines(as.character(as.numeric(Sys.time()) + ttl), file.path(tickets, t))
  t
}
# Requests to the page as httpuv hands them to ui(): with a private access token configured.
app$ACCESS_TOKEN <- strrep("ab", 32)
page <- function(query = "", cookie = "") app$ui(list(QUERY_STRING = query, HTTP_COOKIE = cookie, PATH_INFO = "/"))
is_main <- function(x) identical(x, app$MAIN_UI)
is_denied <- function(x) !inherits(x, "httpResponse") && !is_main(x) && grepl("Access denied", as.character(x))

assert(is_main(page(paste0("?token=", app$ACCESS_TOKEN))), "the launcher's private URL still opens the app")
assert(is_denied(page("")) && is_denied(page("?token=wrong")), "no token, no access")

# A Hub link: http://localhost:PORT/?project=...&view=...#ticket=... The #fragment never reaches the
# server, so a prefetch only gets the opening page and the ticket is untouched.
ticket <- new_ticket()
first <- page(paste0("?project=", res$project_id, "&view=progress"))
assert(is_denied(first) && grepl("location.hash", as.character(first), fixed = TRUE), "a #ticket link first gets the page that moves the ticket into the address")
assert(file.exists(file.path(tickets, ticket)), "...and a prefetch doesn't use the ticket up")

# The page's script reloads with ?ticket=...: spent once, answered with a redirect that sets the key.
swap <- page(paste0("?project=", res$project_id, "&view=progress&ticket=", ticket))
assert(inherits(swap, "httpResponse") && swap$status == 302, "a fresh ticket is swapped (redirect)")
assert(identical(swap$headers$Location, paste0("/?project=", utils::URLencode(res$project_id, reserved = TRUE), "&view=progress")), "the redirect keeps project/view and drops the ticket")
set_cookie <- swap$headers$`Set-Cookie`
assert(grepl("HttpOnly", set_cookie) && grepl("SameSite=Lax", set_cookie) && startsWith(set_cookie, app$BROWSER_KEY_COOKIE), "the key is an HttpOnly cookie unique to this launch")
assert(!file.exists(file.path(tickets, ticket)) && !length(list.files(tickets, pattern = "spent")), "the ticket file is gone")
cookie <- sub(";.*$", "", set_cookie)
assert(is_main(page(paste0("?project=", res$project_id, "&view=progress"), cookie)), "with the key the page opens")
assert(is_main(page("", paste0("other=1; ", cookie))) && app$access_token_valid("", cookie), "reloads and live sessions (reconnects) work with the key")
assert(is_denied(page(paste0("?ticket=", ticket))), "a used ticket is refused in a fresh browser")
assert(is_denied(page("", paste0(app$BROWSER_KEY_COOKIE, "=", strrep("0", 64)))), "a made-up key is refused")
assert(is_denied(page("", paste0("csa_key_000000000000=", sub("^[^=]*=", "", cookie)))), "another launch's cookie name is refused")
assert(is_denied(page(paste0("?ticket=", new_ticket(-1)))), "an expired ticket is refused")
assert(is_denied(page("?ticket=../../etc/passwd")) && !app$consume_open_ticket("../../etc/passwd"), "ticket names are validated")
app$ACCESS_TOKEN <- ""
pass("deep links and one-time tickets")

# ---- managed shared install ------------------------------------------------------------------
suite <- tempfile("suite_"); dir.create(file.path(suite, "CodeSpringLab", "scripts_DoNotTouch"), recursive = TRUE)
assert(!app$is_managed_install(file.path(suite, "CodeSpringLab")), "unmarked folders are not managed installs")
file.create(file.path(suite, "CODESPRING_MANAGED_INSTALL"))
assert(app$is_managed_install(file.path(suite, "CodeSpringLab")), "the marker file identifies a managed install")
unlink(suite, recursive = TRUE)
pass("managed install marker")

# ---- optional steps, modes, notifications: planning ----------------------------------------------
cat2 <- app$hub_catalog()$rnaseq
assert(setequal(names(cat2$quantifiers), c("kallisto", "rsem")), "the catalog offers the app's own RSEM and Kallisto")
assert(identical(cat2$batching$max_jobs_at_once, 5L), "the Slurm job limit is read (here from CSL_HUB_MAX_SUBMIT)")
refused(base_spec(project_name = "hubtest_opt", notify = list(email = TRUE, to = "someone@else.org")), "own address", "emails can't go to another address")
refused(base_spec(project_name = "hubtest_opt", notify = list(email = TRUE, every_hours = 0.1)), "every_hours", "no email storms")
chk <- app$hub_check_rnaseq(base_spec(project_name = "hubtest_opt", quantifiers = list("rsem", "kallisto"), run_mode = "step",
                                      notify = list(email = TRUE, every_hours = 6)))
assert(isTRUE(chk$ok), paste(chk$errors, collapse = "; "))
pl <- app$hub_plan(chk)
assert(identical(pl$run_mode, "step") && identical(pl$notify$email_to, app$hub_email_address()), "plan shows the run mode and the person's own email")
assert(grepl("waves of at most 3", pl$batching$note), "plan explains the waves (limit 5 minus 2 kept free)")
pass("email choices and waves are checked and planned")

# ---- waves under the Slurm limit ---------------------------------------------------------------
waves <- list()
fake_entry <- list(run = function(s) { waves[[length(waves) + 1]] <<- s; "Submitted batch job 1" }, done = function() TRUE)
rs <- app$hub_run_step(p, "FastQC", fake_entry, list(notify = list(email = FALSE)))
assert(isTRUE(rs$ok), "a step completes")
assert(identical(lengths(waves), c(3L, 3L)) && setequal(unlist(waves), app$project_samples(p)), "six samples go out as two waves of three")
pass("per-sample jobs are submitted in waves")

# ---- pausing, continuing, adding and skipping steps ----------------------------------------------
raw <- app$hub_read_json(app$hub_status_path(p))
raw$state <- "paused"; raw$paused_after <- "FastQC"; raw$next_step <- "STAR"
for (i in seq_along(raw$steps)) raw$steps[[i]]$state <- if (raw$steps[[i]]$step %in% c("Cutadapt", "FastQC")) "complete" else "waiting"
app$hub_write_json(raw, app$hub_status_path(p))
st <- app$hub_status(res$project_id)
assert(identical(st$state, "paused") && any(grepl("continue", st$next_actions)), "a paused run says what can be done next")
bad <- app$hub_continue(res$project_id, add = "salmon")
assert(!isTRUE(bad$ok) && grepl("quantifiers must be", bad$errors), "can't add an unknown step")
bad <- app$hub_continue(res$project_id, skip = "STAR")
assert(!isTRUE(bad$ok) && grepl("can be skipped", bad$errors), "required steps can't be skipped")
n_started <- length(started_commands)
cont <- app$hub_continue(res$project_id, add = c("kallisto", "rsem"), run_mode = "all", starter = starter)
assert(isTRUE(cont$ok) && length(started_commands) == n_started + 1, "continuing relaunches the controller")
order_now <- vapply(cont$steps, `[[`, "", "step")
assert(identical(order_now, c("Cutadapt", "FastQC", "STAR", "featureCounts", "DESeq2", "RSEM", "Kallisto")), "added steps take their place in BSR's order (at the end)")
assert(identical(cont$next_step, "STAR"), "continue names the next step")
spec_now <- app$hub_read_json(app$hub_spec_path(p))
assert(identical(spec_now$run_mode, "all") && identical(unlist(spec_now$quantifiers), c("rsem", "kallisto")), "the stored request records the change")
again <- app$hub_continue(res$project_id, add = "kallisto")
assert(!isTRUE(again$ok) && grepl("Only a paused run", again$errors), "a running (queued) run can't be continued twice")
pass("pause, continue, add and skip")

# ---- cancelling ------------------------------------------------------------------------------------
cx <- app$hub_cancel(res$project_id)
assert(isTRUE(cx$ok) && identical(cx$state, "cancelled"), "a run can be cancelled")
assert(identical(app$hub_status(res$project_id)$state, "cancelled"), "status says cancelled")
assert(!isTRUE(app$hub_cancel(res$project_id)$ok), "cancelling twice is refused")
assert(any(grepl("resume", app$hub_status(res$project_id)$next_actions)), "a cancelled run can be resumed")
assert(isTRUE(app$hub_resume(res$project_id, starter = starter)$ok) && !file.exists(app$hub_cancel_flag(p)), "resuming clears the cancel request")
pass("cancel and resume")

# ---- QC read from each tool's own report --------------------------------------------------------
s1 <- app$project_samples(p)[[1]]
star_dir <- file.path(p$data_dir, "star", s1); dir.create(star_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(c("                          Number of input reads |\t1000000", "                   Uniquely mapped reads % |\t91.50%",
             "        % of reads mapped to multiple loci |\t4.20%", "                 % of reads unmapped: too short |\t3.10%"),
           file.path(star_dir, paste0(s1, "Log.final.out")))
fq_dir <- file.path(p$data_dir, "fastqc_cutadapt"); dir.create(fq_dir, recursive = TRUE, showWarnings = FALSE)
pairs <- app$sample_fastq_pairs(p, TRUE)
r1 <- basename(app$split_fastq_path_list(pairs$r1[pairs$sample == s1])[1])
inner_name <- sub(app$fastq_suffix_regex, "_fastqc", r1, ignore.case = TRUE)
zroot <- tempfile("fq_"); dir.create(file.path(zroot, inner_name), recursive = TRUE)
writeLines(c(paste("PASS", "Basic Statistics", r1, sep = "\t"), paste("FAIL", "Per base sequence content", r1, sep = "\t"),
             paste("WARN", "Sequence Duplication Levels", r1, sep = "\t")), file.path(zroot, inner_name, "summary.txt"))
old_wd <- setwd(zroot); utils::zip(file.path(fq_dir, paste0(inner_name, ".zip")), file.path(inner_name, "summary.txt"), flags = "-q"); setwd(old_wd)
qc <- app$hub_qc(res$project_id)
c1 <- Filter(function(r) identical(r$sample, s1), qc$samples)[[1]]
assert(identical(c1$star$uniquely_mapped_pct, 91.5) && identical(c1$star$input_reads, 1e6), "STAR numbers come from Log.final.out")
assert(identical(c1$fastqc$fail, 1L) && identical(as.character(c1$fastqc$failed_modules), "Per base sequence content"), "FastQC's own FAIL flags are read from its zip")
assert(any(grepl(s1, qc$summary$fastqc$samples_with_a_fail)) && grepl("91.5%", qc$summary$star$uniquely_mapped), "the summary names failing samples and ranges")
only_star <- app$hub_qc(res$project_id, "STAR")$samples
assert(all(vapply(only_star, function(r) is.null(r$fastqc), logical(1))) && !is.null(Filter(function(r) identical(r$sample, s1), only_star)[[1]]$star), "one step's QC can be asked for")
pass("QC from FastQC zips and STAR logs")

# ---- emails: only if asked, only to the person ---------------------------------------------------
unlink(mail_log)
app$hub_notify(p, list(notify = list(email = FALSE)), "test", "body")
assert(!file.exists(mail_log), "no email unless the person asked for it")
app$hub_notify(p, list(notify = list(email = TRUE, every_hours = 0)), "STAR finished", app$hub_status_text(p))
sent <- readLines(mail_log)
assert(any(grepl(app$hub_email_address(), sent, fixed = TRUE)) && any(grepl("STAR finished", sent)), "email goes to the person's own address")
assert(endsWith(app$hub_email_address(), "@cshl.edu"), "address is the person's CSHL account")
pass("optional email updates")

# ---- review fixes: no signals to reused pids, nothing from the status file reaches a shell ----------
spec_file <- app$hub_spec_path(p)
assert(!app$hub_controller_alive(list(controller_pid = as.character(Sys.getpid()), controller_host = Sys.info()[["nodename"]]), spec_file),
       "a live process that isn't this run's controller is never treated as it (pid reuse)")
assert(!app$hub_controller_alive(list(controller_pid = "1; touch /tmp/x", controller_host = Sys.info()[["nodename"]]), spec_file), "non-numeric pids are refused")
raw <- app$hub_read_json(app$hub_status_path(p))
canary <- file.path(test_home, "injected")
raw$state <- "running"; raw$controller_pid <- paste0("1; touch ", canary); raw$controller_host <- Sys.info()[["nodename"]]
app$hub_write_json(raw, app$hub_status_path(p))
st <- app$hub_status(res$project_id)
assert(!file.exists(canary) && identical(st$state, "stopped"), "a tampered controller_pid runs nothing and the run shows as stopped")
cx <- app$hub_cancel(res$project_id)
assert(isTRUE(cx$ok) && identical(cx$controller, "not running") && is.null(app$hub_read_json(app$hub_status_path(p))$controller_pid),
       "cancel signals nothing it can't verify, and forgets the pid")
pass("cancel can't hit reused pids; status file can't inject commands")

# ---- the post-run report packet ------------------------------------------------------------------
de <- file.path(p$data_dir, "deseq2")
dir.create(de, recursive = TRUE, showWarnings = FALSE)
cmp <- "treated_vs_control(ref)"
writeLines(c("baseMean\tlog2FoldChange\tlfcSE\tpvalue\tpadj", "GeneUp\t100\t2\t0.1\t1e-10\t1e-8",
             "<script>alert(1)</script>\t50\t-1.5\t0.1\t1e-6\t1e-4"), file.path(de, paste0("DEG_fdr0.05_", cmp, ".txt")))
writeLines(c("baseMean\tlog2FoldChange\tlfcSE\tpvalue\tpadj", "GeneUp\t100\t2\t0.1\t1e-10\t1e-8"), file.path(de, paste0("DEG_", cmp, ".txt")))
writeLines(rep(paste(rep("123.456", 40), collapse = "\t"), 200), file.path(de, paste0("normalized_counts_", cmp, ".txt")))  # > the 10 KB test cap
grDevices::png(file.path(de, paste0("volcano_", cmp, ".png")), width = 200, height = 200); plot(1:3); invisible(grDevices::dev.off())
rep_ <- app$hub_report(res$project_id)
html <- paste(readLines(rep_$report, warn = FALSE), collapse = "\n")
assert(isTRUE(rep_$ok) && file.exists(rep_$report) && startsWith(rep_$report, app$project_result_dir(p)), "the report is written inside the run's own results folder")
assert(identical(rep_$deg_summary[[cmp]], list(significant = 2L, up = 1L, down = 1L)), "DEG counts are read from DESeq2's own FDR list")
assert(grepl("data:image/png;base64,", html, fixed = TRUE) && grepl("Download PNG", html, fixed = TRUE), "plots are embedded with a download link")
assert(grepl(paste0("download=\"DEG_fdr0.05_", cmp, ".txt\""), html, fixed = TRUE), "the DEG list is downloadable from the page")
assert(any(vapply(rep_$left_out, function(x) startsWith(x$label, "normalized_counts_"), logical(1))) && !grepl("download=\"normalized_counts_", html, fixed = TRUE),
       "files over the size cap are listed with their path, not embedded")
assert(!grepl("<script>alert", html, fixed = TRUE) && grepl("&lt;script&gt;alert", html, fixed = TRUE), "text from result files is HTML-escaped")
assert(!grepl("base64,[A-Za-z0-9+/=]*\n", html), "embedded data has no line breaks")
assert(grepl("hasn't finished", html, fixed = TRUE), "a report made mid-run says so")
assert(identical(app$hub_status(res$project_id)$report, rep_$report), "status points to the report")
pass("report packet: app plots and small files embedded, big ones listed, text escaped")

unlink(file.path(Sys.getenv("CSL_RESULTS_ROOT")), recursive = TRUE)
cat("All hub bridge tests passed.\n")
