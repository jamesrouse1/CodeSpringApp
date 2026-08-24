args <- commandArgs(trailingOnly = TRUE)
repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), ".."), mustWork = TRUE)
lab_root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(file.path(repo_root, "..", "CodeSpringLab-fix"), mustWork = TRUE)
Sys.setenv(CSL_CODESPRINGLAB_ROOT = lab_root)

app_env <- new.env(parent = globalenv())
sys.source(file.path(repo_root, "app.R"), envir = app_env)

assert <- function(value, message) if (!isTRUE(value)) stop("ASSERTION FAILED: ", message, call. = FALSE)

recursive_fastq_root <- tempfile("recursive_fastq_pool_")
dir.create(file.path(recursive_fastq_root, "batch_2025", "ECO-02"), recursive = TRUE)
dir.create(file.path(recursive_fastq_root, "batch_2026", "Sample_ECO-18", "fastq"), recursive = TRUE)
recursive_fastqs <- c(
  file.path(recursive_fastq_root, "batch_2025", "ECO-02", "ECO-02_batch2025_R1.fastq.gz"),
  file.path(recursive_fastq_root, "batch_2025", "ECO-02", "ECO-02_batch2025_R2.fastq.gz"),
  file.path(recursive_fastq_root, "batch_2026", "Sample_ECO-18", "fastq", "ECO-18_L007_001.R1.fastq.gz"),
  file.path(recursive_fastq_root, "batch_2026", "Sample_ECO-18", "fastq", "ECO-18_L007_001.R2.fastq.gz"),
  file.path(recursive_fastq_root, "batch_2026", "Sample_ECO-18", "fastq", "ECO-18_L008_001.R1.fastq.gz"),
  file.path(recursive_fastq_root, "batch_2026", "Sample_ECO-18", "fastq", "ECO-18_L008_001.R2.fastq.gz"),
  file.path(recursive_fastq_root, "batch_2025", "ECO-02", "ECO-02_batch2025_singletons.fastq.gz")
)
file.create(recursive_fastqs)
recursive_scan <- app_env$scan_fastq_dirs(
  c(file.path(recursive_fastq_root, "batch_2025"), file.path(recursive_fastq_root, "batch_2026")),
  paired = TRUE,
  metadata_cols = "treatment",
  infer_samples = TRUE
)
assert(NROW(recursive_scan) == 2L && all(grepl("^paired", recursive_scan$status)), "FASTQ discovery pairs reads recursively across multiple parent folders")
assert(all(startsWith(as.character(recursive_scan$filename), "/")), "recursive FASTQ discovery preserves absolute source paths")
assert(all(nzchar(as.character(recursive_scan$sample))), "recursive FASTQ discovery automatically infers sample identifiers")
assert(any(grepl(";", recursive_scan$filename, fixed = TRUE)), "recursive FASTQ discovery pools multiple lanes for one inferred sample")
assert(!any(grepl("singleton|orphan|unpaired", recursive_scan$filename, ignore.case = TRUE)), "automatic FASTQ discovery excludes singleton and orphan recovery files")
unlink(recursive_fastq_root, recursive = TRUE, force = TRUE)
runtime_files <- c(
  file.path(repo_root, "app.R"),
  file.path(repo_root, "run_codespringweb.sh"),
  list.files(file.path(lab_root, "scripts_DoNotTouch"), pattern = "\\.(R|r|py|sh)$", recursive = TRUE, full.names = TRUE)
)
runtime_text <- unlist(lapply(runtime_files[file.exists(runtime_files)], readLines, warn = FALSE), use.names = FALSE)
app_text <- paste(readLines(file.path(repo_root, "app.R"), warn = FALSE), collapse = "\n")
launcher_text <- paste(readLines(file.path(repo_root, "run_codespringweb.sh"), warn = FALSE), collapse = "\n")
server_source <- paste(deparse(body(app_env$MAIN_SERVER), width.cutoff = 500L), collapse = "\n")
cutrun_batch_status_source <- paste(deparse(body(app_env$cutrun_diffbind_batch_status_ui)), collapse = "\n")
owner_path_pattern <- "(/grid/bsr/home/rouse|/home/rouse|/Users/rouse|rouse@bamdev)"
assert(!any(grepl(owner_path_pattern, runtime_text, ignore.case = TRUE)), "runtime code contains no hardcoded rouse home, login, or server path")
scrna_step_meta <- app_env$run_step_meta(list(analysis_key = "scrna", analysis = "scRNA-seq", counts_only = FALSE))
assert(
  identical(as.character(scrna_step_meta$step), c("Input inspection", "QC & doublets", "Normalize & PCA", "UMAP & clustering", "Annotate & markers", "Signature scoring", "Differential expression", "Pathway analysis")) && NROW(scrna_step_meta) == 8L,
  "single-cell run-step metadata exposes the complete checkpointed analysis workflow"
)
assert(
  identical(app_env$scrna_stage_step("preprocess"), "Normalize & PCA") &&
    identical(app_env$scrna_stage_step("annotate"), "Annotate & markers"),
  "single-cell stage identifiers map to their visible workflow steps"
)
assert(
  !grepl("scrna_runtime_executable", app_text, fixed = TRUE) &&
    grepl("Annotation files are requested only at the annotation step", app_text, fixed = TRUE) &&
    grepl("Upload one Seurat/Scanpy object from laptop", app_text, fixed = TRUE),
  "single-cell setup exposes annotation and upload controls in the appropriate workflow steps"
)
assert(
  !app_env$scrna_uses_input_manifest(list(analysis_key = "scrna", scrna_input_mode = "single")) &&
    app_env$scrna_uses_input_manifest(list(analysis_key = "scrna", scrna_input_mode = "multiple")) &&
    grepl("Read a processed Seurat or Scanpy object", app_text, fixed = TRUE) &&
    grepl("One folder containing all samples", app_text, fixed = TRUE) &&
    grepl("new_scrna_inputs_editor_ui", app_text, fixed = TRUE) &&
    grepl("Use an existing sample-design file instead", app_text, fixed = TRUE),
  "scRNA setup separates processed objects from new single- or multi-sample analyses"
)
stale_scrna_manifest <- tempfile(fileext = ".tsv")
utils::write.table(
  data.frame(sample_id = paste0("sample_", seq_len(8)), input_path = paste0("/input/sample_", seq_len(8)), stringsAsFactors = FALSE),
  stale_scrna_manifest, sep = "\t", row.names = FALSE, quote = FALSE
)
assert(
  app_env$scrna_uses_input_manifest(list(analysis_key = "scrna", scrna_input_mode = "single", scrna_input_manifest = stale_scrna_manifest)),
  "an eight-row scRNA manifest overrides a stale one-input project flag"
)
unlink(stale_scrna_manifest)
assert(
  grepl("scrna_preintegration_umap_ui", app_text, fixed = TRUE) &&
    grepl("02_preintegration_umap_", app_text, fixed = TRUE) &&
    grepl("scrna_harmony_theta", app_text, fixed = TRUE) &&
    grepl("scrna_harmony_lambda", app_text, fixed = TRUE) &&
    grepl("scrna_harmony_max_iter", app_text, fixed = TRUE),
  "single-cell integration displays the uncorrected UMAP and records reproducible Harmony controls"
)
assert(
  grepl('NODE_HOST="$(hostname -s', launcher_text, fixed = TRUE) &&
    grepl("CSL_WEB_SSH_HOST", launcher_text, fixed = TRUE) &&
    !grepl("@bamdev1", launcher_text, fixed = TRUE),
  "launcher generates its tunnel command for the node that started the app rather than bamdev1"
)
assert(grepl("observeEvent(input$genome_browser_comparison", server_source, fixed = TRUE) && grepl("samples_override = available", server_source, fixed = TRUE), "changing a genome-browser comparison resets its samples and reloads the selected comparison")
assert(
  (
    grepl("eligible comparisons selected automatically", cutrun_batch_status_source, fixed = TRUE) ||
      grepl("eligible source × comparison jobs", cutrun_batch_status_source, fixed = TRUE)
  ) &&
    grepl("Submit all eligible comparisons", server_source, fixed = TRUE) &&
    grepl("as.character(plan$id[plan$eligible])", server_source, fixed = TRUE),
  "CUT&RUN DiffBind visibly selects and submits every eligible comparison"
)
assert(
  grepl("cutrun_diffbind_batch_status_ui(p, plan, jobs)", server_source, fixed = TRUE) &&
    !grepl("cutrun_diffbind_job_status", server_source, fixed = TRUE),
  "CUT&RUN DiffBind uses the RNA-style live comparison status summary"
)
assert(
  grepl("pca_differential_peaks.png", server_source, fixed = TRUE) &&
    grepl("PCA using differential peaks", server_source, fixed = TRUE),
  "CUT&RUN DiffBind Results Explorer exposes the contrast-specific differential-peak PCA"
)
assert(
  grepl("observeEvent(input$genome_browser_ready, {", server_source, fixed = TRUE) &&
    grepl("is.null(input$genome_browser_sample_peak_sample)", server_source, fixed = TRUE) &&
    grepl("send_genome_browser()", server_source, fixed = TRUE) &&
    grepl("genome_browser_mode_state", server_source, fixed = TRUE),
  "genome browser loads immediately after lightweight ATAC controls while preserving the selected browser mode"
)
assert(
  any(grepl("window.codespringIgvSignature === signature", runtime_text, fixed = TRUE)) &&
    any(grepl("signature = browser_signature", runtime_text, fixed = TRUE)),
  "genome browser reuses an unchanged IGV instance instead of repeatedly reloading the same large peak tracks"
)
igv_clear_position <- regexpr("if (previous && window.codespringIgvBrowser === previous)", app_text, fixed = TRUE)[[1]]
igv_stale_position <- regexpr("if (requestId !== window.codespringIgvLoadRequestId) return null;", substring(app_text, igv_clear_position), fixed = TRUE)[[1]]
assert(
  igv_clear_position > 0L && igv_stale_position > 0L,
  "superseded IGV reloads clear the removed browser before the next queued request, preventing a grey double-removal state"
)
browser_controls_source <- sub(
  "^[\\s\\S]*output\\$genome_browser_controls_ui <- renderUI\\(\\{",
  "",
  sub("\\n  observeEvent\\(input\\$genome_browser_mode,[\\s\\S]*$", "", app_text, perl = TRUE),
  perl = TRUE
)
assert(
  !grepl("progress_refresh()", browser_controls_source, fixed = TRUE) &&
    grepl("if (nzchar(remembered_mode)) remembered_mode else isolate(input$genome_browser_mode)", browser_controls_source, fixed = TRUE) &&
    grepl('freezeReactiveValue(input, "genome_browser_mode")', browser_controls_source, fixed = TRUE),
  "genome-browser controls avoid timer rebuilds and freeze programmatic mode changes"
)
assert(
  grepl('observeEvent(input$genome_browser_mode,', app_text, fixed = TRUE) &&
    grepl('}, ignoreInit = TRUE, priority = 1000)', app_text, fixed = TRUE),
  "only deliberate post-initialization browser-mode changes update the saved mode"
)
assert(
  grepl('send_genome_browser(', server_source, fixed = TRUE) &&
    grepl('mode_override = "comparison"', server_source, fixed = TRUE) &&
    grepl('send_genome_browser(mode_override = mode)', server_source, fixed = TRUE) &&
    grepl('}, ignoreInit = TRUE, priority = -100)', server_source, fixed = TRUE) &&
    grepl('isolate(tryCatch({', server_source, fixed = TRUE) &&
    grepl('Could not switch genome-browser tracks:', server_source, fixed = TRUE) &&
    grepl('if (is.null(input$genome_browser_show_differential_peaks)) TRUE', server_source, fixed = TRUE),
  "every browser-mode change explicitly replaces IGV tracks and comparison mode defaults to its differential BED"
)
assert(
  grepl("genome_browser_comparison_show_peaks", server_source, fixed = TRUE),
  "comparison views use a separate opt-in for potentially large individual sample peak tracks"
)
assert(
  grepl('c("Individual sample peaks" = "sample_peak")', app_text, fixed = TRUE) &&
    grepl('"genome_browser_sample_peak_sample", "Sample"', app_text, fixed = TRUE) &&
    grepl('"genome_browser_sample_peak_file", "Peak file"', app_text, fixed = TRUE) &&
    grepl('selectInput("genome_browser_comparison", "Differential comparison"', app_text, fixed = TRUE),
  "ATAC genome browser exposes sample-peak and completed-comparison selectors"
)
assert(
  grepl("new_counts_source_mode", server_source, fixed = TRUE) &&
    grepl("browse_new_counts_server_file", server_source, fixed = TRUE) &&
    grepl("open_server_browser(\"new_counts_server_file\", \"file\"", server_source, fixed = TRUE),
  "counts-only projects support laptop uploads and server-side file browsing"
)
assert(
  grepl("CODE_SPRING_UPLOAD_LIMIT_BYTES <- 2 * 1024^3", app_text, fixed = TRUE) &&
    grepl("shiny.maxRequestSize", app_text, fixed = TRUE),
  "large local uploads, including Seurat references, support files up to 2 GB"
)
assert(
  grepl("Seurat reference label transfer", app_text, fixed = TRUE) &&
    grepl("scrna_reference_file", app_text, fixed = TRUE) &&
    grepl("Inspect reference labels", app_text, fixed = TRUE) &&
    grepl("scrna_reference_label_selector_ui", app_text, fixed = TRUE) &&
    grepl("c(\"rda\", \"rds\")", app_text, fixed = TRUE),
  "Seurat projects inspect server-side or uploaded references and select a valid label source"
)
reference_choice_file <- tempfile(fileext = ".tsv")
utils::write.table(
  data.frame(value = c("", "cell_type"), source = c("Active identities", "cell_type"), label_count = c(4, 4), non_missing_cells = c(100, 100), stringsAsFactors = FALSE),
  reference_choice_file, sep = "\t", row.names = FALSE, quote = FALSE
)
reference_choices <- app_env$read_scrna_reference_label_choices(reference_choice_file)
assert(NROW(reference_choices) == 2L && identical(reference_choices$value, c("", "cell_type")), "reference label-choice files preserve active identities and valid metadata fields")
unlink(reference_choice_file)
assert(
  grepl("Use selected file", server_source, fixed = TRUE) &&
    grepl("file.access(value, mode = 4)", server_source, fixed = TRUE),
  "the server browser selects readable files rather than only folders"
)
assert(
  grepl("codespring-igv-locus", server_source, fixed = TRUE),
  "CUT&RUN peak selection navigates the established IGV instance"
)
assert(any(grepl("codespringIgvLoadPromise", runtime_text, fixed = TRUE)), "IGV replacements are serialized so repeated reload events cannot create duplicate browsers")
assert(grepl("comparison_default_locus", server_source, fixed = TRUE) && grepl("locus_override = top_peak", server_source, fixed = TRUE), "each differential comparison defaults IGV to its most significant ranked peak")
assert(app_env$path_is_within(app_env$APP_HOME, app_env$CURRENT_HOME), "private app state is derived from the effective Unix user's home")
assert(identical(app_env$DEFAULT_RESULTS_ROOT, normalizePath(file.path(app_env$CURRENT_HOME, "csl_results"), winslash = "/", mustWork = FALSE)), "default results root is derived from the effective Unix user's home")
blank_new_project <- app_env$new_project_from_inputs(list(
  new_project_analysis = "RNA-seq",
  new_project_name = "blank_setup_test",
  new_project_mode = "new",
  new_results_root = "",
  new_design_matrix_path = "",
  new_fastq_location_mode = "one",
  new_fastq_dir = "",
  new_paired_end = "paired"
))
assert(identical(blank_new_project$results_root, app_env$DEFAULT_RESULTS_ROOT), "a blank results field safely falls back to the current user's results root")
assert(!app_env$is_bundled_example_design(blank_new_project$design_matrix_path), "a normal new project does not inherit a bundled example design matrix")
assert(!nzchar(blank_new_project$fastq_dir), "a normal new project does not inherit bundled example FASTQs")
assert(
  grepl('default_fastq_dir <- ""', app_text, fixed = TRUE) &&
    grepl('default_design_dir <- ""', app_text, fixed = TRUE) &&
    grepl('updateTextInput(session, "new_results_root", value = DEFAULT_RESULTS_ROOT)', app_text, fixed = TRUE),
  "bundled example paths are loaded only by the explicit example-data action"
)
assert(
  grepl("setup-path-preview", app_text, fixed = TRUE) &&
    grepl('`data-path-source` = "new_design_matrix_path"', app_text, fixed = TRUE) &&
    grepl('`data-path-source` = "new_results_root"', app_text, fixed = TRUE) &&
    grepl("setup-sidebar-wide", app_text, fixed = TRUE) &&
    grepl("cslSyncFastqPathUi", app_text, fixed = TRUE) &&
    grepl("var setupTabActive", app_text, fixed = TRUE) &&
    grepl("showWideSidebar = setupTabActive && $('#new_project_ui:visible').length > 0", app_text, fixed = TRUE),
  "new-project setup exposes scrollable full-path previews and widens the sidebar only on the Setup tab"
)
assert(
  grepl('textInput("new_results_root", "Results output storage location", value = DEFAULT_RESULTS_ROOT', app_text, fixed = TRUE) &&
    grepl('choices = c("Seurat" = "seurat", "Scanpy" = "scanpy")', app_text, fixed = TRUE) &&
    grepl('updateSelectInput(session, "new_scrna_engine", selected = "seurat")', app_text, fixed = TRUE),
  "single-cell setup visibly defaults to the current user's results folder and the PBMC example defaults to Seurat"
)
pbmc_defaults <- app_env$scrna_pbmc3k_tutorial_settings()
assert(
  identical(pbmc_defaults$max_features, 2500L) &&
    identical(pbmc_defaults$doublet_method, "auto") &&
    isTRUE(pbmc_defaults$remove_doublets) &&
    grepl('numericInput("scrna_max_features", "Maximum detected genes per cell (0 = disabled)"', app_text, fixed = TRUE) &&
    !grepl('tags$summary("Advanced QC setting")', app_text, fixed = TRUE) &&
    grepl('qc_preset = if (pbmc_example) "pbmc3k" else ""', app_text, fixed = TRUE) &&
    grepl('tags$strong("Seurat PBMC 3K tutorial cutoffs")', app_text, fixed = TRUE),
  "the PBMC example exposes its maximum-gene filter directly and removes predicted doublets by default"
)
pbmc_manifest_dir <- tempfile("pbmc_example_")
dir.create(file.path(pbmc_manifest_dir, "manifest"), recursive = TRUE)
pbmc_manifest <- file.path(pbmc_manifest_dir, "manifest", "scrna_samples.tsv")
utils::write.table(data.frame(sample_id = "pbmc3k", input_path = app_env$SCRNA_PBMC3K_EXAMPLE_DIR), pbmc_manifest, sep = "\t", row.names = FALSE, quote = FALSE)
pbmc_project <- list(data_dir = pbmc_manifest_dir, design_matrix_path = pbmc_manifest, scrna_input_manifest = pbmc_manifest)
pbmc_local_umap <- app_env$scrna_umap_focus_settings(pbmc_project, "local")
pbmc_global_umap <- app_env$scrna_umap_focus_settings(pbmc_project, "global")
assert(
  identical(pbmc_local_umap, list(n_neighbors = 20L, min_dist = 0.3, n_pcs = 10L)) &&
    identical(pbmc_global_umap, list(n_neighbors = 30L, min_dist = 0.6, n_pcs = 10L)) &&
    grepl('if (is.null(tutorial)) fluidRow(', app_text, fixed = TRUE) &&
    grepl('tutorial_umap$n_neighbors %||% input$scrna_n_neighbors', app_text, fixed = TRUE),
  "the PBMC example hides manual embedding fields and applies distinct Local/Global presets"
)
dir.create(file.path(pbmc_manifest_dir, "scrna", "tables"), recursive = TRUE)
utils::write.table(
  data.frame(cluster = c("0", "1"), `Naive CD4 T` = c(0.8, 0.1), `B cells` = c(0.2, 0.9), check.names = FALSE),
  file.path(pbmc_manifest_dir, "scrna", "tables", "marker_annotation_cluster_scores.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
suggestions <- app_env$scrna_marker_suggestion_table(pbmc_project)
assert(
  identical(suggestions$top_suggestion, c("Naive CD4 T", "B cells")) &&
    grepl("Automatically use the highest normalized average per cluster", app_text, fixed = TRUE) &&
    grepl("Review suggestions and label clusters manually", app_text, fixed = TRUE) &&
    grepl('marker_species = "auto"', app_text, fixed = TRUE),
  "the PBMC example supports automatic species detection and ranked automatic/manual cluster annotation"
)
utils::write.table(
  data.frame(cell = c("cell_A", "cell_B"), UMAP_1 = c(0, 1), UMAP_2 = c(1, 0), cluster = c("0", "1")),
  file.path(pbmc_manifest_dir, "scrna", "tables", "umap_coordinates.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
utils::write.table(
  data.frame(cell = c("cell_A", "cell_B"), cell_type = c("Naive CD4 T", "B cells")),
  file.path(pbmc_manifest_dir, "scrna", "tables", "cell_metadata.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
assert(
  "cell_type" %in% app_env$scrna_all_metadata_columns(pbmc_project) &&
    "cell_type" %in% app_env$scrna_de_population_fields(pbmc_project, app_env$scrna_all_metadata_columns(pbmc_project)) &&
    grepl('signature_species = "auto"', app_text, fixed = TRUE) &&
    !grepl('selectInput("scrna_signature_species"', app_text, fixed = TRUE) &&
    grepl('Group marker figures by', app_text, fixed = TRUE) &&
    grepl('(_panel_[0-9]+)?\\\\.png', app_text, fixed = TRUE),
  "completed annotations feed cell-type differential expression and signature species is detected automatically"
)
manual_map <- app_env$write_scrna_manual_cluster_mapping(pbmc_project, "cell_type", c(`0` = "Naive CD4 T", `1` = "B cells"))
manual_values <- utils::read.delim(manual_map, stringsAsFactors = FALSE)
assert(identical(manual_values$cell_type, c("Naive CD4 T", "B cells")), "manual PBMC cluster choices are expanded to a pipeline-ready per-cell mapping")
de_dir <- file.path(pbmc_manifest_dir, "scrna", "differential_expression", "B_cells_vs_Naive_CD4_T", "global")
dir.create(de_dir, recursive = TRUE)
writeLines("gene\tavg_log2FC\tp_val\nCD79A\t1.2\t0.001", file.path(de_dir, "cell_level_Wilcoxon__B_cells_vs_Naive_CD4_T__global.tsv"))
de_catalog <- app_env$scrna_completed_de_comparison_catalog(pbmc_project)
assert(
  NROW(de_catalog) == 1L && identical(de_catalog$method[[1]], "Cell-level Wilcoxon") &&
    grepl('selectInput("scrna_pathway_de_file", "Differential-expression comparison"', app_text, fixed = TRUE),
  "pathway analysis lists every completed cell-level or pseudobulk comparison"
)
unlink(pbmc_manifest_dir, recursive = TRUE)
assert(
  grepl("rna_overview_sample_progress_ui", server_source, fixed = TRUE) &&
    grepl('"Results Explorer"', server_source, fixed = TRUE),
  "RNA-seq Results Explorer reuses the live Progress sample matrix and refreshes it while visible"
)
assert(
  grepl("observeEvent(input$refresh_rna_results", server_source, fixed = TRUE) &&
    grepl('native_registered_id("")', server_source, fixed = TRUE),
  "RNA-seq Results Explorer can rescan newly completed DESeq2 outputs without restarting the app"
)
assert(
  grepl("Automatic (Harmony when a technical batch is selected)", app_text, fixed = TRUE) &&
    grepl("RPCA (anchor-based; smaller datasets)", app_text, fixed = TRUE),
  "Seurat automatic integration uses scalable Harmony while retaining explicit anchor methods for smaller datasets"
)
embedding_test_root <- tempfile("scrna_embedding_views_")
dir.create(file.path(embedding_test_root, "scrna", "tables"), recursive = TRUE)
embedding_project <- list(data_dir = embedding_test_root, analysis_key = "scrna", analysis = "scRNA-seq")
utils::write.table(
  data.frame(cell = c("c1", "c2"), UMAP_1 = c(10, 20), UMAP_2 = c(30, 40), cluster = c("0", "1"), condition = c("control", "treated")),
  file.path(embedding_test_root, "scrna", "tables", "umap_coordinates.tsv"), sep = "\t", row.names = FALSE, quote = FALSE
)
utils::write.table(
  data.frame(cell = c("c1", "c2"), UMAP_1 = c(-1, -2), UMAP_2 = c(-3, -4), sample_id = c("s1", "s2")),
  file.path(embedding_test_root, "scrna", "tables", "preintegration_umap_coordinates.tsv"), sep = "\t", row.names = FALSE, quote = FALSE
)
utils::write.table(
  data.frame(cell = c("c1", "c2"), cell_type = c("Monocyte", "HSC"), cell_type_prediction_score = c(0.91, 0.87)),
  file.path(embedding_test_root, "scrna", "tables", "reference_transfer_per_cell__cell_type.tsv"), sep = "\t", row.names = FALSE, quote = FALSE
)
embedding_views <- app_env$scrna_embedding_view_choices(embedding_project)
unintegrated_embedding <- app_env$scrna_embedding_table(embedding_project, columns = c("cluster", "condition", "sample_id", "cell_type"), max_points = Inf, view = "unintegrated")
assert(identical(unname(embedding_views), c("integrated", "unintegrated")), "interactive UMAP offers integrated and unintegrated coordinates when both tables exist")
assert(grepl('itemsizing = "constant"', app_text, fixed = TRUE) && grepl('legend = list(x = 1.02', app_text, fixed = TRUE), "interactive UMAP keeps categorical legend symbols visible and reserves a fixed legend column")
assert(grepl("After integration / final clustering UMAP", app_text, fixed = TRUE), "run section labels the final post-integration UMAP")
assert(identical(unintegrated_embedding$UMAP_1, c(-1, -2)) && identical(as.character(unintegrated_embedding$cluster), c("0", "1")), "unintegrated UMAP retains its coordinates and joins final annotations by exact cell ID")
assert(identical(as.character(unintegrated_embedding$cell_type), c("Monocyte", "HSC")) && "cell_type" %in% app_env$scrna_embedding_color_choices(embedding_project, "unintegrated"), "interactive UMAP discovers completed reference-transfer labels without requiring annotation to be rerun")
unlink(file.path(embedding_test_root, "scrna", "tables", "umap_coordinates.tsv"))
assert(identical(app_env$scrna_selected_embedding_view(embedding_project, ""), "unintegrated"), "interactive UMAP displays the unintegrated coordinates while integration is still running")
assert("sample_id" %in% app_env$scrna_embedding_color_choices(embedding_project, "unintegrated"), "the running integration UMAP can be colored by sample")
unlink(embedding_test_root, recursive = TRUE, force = TRUE)
assert(
  grepl('"Each input sample (recommended default)" = "sample_id"', app_text, fixed = TRUE) &&
    grepl('batch_column = if (is.null(input$scrna_batch_column)) "sample_id"', app_text, fixed = TRUE),
  "multi-sample integration defaults to one integration group per input sample"
)
assert(identical(unname(app_env$project_analysis_choices()), c("RNA-seq", "scRNA-seq", "ATAC-seq", "CUT&RUN", "ChIP-seq")), "biological project selectors use one canonical order and spelling")
assert(
  grepl('unset = file.path(DEFAULT_RESULTS_ROOT, "fetchngs")', app_text, fixed = TRUE) &&
    grepl('unset = file.path(CURRENT_HOME, ".codespringflow")', app_text, fixed = TRUE) &&
    !identical(app_env$FETCHNGS_RESULTS_ROOT, app_env$FETCHNGS_RUNTIME_ROOT),
  "FetchNGS results use the user's csl_results folder while cache and work data use a separate private runtime folder"
)
main_analysis_choices <- unname(app_env$analysis_choices())
project_creation_choices <- unname(app_env$project_analysis_choices())
established_project_choices <- c("RNA-seq", "scRNA-seq", "ATAC-seq", "CUT&RUN", "ChIP-seq")
assert(
  sum(main_analysis_choices == "FetchNGS") == 1L &&
    !"FetchNGS" %in% project_creation_choices &&
    all(established_project_choices %in% main_analysis_choices) &&
    all(established_project_choices %in% project_creation_choices) &&
    identical(app_env$analysis_key("FetchNGS"), "fetchngs") &&
    identical(app_env$analysis_label("fetchngs"), "FetchNGS"),
  "FetchNGS is a distinct main analysis choice but is excluded from biological project creation"
)
assert(
  grepl("input.analysis == 'FetchNGS'", app_text, fixed = TRUE) &&
    grepl('fetchngs_tabs <- c("FetchNGS", "FetchNGS Outputs")', server_source, fixed = TRUE) &&
      grepl('analysis_tabs <- c("Setup", "Samples & Design", "Run Pipeline", "Progress", "Results Explorer", "Logs", "Methods")', app_text, fixed = TRUE) &&
    grepl('updateTabsetPanel(session, "web_main_tabs", selected = "FetchNGS")', server_source, fixed = TRUE) &&
    grepl('updateTabsetPanel(session, "web_main_tabs", selected = "Setup")', server_source, fixed = TRUE),
  "FetchNGS selection exposes its run and output tabs while returning to an analysis restores the project tabs"
)
show_fetch_tabs_pos <- regexpr('lapply(fetchngs_tabs, function(tab) showTab("web_main_tabs", tab, session = session))', server_source, fixed = TRUE)[[1]]
hide_analysis_tabs_pos <- regexpr('lapply(analysis_tabs, function(tab) hideTab("web_main_tabs", tab, session = session))', server_source, fixed = TRUE)[[1]]
show_analysis_tabs_pos <- regexpr('lapply(non_viewer_tabs, function(tab) showTab("web_main_tabs", tab, session = session))', server_source, fixed = TRUE)[[1]]
hide_fetch_tabs_pos <- regexpr('lapply(fetchngs_tabs, function(tab) hideTab("web_main_tabs", tab, session = session))', server_source, fixed = TRUE)[[1]]
assert(
  all(c(show_fetch_tabs_pos, hide_analysis_tabs_pos, show_analysis_tabs_pos, hide_fetch_tabs_pos) > 0L) &&
    show_fetch_tabs_pos < hide_analysis_tabs_pos &&
    show_analysis_tabs_pos < hide_fetch_tabs_pos,
  "analysis switching reveals and selects a destination before hiding the active source tab"
)
for (key in c("rna", "atac", "cutrun", "chip")) {
  tabs <- app_env$results_explorer_tabs(key)
  assert(identical(tabs[[1]], "Overview") && identical(tail(tabs, 1), "Files") && "QC" %in% tabs, paste(key, "follows the shared Results Explorer navigation contract"))
  assert(nzchar(app_env$analysis_description(key)), paste(key, "has a setup description"))
}
assert(app_env$is_codespring_process_command("Rscript -e shiny::runApp('/home/user/CodeSpringWeb', port=8601)"), "CodeSpringApp process command recognized")
assert(!app_env$is_codespring_process_command("Rscript unrelated_analysis.R"), "unrelated Rscript process is not treated as CodeSpringApp")
assert(!app_env$is_codespring_process_command("Rscript -e shiny::runApp('/home/user/another_app')"), "unrelated Shiny app is not treated as CodeSpringApp")
root <- tempfile("codespring-app-smoke-")
dir.create(root, recursive = TRUE)
on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

custom_fetchngs_root <- file.path(root, "custom_fetchngs")
assert(
  identical(
    app_env$resolve_fetchngs_results_root("default", "", app_env$FETCHNGS_RESULTS_ROOT),
    normalizePath(app_env$FETCHNGS_RESULTS_ROOT, winslash = "/", mustWork = FALSE)
  ) &&
    identical(
      app_env$resolve_fetchngs_results_root("custom", custom_fetchngs_root, app_env$FETCHNGS_RESULTS_ROOT),
      normalizePath(custom_fetchngs_root, winslash = "/", mustWork = FALSE)
    ) &&
    identical(app_env$ensure_fetchngs_results_root(custom_fetchngs_root), normalizePath(custom_fetchngs_root, winslash = "/", mustWork = TRUE)) &&
    dir.exists(custom_fetchngs_root),
  "FetchNGS keeps the home-folder default and creates an explicitly selected writable custom root"
)
relative_root_error <- tryCatch({ app_env$resolve_fetchngs_results_root("custom", "relative/fetchngs"); "" }, error = conditionMessage)
filesystem_root_error <- tryCatch({ app_env$resolve_fetchngs_results_root("custom", "/"); "" }, error = conditionMessage)
assert(
  grepl("absolute server path", relative_root_error, fixed = TRUE) &&
    grepl("filesystem root", filesystem_root_error, fixed = TRUE),
  "FetchNGS rejects relative custom roots and the filesystem root"
)

fetchngs_accessions <- app_env$parse_fetchngs_accessions("SRR14593545\nERR1160846, SRR14593545")
assert(
  identical(fetchngs_accessions, c("SRR14593545", "ERR1160846")),
  "FetchNGS pasted accessions are split and deduplicated"
)
accepted_fetchngs_accessions <- c(
  "SRR11605097", "ERR4007730", "DRR171822",
  "SRX8171613", "ERX4009132", "DRX162434",
  "SRS6531847", "ERS4399630", "DRS090921",
  "SAMN14689442", "SAMEA6638373", "SAMD00114846",
  "SRP256957", "ERP120836", "DRP004793",
  "SRA1068758", "ERA2420837", "DRA008156",
  "PRJNA625551", "PRJEB37513", "PRJDB4176",
  "GSM4432381", "GSE147507"
)
accepted_fetchngs_types <- vapply(accepted_fetchngs_accessions, app_env$fetchngs_accession_type, character(1))
invalid_fetchngs_accession_error <- tryCatch({ app_env$parse_fetchngs_accessions("NOT_AN_ACCESSION"); "" }, error = conditionMessage)
assert(
  all(nzchar(accepted_fetchngs_types)) &&
    identical(app_env$parse_fetchngs_accessions("srr14593545"), "SRR14593545") &&
    grepl("Unrecognized accession format", invalid_fetchngs_accession_error, fixed = TRUE),
  "FetchNGS recognizes only documented accession families and normalizes them to uppercase"
)
too_many_run_accessions <- paste0("SRR", seq_len(app_env$FETCHNGS_MAX_RUNS + 1L))
too_many_run_error <- tryCatch({ app_env$fetchngs_validate_run_accession_limit(too_many_run_accessions); "" }, error = conditionMessage)
assert(
  grepl(paste0("limit is ", app_env$FETCHNGS_MAX_RUNS), too_many_run_error, fixed = TRUE),
  "FetchNGS rejects more direct run accessions than the administrator limit"
)
fetchngs_pasted_input <- app_env$write_fetchngs_accession_input(
  "pasted_app_helper",
  paste(fetchngs_accessions, collapse = "\n"),
  app_home = root
)
assert(
  identical(tolower(tools::file_ext(fetchngs_pasted_input)), "csv") &&
    identical(readLines(fetchngs_pasted_input, warn = FALSE), fetchngs_accessions),
  "Shiny writes pasted FetchNGS accessions one per line with a validator-compatible CSV extension"
)
fetchngs_input <- file.path(root, "fetchngs_accessions.txt")
writeLines(fetchngs_accessions, fetchngs_input)
fetchngs_root <- file.path(root, "csl_results", "fetchngs")
fetchngs_runtime_root <- file.path(root, ".codespringflow")
fetchngs_output <- app_env$run_fetchngs_cli(
  c("fetchngs", "--input", fetchngs_input, "--name", "app_helper_smoke", "--metadata-only", "--dry-run"),
  results_root = fetchngs_root,
  runtime_root = fetchngs_runtime_root
)
fetchngs_run <- app_env$fetchngs_run_dir("app_helper_smoke", fetchngs_root)
fetchngs_summary <- app_env$fetchngs_run_summary("app_helper_smoke", fetchngs_root, query_scheduler = FALSE)
manifest_fixture_dir <- file.path(root, "fetchngs_manifest_fixture")
dir.create(manifest_fixture_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  c("field\tvalue", "run_name\tmanifest_fixture", "work_namespace\t", "metadata_only\ttrue"),
  file.path(manifest_fixture_dir, "run_manifest.tsv")
)
manifest_fixture <- app_env$fetchngs_read_manifest(manifest_fixture_dir)
manifestless_run <- app_env$fetchngs_run_dir("legacy_manifestless", fetchngs_root)
dir.create(manifestless_run, recursive = TRUE, showWarnings = FALSE)
manifestless_summary <- app_env$fetchngs_run_summary("legacy_manifestless", fetchngs_root, query_scheduler = FALSE)
fetchngs_history_path <- file.path(root, "fetchngs_results_roots.tsv")
app_env$remember_fetchngs_results_root(fetchngs_root, fetchngs_history_path)
known_fetchngs_roots <- app_env$fetchngs_known_results_roots(custom_fetchngs_root, fetchngs_history_path)
assert(
  grepl("Dry run only", fetchngs_output, fixed = TRUE) &&
    file.exists(file.path(fetchngs_run, "run.sbatch")) &&
    file.exists(file.path(fetchngs_run, "input", "accessions.csv")) &&
    identical(fetchngs_summary$Status[[1]], "Bundle only") &&
    identical(fetchngs_summary$`Metadata only`[[1]], "Yes") &&
    identical(app_env$fetchngs_manifest_value(manifest_fixture, "metadata_only", "false"), "true") &&
    identical(app_env$fetchngs_manifest_value(manifest_fixture, "work_namespace", "default"), "default") &&
    identical(app_env$fetchngs_manifest_value(manifest_fixture, "missing_field", "fallback"), "fallback") &&
    identical(manifestless_summary$Status[[1]], "Bundle only") &&
    identical(manifestless_summary$`Metadata only`[[1]], "No") &&
    identical(manifestless_summary$`FetchNGS version`[[1]], app_env$FETCHNGS_DEFAULT_VERSION),
  "Shiny FetchNGS helpers create a validator-compatible metadata-only bundle without submitting Slurm"
)
assert(
  all(normalizePath(c(custom_fetchngs_root, fetchngs_root), winslash = "/", mustWork = TRUE) %in% known_fetchngs_roots),
  "FetchNGS remembers prior readable output roots so older runs remain selectable"
)
fake_ena_reports <- list(
  SRR14593545 = data.frame(run_accession = "SRR14593545", fastq_bytes = "1000;2000", stringsAsFactors = FALSE),
  ERR1160846 = data.frame(run_accession = "ERR1160846", fastq_bytes = "4000", stringsAsFactors = FALSE)
)
fake_ena_report_fetcher <- function(accession, max_rows) fake_ena_reports[[accession]]
fake_shared_space <- function(path) list(filesystem = "testfs", available_bytes = 1000000, checked_path = path, mount = "/test")
fetchngs_preflight <- app_env$fetchngs_download_preflight(
  fetchngs_accessions,
  results_root = fetchngs_root,
  runtime_root = fetchngs_runtime_root,
  max_runs = 20L,
  max_download_bytes = 10000,
  storage_multiplier = 4,
  report_fetcher = fake_ena_report_fetcher,
  space_fetcher = fake_shared_space
)
assert(
  identical(fetchngs_preflight$run_count, 2L) &&
    identical(fetchngs_preflight$estimated_bytes, 7000) &&
    identical(fetchngs_preflight$required_storage_bytes, 28000),
  "FetchNGS preflight deduplicates ENA runs, sums paired FASTQ bytes, and applies the storage multiplier"
)
oversize_preflight_error <- tryCatch({
  app_env$fetchngs_download_preflight(
    fetchngs_accessions, fetchngs_root, fetchngs_runtime_root,
    max_runs = 20L, max_download_bytes = 6000, storage_multiplier = 4,
    report_fetcher = fake_ena_report_fetcher, space_fetcher = fake_shared_space
  )
  ""
}, error = conditionMessage)
missing_size_fetcher <- function(accession, max_rows) data.frame(run_accession = accession, fastq_bytes = "", stringsAsFactors = FALSE)
missing_size_preflight <- app_env$fetchngs_download_preflight(
  c("SRR14593545", "ERR1160846"), fetchngs_root, fetchngs_runtime_root,
  report_fetcher = function(accession, max_rows) {
    if (identical(accession, "SRR14593545")) return(missing_size_fetcher(accession, max_rows))
    fake_ena_reports[[accession]]
  },
  space_fetcher = fake_shared_space
)
missing_size_skip <- app_env$fetchngs_preflight_accessions(missing_size_preflight, "skip")
missing_size_continue <- app_env$fetchngs_preflight_accessions(missing_size_preflight, "continue")
all_unknown_skip_error <- tryCatch({
  all_unknown <- app_env$fetchngs_download_preflight(
    "SRR14593545", fetchngs_root, fetchngs_runtime_root,
    report_fetcher = missing_size_fetcher, space_fetcher = fake_shared_space
  )
  app_env$fetchngs_preflight_accessions(all_unknown, "skip")
  ""
}, error = conditionMessage)
lookup_failure_preflight <- app_env$fetchngs_download_preflight(
  "SRP256957", fetchngs_root, fetchngs_runtime_root,
  report_fetcher = function(accession, max_rows) stop("simulated ENA timeout"),
  space_fetcher = fake_shared_space
)
too_many_resolved_fetcher <- function(accession, max_rows) {
  data.frame(
    run_accession = paste0("SRR", seq_len(max_rows)),
    fastq_bytes = rep("100", max_rows),
    stringsAsFactors = FALSE
  )
}
too_many_resolved_error <- tryCatch({
  app_env$fetchngs_download_preflight(
    "SRP256957", fetchngs_root, fetchngs_runtime_root,
    max_runs = 20L, report_fetcher = too_many_resolved_fetcher, space_fetcher = fake_shared_space
  )
  ""
}, error = conditionMessage)
fake_low_space <- function(path) list(filesystem = "testfs", available_bytes = 20000, checked_path = path, mount = "/test")
low_space_preflight_error <- tryCatch({
  app_env$fetchngs_download_preflight(
    fetchngs_accessions, fetchngs_root, fetchngs_runtime_root,
    max_runs = 20L, max_download_bytes = 10000, storage_multiplier = 4,
    report_fetcher = fake_ena_report_fetcher, space_fetcher = fake_low_space
  )
  ""
}, error = conditionMessage)
assert(
  grepl("above the administrator limit", oversize_preflight_error, fixed = TRUE) &&
    grepl("more than the administrator limit", too_many_resolved_error, fixed = TRUE) &&
    grepl("Insufficient free storage", low_space_preflight_error, fixed = TRUE),
  "FetchNGS blocks oversized, over-20-run, and insufficient-storage downloads before submission"
)
assert(
  isTRUE(missing_size_preflight$has_unknown_size) &&
    identical(missing_size_preflight$unknown_runs, "SRR14593545") &&
    identical(missing_size_skip, "ERR1160846") &&
    identical(missing_size_continue, c("SRR14593545", "ERR1160846")) &&
    grepl("nothing left to submit", all_unknown_skip_error, fixed = TRUE) &&
    isTRUE(lookup_failure_preflight$has_unknown_size) &&
    identical(lookup_failure_preflight$unknown_accessions, "SRP256957"),
  "FetchNGS reports unknown sizes and supports explicit skip or acknowledged-continue decisions"
)
fetchngs_bundle_input <- file.path(fetchngs_run, "input", "accessions.csv")
fetchngs_bundle_backup <- app_env$replace_fetchngs_run_accession_input(
  fetchngs_bundle_input, missing_size_skip, fetchngs_run
)
assert(
  file.exists(fetchngs_bundle_backup) &&
    identical(readLines(fetchngs_bundle_backup, warn = FALSE), fetchngs_accessions) &&
    identical(readLines(fetchngs_bundle_input, warn = FALSE), missing_size_skip),
  "FetchNGS resume skip safely backs up and replaces the bundle accession list"
)
fetchngs_metadata_dir <- file.path(fetchngs_run, "results", "metadata")
fetchngs_fastq_dir <- file.path(fetchngs_run, "results", "fastq")
dir.create(fetchngs_metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fetchngs_fastq_dir, recursive = TRUE, showWarnings = FALSE)
fetchngs_metadata_file <- file.path(fetchngs_metadata_dir, "runinfo.tsv")
fetchngs_yaml_file <- file.path(fetchngs_run, "results", "software_versions.yml")
fetchngs_fastq_file <- file.path(fetchngs_fastq_dir, "sample.fastq.gz")
writeLines(c("run\tlibrary_layout", "SRR14593545\tPAIRED"), fetchngs_metadata_file)
writeLines(c("fetchngs: 1.12.0", "nextflow: 24.04.4"), fetchngs_yaml_file)
writeLines("compressed-placeholder", fetchngs_fastq_file)
fetchngs_catalog <- app_env$fetchngs_output_file_catalog("app_helper_smoke", fetchngs_root)
limited_fetchngs_files <- app_env$fetchngs_output_files("app_helper_smoke", fetchngs_root, max_files = 2L)
assert(
  NROW(fetchngs_catalog) == 3L &&
    all(c("File", "Folder", "Type", "Size", "Modified", "path") %in% names(fetchngs_catalog)) &&
    identical(app_env$validated_fetchngs_output_path("app_helper_smoke", fetchngs_metadata_file, fetchngs_root), normalizePath(fetchngs_metadata_file)) &&
    identical(app_env$validated_fetchngs_output_path("app_helper_smoke", fetchngs_input, fetchngs_root), "") &&
    identical(app_env$fetchngs_output_preview_kind(fetchngs_metadata_file), "table") &&
    identical(app_env$fetchngs_output_preview_kind(fetchngs_yaml_file), "text") &&
    identical(app_env$fetchngs_output_preview_kind(fetchngs_fastq_file), "binary"),
  "FetchNGS output viewing catalogs only selected-run result files and chooses safe preview modes"
)
assert(
  length(limited_fetchngs_files) == 2L && isTRUE(attr(limited_fetchngs_files, "truncated")),
  "FetchNGS output scanning is bounded so large result trees cannot block the app indefinitely"
)
fetchngs_log_dir <- file.path(fetchngs_run, "logs")
dir.create(fetchngs_log_dir, recursive = TRUE, showWarnings = FALSE)
writeLines("900002", file.path(fetchngs_run, "job_id.txt"))
writeLines("OLD FAILURE THAT MUST NOT BE SHOWN", file.path(fetchngs_log_dir, "controller-900001.err"))
waiting_fetchngs_log <- app_env$fetchngs_latest_log("app_helper_smoke", fetchngs_root)
waiting_fetchngs_summary <- app_env$fetchngs_run_summary("app_helper_smoke", fetchngs_root, query_scheduler = FALSE)
assert(
  grepl("SLURM job 900002 is the current submission", waiting_fetchngs_log, fixed = TRUE) &&
    !grepl("OLD FAILURE", waiting_fetchngs_log, fixed = TRUE) &&
    identical(waiting_fetchngs_summary$Status[[1]], "Submitted; waiting for Slurm"),
  "FetchNGS follows a newly recorded resume job while Slurm creates its current controller log"
)
writeLines("CURRENT RESUME LOG", file.path(fetchngs_log_dir, "controller-900002.out"))
current_fetchngs_log <- app_env$fetchngs_latest_log("app_helper_smoke", fetchngs_root)
assert(
  grepl("Current SLURM job: 900002", current_fetchngs_log, fixed = TRUE) &&
    grepl("CURRENT RESUME LOG", current_fetchngs_log, fixed = TRUE) &&
    !grepl("OLD FAILURE", current_fetchngs_log, fixed = TRUE),
  "FetchNGS displays only the current resumed job's controller logs"
)
assert(
  app_env$fetchngs_job_is_active("RUNNING") &&
    app_env$fetchngs_job_is_active("PENDING") &&
    !app_env$fetchngs_job_is_active("FAILED") &&
    !app_env$fetchngs_job_is_active("COMPLETED") &&
    app_env$fetchngs_job_is_terminal("FAILED") &&
    app_env$fetchngs_job_is_terminal("COMPLETED") &&
    !app_env$fetchngs_job_is_terminal("RUNNING") &&
    !app_env$fetchngs_job_is_terminal("UNKNOWN_TRANSITION"),
  "FetchNGS deletion distinguishes active jobs from terminal jobs"
)
# Simulate a successful resume replacing the previous job ID.
writeLines(c("111111", "222222"), file.path(fetchngs_run, "job_history.txt"))
writeLines("222222", file.path(fetchngs_run, "job_id.txt"))
dir.create(file.path(fetchngs_run, "logs"), recursive = TRUE, showWarnings = FALSE)
writeLines("resumed controller log", file.path(fetchngs_run, "logs", "controller-222222.out"))
original_fetchngs_scheduler_state <- app_env$fetchngs_scheduler_state
app_env$fetchngs_scheduler_state <- function(job_id) if (identical(job_id, "222222")) "RUNNING" else ""
resumed_fetchngs_summary <- app_env$fetchngs_run_summary("app_helper_smoke", fetchngs_root, query_scheduler = TRUE)
resumed_fetchngs_log <- app_env$fetchngs_latest_log("app_helper_smoke", fetchngs_root)
app_env$fetchngs_scheduler_state <- original_fetchngs_scheduler_state
assert(
  identical(resumed_fetchngs_summary$`Job ID`[[1]], "222222") &&
    identical(resumed_fetchngs_summary$Status[[1]], "Running") &&
    grepl("controller-222222.out", resumed_fetchngs_log, fixed = TRUE) &&
    grepl("resumed controller log", resumed_fetchngs_log, fixed = TRUE) &&
    isTRUE(app_env$fetchngs_job_submission_is_recent(fetchngs_run)),
  "FetchNGS resume synchronization replaces the old job ID and follows the new status and controller log"
)
Sys.setFileTime(file.path(fetchngs_run, "job_id.txt"), Sys.time() - 600)
assert(
  !app_env$fetchngs_job_submission_is_recent(fetchngs_run),
  "FetchNGS recent-submission protection expires after the scheduler registration window"
)

fetchngs_delete_message <- app_env$delete_fetchngs_run(
  "app_helper_smoke",
  results_root = fetchngs_root,
  runtime_root = fetchngs_runtime_root,
  app_home = root,
  query_scheduler = FALSE
)
assert(
  grepl("Deleted FetchNGS run", fetchngs_delete_message, fixed = TRUE) &&
    !dir.exists(fetchngs_run) &&
    !dir.exists(file.path(fetchngs_runtime_root, "work", "fetchngs", "app_helper_smoke")) &&
    dir.exists(fetchngs_root),
  "FetchNGS deletion removes only the selected result and work folders while keeping the results root"
)
initial_fetchngs_ui <- as.character(htmltools::renderTags(app_env$MAIN_UI)$html)
assert(
    grepl('tabPanel("FetchNGS"', app_text, fixed = TRUE) &&
    grepl('tabPanel("FetchNGS Outputs"', app_text, fixed = TRUE) &&
    grepl('textInput("fetchngs_older_results_root"', app_text, fixed = TRUE) &&
    grepl('uiOutput("fetchngs_output_root_ui")', app_text, fixed = TRUE) &&
    grepl('"fetchngs_results_mode", "FetchNGS results folder"', app_text, fixed = TRUE) &&
    grepl('actionButton("browse_fetchngs_results_root"', app_text, fixed = TRUE) &&
    grepl('actionButton("submit_fetchngs"', app_text, fixed = TRUE) &&
    grepl('uiOutput("fetchngs_run_actions_ui")', app_text, fixed = TRUE) &&
    grepl('selection = "single"', app_text, fixed = TRUE) &&
    grepl("fetchngs_runs_table_rows_selected", app_text, fixed = TRUE) &&
    grepl('updateSelectInput(session, "fetchngs_selected_run"', app_text, fixed = TRUE) &&
    grepl('"resume_fetchngs", "Resume selected run"', server_source, fixed = TRUE) &&
    grepl('"delete_fetchngs", "Delete selected run"', server_source, fixed = TRUE) &&
    grepl('completed <- identical(toupper(trimws(state)), "COMPLETED")', app_text, fixed = TRUE) &&
    grepl("resume_blocked <- active || recent || completed", app_text, fixed = TRUE) &&
    grepl("delete_blocked <- active || recent", app_text, fixed = TRUE) &&
    grepl('disabled = if (resume_blocked) "disabled" else NULL', app_text, fixed = TRUE) &&
    grepl('disabled = if (delete_blocked) "disabled" else NULL', app_text, fixed = TRUE) &&
    grepl("fetchngs_job_submission_is_recent(run_dir)", server_source, fixed = TRUE) &&
    grepl("fetchngs_run_actions_ui", initial_fetchngs_ui, fixed = TRUE) &&
    !grepl("resume_fetchngs", initial_fetchngs_ui, fixed = TRUE) &&
    !grepl("delete_fetchngs", initial_fetchngs_ui, fixed = TRUE) &&
    grepl('actionButton("confirm_delete_fetchngs"', server_source, fixed = TRUE) &&
    grepl('tags$strong("Accepted accessions")', app_text, fixed = TRUE) &&
    grepl('tags$strong("Download safeguards")', app_text, fixed = TRUE) &&
    grepl('uiOutput("fetchngs_message")', app_text, fixed = TRUE) &&
    grepl('actionButton("confirm_fetchngs_unknown_size"', server_source, fixed = TRUE) &&
    grepl('actionButton("cancel_fetchngs_unknown_size"', server_source, fixed = TRUE) &&
    grepl('fetchngs_unknown_size_ack', server_source, fixed = TRUE) &&
    grepl("fetchngs_download_preflight(accessions, results_root = results_root)", server_source, fixed = TRUE) &&
    grepl("Rechecking accession and download safety before resuming", server_source, fixed = TRUE) &&
    grepl("manifest <- fetchngs_read_manifest(run_dir)", server_source, fixed = TRUE) &&
    grepl("Resume submitted for", server_source, fixed = TRUE) &&
    grepl("job_id <- fetchngs_latest_job_id(context$run_dir)", server_source, fixed = TRUE) &&
    grepl("body:not(:has(#analysis option[value='FetchNGS']:checked))", app_text, fixed = TRUE) &&
    !grepl("codespring-tabs-initializing", app_text, fixed = TRUE) &&
    !grepl("codespring-tabs-ready", app_text, fixed = TRUE),
  "the app exposes FetchNGS location, creation, status, output viewing, resume, and confirmed deletion controls"
)

shared_scanpy_sif <- "/grid/bsr/data/data/bsr_readable_data/containers/scanpy/codespring-scanpy_1.0.0.sif"
assert(
  shared_scanpy_sif %in% app_env$scanpy_container_candidates(),
  "the shared BSR Scanpy image is a built-in runtime candidate for every user"
)

scanpy_sif <- file.path(root, "codespring-scanpy_1.0.0.sif")
file.create(scanpy_sif)
Sys.setenv(CSL_SCANPY_SIF = scanpy_sif)
scanpy_container <- app_env$scanpy_container_check()
assert(
  isTRUE(scanpy_container$ready) && identical(scanpy_container$path, normalizePath(scanpy_sif, winslash = "/", mustWork = FALSE)),
  "an explicit shared Scanpy SIF is detected without creating a per-user Python runtime"
)
Sys.unsetenv("CSL_SCANPY_SIF")

design_path <- file.path(root, "design_matrix.txt")
design <- data.frame(
  sample = c("A1", "I1", "A2", "I2", "B1", "I3", "B2", "I4"),
  treatment = rep(c("A", "A", "B", "B"), each = 2),
  reference = rep(c("chip", "input"), 4),
  condition = rep(c("A", "A", "B", "B"), each = 2),
  replicate = rep(c(1, 1, 2, 2), each = 2),
  control_sample = c("I1", "", "I2", "", "I3", "", "I4", ""),
  filename = paste0(c("A1", "I1", "A2", "I2", "B1", "I3", "B2", "I4"), ".fastq.gz"),
  stringsAsFactors = FALSE
)
write.table(design, design_path, sep = "\t", row.names = FALSE, quote = FALSE)

chip_project <- list(
  id = "fake-chip", name = "fake-chip", analysis_key = "chip", analysis = "ChIP-seq",
  design_matrix_path = design_path, data_dir = root, results_root = dirname(root),
  fastq_dir = root, fastq_dirs = root, paired_end = FALSE, genome = "mouse"
)
assert(identical(app_env$chip_control_sample_for(chip_project, "A1"), "I1"), "explicit ChIP control resolution")
assert(nrow(app_env$chip_target_design(chip_project)) == 4L, "input rows excluded from ChIP targets")
mouse_chip_ref <- app_env$chip_reference_resources(chip_project)
human_chip_project <- chip_project
human_chip_project$genome <- "human"
human_chip_ref <- app_env$chip_reference_resources(human_chip_project)
assert(identical(mouse_chip_ref$genome_version, "mouse_gencodeM39") && grepl("mouse_gencodeM39", mouse_chip_ref$bowtie2_index), "ChIP mouse reference uses GRCm39/GENCODE M39")
assert(identical(human_chip_ref$genome_version, "human_gencode50") && grepl("human_gencode50", human_chip_ref$bowtie2_index), "ChIP human reference uses GRCh38/GENCODE v50")
assert(length(app_env$genome_reference_choices("mouse", "ChIP-seq")) == 1L, "ChIP setup offers only the current mouse reference")
assert(length(app_env$genome_reference_choices("human", "ChIP-seq")) == 1L, "ChIP setup offers only the current human reference")
chip_example <- app_env$example_dataset_paths("chip")
assert(identical(chip_example$species, "human") && identical(chip_example$paired_end, "single"), "bundled human ChIP example selects GRCh38 and single-end reads")

dir.create(file.path(root, "fastqc"), recursive = TRUE, showWarnings = FALSE)
screen_path <- file.path(root, "fastqc", "A1_screen.txt")
write.table(
  data.frame(Genome = c("Human", "Mouse"), `%Unmapped` = c(2, 98), check.names = FALSE),
  screen_path, sep = "\t", row.names = FALSE, quote = FALSE
)
screen_pairs <- data.frame(sample = "A1", r1 = file.path(root, "A1.fastq.gz"), r2 = file.path(root, "A1.fastq.gz"), stringsAsFactors = FALSE)
assert(nzchar(app_env$fastq_screen_species_mismatch(chip_project, screen_pairs)), "ChIP preflight blocks a strong human-versus-mouse FastQ Screen mismatch")
assert(!nzchar(app_env$fastq_screen_species_mismatch(human_chip_project, screen_pairs)), "ChIP preflight accepts the matching human reference")
maize_choices <- app_env$genome_reference_choices("maize", "RNA-seq")
assert(
  identical(unname(maize_choices), c("maize_b73_nam5", "maize_nc350_nam1", "maize_w22_nrgene2")),
  "RNA-seq setup offers B73, NC350, and W22 maize references"
)
maize_rna_project <- within(chip_project, {
  analysis_key <- "rna"; analysis <- "RNA-seq"; genome <- "maize"; genome_version <- "maize_nc350_nam1"
})
maize_resources <- app_env$genome_resources(maize_rna_project)
assert(
  identical(app_env$genome_species(maize_rna_project), "maize") &&
    grepl("STAR_index/NC350$", maize_resources$star_index) &&
    grepl("Zm-NC350-REFERENCE-NAM-1.0.gtf$", maize_resources$gtf) &&
    grepl("Zm-NC350-REFERENCE-NAM-1.0.annotation_forStrandDetect_geneID.bed$", maize_resources$strand_bed),
  "maize RNA-seq resolves the selected variety's matching STAR index, GTF, and strand BED"
)
maize_reference_files <- list(
  maize_b73_nam5 = c(
    variety = "B73", star = "B73", gtf = "Zm-B73-REFERENCE-NAM-5.0.gtf",
    bed = "Zm-B73-REFERENCE-NAM-5.0.annotation_forStrandDetect_geneID.bed"
  ),
  maize_nc350_nam1 = c(
    variety = "NC350", star = "NC350", gtf = "Zm-NC350-REFERENCE-NAM-1.0.gtf",
    bed = "Zm-NC350-REFERENCE-NAM-1.0.annotation_forStrandDetect_geneID.bed"
  ),
  maize_w22_nrgene2 = c(
    variety = "W22", star = "W22", gtf = "Zm-W22-REFERENCE-NRGENE-2.0.clean.gtf",
    bed = "Zm-W22-REFERENCE-NRGENE-2.0.annotation_forStrandDetect_geneID.bed"
  )
)
for (reference_key in names(maize_reference_files)) {
  expected <- maize_reference_files[[reference_key]]
  selected_project <- maize_rna_project
  selected_project$genome_version <- reference_key
  resources <- app_env$genome_resources(selected_project)
  assert(
    identical(resources$variety, unname(expected[["variety"]])) &&
      identical(basename(resources$star_index), unname(expected[["star"]])) &&
      identical(basename(resources$gtf), unname(expected[["gtf"]])) &&
      identical(basename(resources$strand_bed), unname(expected[["bed"]])),
    paste(reference_key, "uses only its own variety-matched STAR index, GTF, and strand BED")
  )
}
assert(
  !app_env$rna_optional_quantifiers_available(maize_rna_project) &&
    !any(c("RSEM (optional)", "Kallisto (optional)") %in% app_env$sample_level_steps_for_project(maize_rna_project)),
  "maize RNA-seq excludes RSEM and Kallisto from available pipeline steps"
)
maize_atac_project <- within(maize_rna_project, { analysis_key <- "atac"; analysis <- "ATAC-seq" })
assert(
  identical(app_env$genome_species(maize_atac_project), "mouse") &&
    length(app_env$genome_reference_choices("maize", "ATAC-seq")) == 1L,
  "maize is restricted to RNA-seq and cannot become an ATAC reference"
)
rna_adapter_project <- within(chip_project, { analysis_key <- "rna"; analysis <- "RNA-seq" })
atac_adapter_project <- within(chip_project, { analysis_key <- "atac"; analysis <- "ATAC-seq" })
assert(identical(unname(app_env$default_adapter_pair(rna_adapter_project)), c("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")), "RNA-seq defaults to paired Illumina TruSeq adapters")
assert(identical(unname(app_env$default_adapter_pair(atac_adapter_project)), c("CTGTCTCTTATACACATCTCCGAGCCCACGAGAC", "CTGTCTCTTATACACATCTGACGCTGCCGACGA")), "ATAC-seq defaults to paired Nextera transposase adapters")
assert(app_env$adapter_input_prefix(rna_adapter_project) != app_env$adapter_input_prefix(atac_adapter_project), "adapter input state is isolated by analysis")
assert(identical(app_env$numeric_sort_kind(c("900 KB", "1.2 GB", "14 B")), "bytes"), "human-readable file sizes receive numeric table sorting")
assert(identical(app_env$numeric_sort_kind(c("9", "100", "2.5")), "numeric"), "numeric text receives numeric table sorting")
assert(identical(app_env$numeric_sort_kind(c("00:09:00", "01:00:00")), "duration"), "elapsed times receive duration sorting")
size_defs <- app_env$smart_table_column_defs(data.frame(size = c("900 KB", "1.2 GB"), stringsAsFactors = FALSE))
assert(any(vapply(size_defs, function(def) identical(if (is.null(def$type)) "" else def$type, "num") && !is.null(def$render) && grepl("Math.pow(1024", as.character(def$render), fixed = TRUE), logical(1))), "table renderer uses byte-aware numeric sort values")
rough_table <- structure(list(c(1, 2), I(list(c("a", "b"), "c"))), names = c("", ""), class = "data.frame", row.names = c(NA_integer_, -2L))
normalized_table <- app_env$normalize_csl_table_data(rough_table)
assert(!anyDuplicated(names(normalized_table)) && all(nzchar(names(normalized_table))) && !any(vapply(normalized_table, is.list, logical(1))), "table renderer sanitizes blank/duplicate names and list columns before DT serialization")
track_status <- data.frame(
  sample = c("complete_sample", "repair_sample", "unaligned_sample"),
  status = c("Complete", "Repair available", "Full Bowtie2 required"),
  stringsAsFactors = FALSE
)
assert(
  identical(app_env$cutrun_track_regeneration_candidates(track_status), c("complete_sample", "repair_sample")),
  "CUT&RUN signal-track regeneration accepts complete and repairable aligned samples without offering unaligned samples"
)
for (project_variant in list(
  RNA = within(chip_project, { analysis_key <- "rna"; analysis <- "RNA-seq" }),
  CUTRUN = within(chip_project, { analysis_key <- "cutrun"; analysis <- "CUT&RUN" }),
  ATAC = within(chip_project, { analysis_key <- "atac"; analysis <- "ATAC-seq" }),
  ChIP = chip_project
)) {
  step_meta <- app_env$run_step_meta(project_variant)
  assert(NROW(step_meta) == length(app_env$pipeline_order(project_variant)), paste(project_variant$analysis, "stepper descriptions match its pipeline"))
}
for (project_variant in list(
  within(chip_project, { analysis_key <- "cutrun"; analysis <- "CUT&RUN" }),
  chip_project
)) {
  assert(identical(tail(app_env$pipeline_order(project_variant), 1), "Peak Annotation"), paste(project_variant$analysis, "ends with Peak Annotation"))
}
cutrun_project <- within(chip_project, { analysis_key <- "cutrun"; analysis <- "CUT&RUN" })
assert(
  all(c("Peak-Calling Summary", "Differential Peak Summary") %in% app_env$cutrun_pipeline_order()),
  "CUT&RUN project summaries are first-class pipeline steps"
)
dir.create(file.path(root, "cutrun_summaries"), recursive = TRUE)
file.create(
  file.path(root, "cutrun_summaries", "peak_calling_summary_COMPLETE"),
  file.path(root, "cutrun_summaries", "differential_summary_COMPLETE")
)
cutrun_summary_status <- app_env$project_status(cutrun_project, jobs = data.frame(), progress = data.frame(), active_states = list())
assert(
  all(cutrun_summary_status$status[cutrun_summary_status$step %in% c("Peak-Calling Summary", "Differential Peak Summary")] == "Complete"),
  "CUT&RUN summary completion markers drive pipeline status"
)
summary_job_html <- as.character(app_env$cutrun_summary_job_ui(
  cutrun_project,
  data.frame(step = "Peak-Calling Summary", slurm_state = "RUNNING", job_id = "12345", elapsed = "00:01:02", stringsAsFactors = FALSE),
  "Peak-Calling Summary",
  "not_complete"
))
assert(
  grepl("summary-job-bar active", summary_job_html, fixed = TRUE) &&
    grepl("Job 12345", summary_job_html, fixed = TRUE),
  "CUT&RUN summary jobs visibly show an active progress bar and SLURM job ID"
)
atac_project_variant <- within(chip_project, { analysis_key <- "atac"; analysis <- "ATAC-seq" })
assert(identical(tail(app_env$pipeline_order(atac_project_variant), 1), "Differential Peaks"), "ATAC annotation runs inside MACS2 and DiffBind rather than as a final step")
assert(identical(app_env$canonical_job_step("peak_annotation"), "Peak Annotation"), "peak annotation job labels canonicalize")
assert(identical(app_env$step_data_paths(chip_project, "Peak Annotation"), file.path(root, "peak_annotation")), "peak annotation cleanup is confined to its project folder")

blank_editor <- app_env$blank_design_matrix_rows(c("condition", "replicate"), rows = 3)
assert(NROW(blank_editor) == 3L && all(!blank_editor$include), "blank design setup provides editable excluded rows")
assert(
  grepl('actionButton("add_design_rows", "Add 1 blank row"', app_text, fixed = TRUE) &&
    grepl('actionButton("add_design_rows", "Add 1 sample row"', app_text, fixed = TRUE) &&
    !grepl("Add 5 blank rows", app_text, fixed = TRUE) &&
    !grepl("Add 5 sample rows", app_text, fixed = TRUE) &&
    grepl("blank_design_matrix_rows(metadata, rows = 1)", server_source, fixed = TRUE),
  "each design-row action adds exactly one editable row"
)
blank_form <- as.character(app_env$design_form_table_ui(blank_editor))
assert(grepl("design_form_1_sample", blank_form, fixed = TRUE) && grepl("design_form_1_filename", blank_form, fixed = TRUE), "blank design setup renders visible text inputs")
provided_editor <- app_env$design_editor_from_project(chip_project)
form_values <- list()
form_values[[app_env$design_form_input_id(1, "treatment")]] <- "edited_treatment"
form_values[[app_env$design_form_input_id(1, "include")]] <- FALSE
edited_design <- app_env$apply_design_form_values(provided_editor, form_values)
assert(identical(edited_design$treatment[[1]], "edited_treatment") && !edited_design$include[[1]], "provided design matrices remain editable through visible form controls")
scrna_editor <- data.frame(
  sample_id = c("sample_1", "sample_2"), input_path = c("/input/one", "/input/two"),
  input_type = "filtered_10x_matrix", condition = c("control", "treated"), stringsAsFactors = FALSE
)
scrna_form <- as.character(app_env$design_form_table_ui(scrna_editor, prefix = "new_scrna_form", columns = names(scrna_editor)))
scrna_values <- list(new_scrna_form_2_condition = "edited_treated")
edited_scrna <- app_env$apply_design_form_values(scrna_editor, scrna_values, prefix = "new_scrna_form")
assert(
  grepl("new_scrna_form_1_sample_id", scrna_form, fixed = TRUE) &&
    grepl("new_scrna_form_2_input_path", scrna_form, fixed = TRUE) &&
    identical(edited_scrna$condition[[2]], "edited_treated"),
  "scRNA manifests use the same visible input-box editor and preserve form edits"
)

duplicate_design <- data.frame(
  include = TRUE, sample = c("sample-A", "sample A"), cell_type = "", condition = c("A", "B"),
  replicate = 1:2, filename = c("a.fastq.gz", "b.fastq.gz"), status = "", stringsAsFactors = FALSE
)
atac_design_project <- chip_project
atac_design_project$analysis_key <- "atac"
atac_design_project$analysis <- "ATAC-seq"
assert(
  identical(app_env$default_metadata_cols(atac_design_project), c("condition", "replicate")),
  "new ATAC designs default to only condition and replicate"
)
duplicate_error <- tryCatch({
  app_env$write_design_matrix(atac_design_project, duplicate_design, c("condition", "replicate"))
  ""
}, error = conditionMessage)
assert(grepl("remain unique", duplicate_error), "filesystem-safe sample collisions rejected")

blank_design <- duplicate_design[1, , drop = FALSE]
blank_design$sample <- "sample1"
blank_design$filename <- ""
blank_error <- tryCatch({ app_env$write_design_matrix(atac_design_project, blank_design, c("condition", "replicate")); "" }, error = conditionMessage)
assert(grepl("FASTQ filename", blank_error), "blank included FASTQ filenames rejected")

valid_design <- duplicate_design[1, , drop = FALSE]
valid_design$sample <- "sample1"
saved_design <- app_env$write_design_matrix(atac_design_project, valid_design, c("condition", "replicate"))
assert(file.exists(saved_design) && file.info(saved_design)$size > 0, "design matrix saved atomically")
saved_table <- app_env$safe_read_table(saved_design)
assert(
  all(c("condition", "replicate") %in% names(saved_table)) && !"cell_type" %in% names(saved_table),
  "new ATAC designs omit the optional cell_type column"
)

legacy_atac_project <- atac_design_project
legacy_atac_project$design_matrix_path <- saved_design
legacy_atac_table <- saved_table
legacy_atac_table$cell_type <- "CD4"
legacy_atac_table <- legacy_atac_table[, c("sample", "cell_type", "condition", "replicate", "filename"), drop = FALSE]
write.table(legacy_atac_table, saved_design, sep = "\t", row.names = FALSE, quote = FALSE)
assert(
  "cell_type" %in% app_env$project_metadata_cols(legacy_atac_project),
  "existing ATAC projects retain an explicitly saved cell_type column"
)

unsafe_design <- valid_design
unsafe_design$condition <- "A\tB"
unsafe_error <- tryCatch({ app_env$write_design_matrix(atac_design_project, unsafe_design, c("condition", "replicate")); "" }, error = conditionMessage)
assert(grepl("tabs or line breaks", unsafe_error), "tab characters rejected before TSV save")

for (key in c("rna", "cutrun", "atac", "chip")) {
  example <- app_env$example_dataset_paths(key)
  example_design <- file.path(example$design_dir, "design_matrix.txt")
  assert(dir.exists(example$fastq_dir) && file.exists(example_design), paste(key, "bundled example paths exist"))
  table <- app_env$safe_read_table(example_design)
  assert(NROW(table) > 0 && !anyDuplicated(table$sample), paste(key, "bundled example design has unique samples"))
  reads <- trimws(unlist(strsplit(as.character(table$filename), "[;,]")))
  read_paths <- file.path(example$fastq_dir, reads[nzchar(reads)])
  assert(all(file.exists(read_paths)), paste(key, "bundled example FASTQs match the design"))
  readable <- vapply(read_paths, function(path) {
    connection <- gzfile(path, open = "rt")
    on.exit(close(connection), add = TRUE)
    length(readLines(connection, n = 4L, warn = FALSE)) == 4L
  }, logical(1))
  assert(all(readable), paste(key, "bundled example FASTQs are readable gzip data"))
}
assert(identical(app_env$example_dataset_paths("rna")$name, "example_dataset"), "RNA example project name matches the CodeSpringLab notebook output folder")
rna_example <- app_env$example_dataset_paths("rna")
rna_example_manifest <- file.path(rna_example$design_dir, "design_matrix.txt")
rna_example_design <- app_env$safe_read_table(rna_example_manifest)
rna_example_results <- file.path(root, "example_dataset", "data")
rna_example_cutadapt <- file.path(rna_example_results, "cutadapt")
dir.create(rna_example_cutadapt, recursive = TRUE, showWarnings = FALSE)
rna_example_reads <- trimws(unlist(strsplit(as.character(rna_example_design$filename), "[;,]")))
for (read in rna_example_reads[nzchar(rna_example_reads)]) {
  writeBin(as.raw(rep(seq_len(100), 2)), file.path(rna_example_cutadapt, basename(read)))
}
rna_example_project <- list(
  id = "rna/example_dataset", name = "example_dataset", analysis_key = "rna", analysis = "RNA-seq",
  design_matrix_path = rna_example_manifest, data_dir = rna_example_results, results_root = dirname(dirname(rna_example_results)),
  fastq_dir = rna_example$fastq_dir, fastq_dirs = rna_example$fastq_dir, paired_end = TRUE, genome = "mouse"
)
rna_example_progress <- app_env$sample_progress(rna_example_project, jobs = data.frame())$table
rna_example_cutadapt_progress <- rna_example_progress[rna_example_progress$step == "Cutadapt", , drop = FALSE]
assert(NROW(rna_example_cutadapt_progress) == NROW(rna_example_design) && all(rna_example_cutadapt_progress$status == "Completed"), "RNA example manifest maps every Cutadapt R1/R2 output")
rna_example_status <- app_env$project_status(rna_example_project, jobs = data.frame(), progress = rna_example_progress)
assert(identical(rna_example_status$status[rna_example_status$step == "Cutadapt"], "Complete"), "RNA example Cutadapt step reports Complete when all manifest outputs exist")
cutrun_example <- app_env$safe_read_table(file.path(app_env$example_dataset_paths("cutrun")$design_dir, "design_matrix.txt"))
assert(sum(cutrun_example$target_class == "control") == 2L, "CUT&RUN example has explicit matched controls")
assert(all(c("cell_type", "mark", "target_class", "condition", "replicate", "control_sample") %in% names(cutrun_example)), "CUT&RUN example contains the editable assay metadata")
cutrun_targets <- cutrun_example$target_class != "control"
assert(all(nzchar(cutrun_example$control_sample[cutrun_targets])), "CUT&RUN example assigns every target to an IgG control")

individual_peak_file <- file.path(root, "individual_sample.stringent.bed")
writeLines(c(
  "chr1\t100\t200\t10",
  "chr1\t300\t450\t30",
  "chr2\t500\t600\t20"
), individual_peak_file)
individual_navigation <- app_env$cutrun_individual_peak_navigation(individual_peak_file, max_peaks = 2L)
assert(
  identical(unname(individual_navigation$peaks[[1]]), "chr1:301-450") &&
    identical(individual_navigation$total, 3L) && identical(individual_navigation$shown, 2L),
  "individual CUT&RUN peak navigation exposes bounded peak choices ordered by signal while retaining total peak count"
)
browser_peak_root <- file.path(root, "seacr", "spikein_non_stringent", "target")
dir.create(browser_peak_root, recursive = TRUE)
browser_peak_catalog <- data.frame(
  sample = c("target", "target"), kind = c("peaks", "peaks"), format = c("bed", "narrowPeak"),
  label = c("SEACR", "MACS2"),
  path = c(file.path(browser_peak_root, "target.stringent.bed"), file.path(root, "macs2", "target", "target_peaks.narrowPeak")),
  stringsAsFactors = FALSE
)
browser_peak_rows <- app_env$cutrun_browser_peak_rows(
  within(chip_project, { analysis_key <- "cutrun"; analysis <- "CUT&RUN" }), browser_peak_catalog, "target"
)
assert(
  identical(sort(browser_peak_rows$tool), c("MACS2", "SEACR")) &&
    any(grepl("Spike-in-scaled", browser_peak_rows$parameters, fixed = TRUE)),
  "CUT&RUN individual browser separates peak-calling tool from the selected parameter set"
)
browser_signal_catalog <- data.frame(
  sample = c("target", "target", "igg", "igg"), kind = "signal", format = "bigwig", label = "signal",
  path = c("/tmp/target_spikein.bw", "/tmp/target_CPM.bw", "/tmp/igg_CPM.bw", "/tmp/igg_raw.bw"),
  stringsAsFactors = FALSE
)
assert(
  identical(app_env$cutrun_browser_signal_modes(browser_signal_catalog, c("target", "igg")), "cpm"),
  "individual browser offers matched target/IgG normalization modes rather than unrelated tracks"
)
assert(
  identical(app_env$cutrun_browser_signal_modes_for_peak_call(c("spikein", "cpm", "raw"), "SEACR", "Raw track / SEACR norm / relaxed"), c("spikein", "cpm", "raw")) &&
    identical(app_env$cutrun_browser_signal_modes_for_peak_call(c("spikein", "cpm", "raw"), "MACS2", "narrow peaks; q ≤ 0.01"), c("spikein", "cpm", "raw")),
  "individual browser keeps CPM available for display regardless of the selected peak caller input"
)
assert(
  grepl("cutrun_peak_mode", server_source, fixed = TRUE) &&
    grepl("cutrun_control_sample_for(p, target_sample)", server_source, fixed = TRUE),
  "CUT&RUN individual browser mode pairs each target with its matched IgG"
)
assert(
  grepl("observeEvent(input$genome_browser_cutrun_sample", server_source, fixed = TRUE) &&
    grepl("send_genome_browser()", server_source, fixed = TRUE),
  "changing a CUT&RUN browser target immediately reloads its target and matched-IgG tracks"
)
assert(
  grepl("Signal normalization for display", server_source, fixed = TRUE) &&
    grepl("updateSelectInput(session, \"genome_browser_cutrun_signal_normalization\",", server_source, fixed = TRUE) &&
    grepl('selected = "cpm"', server_source, fixed = TRUE),
  "CUT&RUN browser resets CPM as the visualization default when peak caller/settings change"
)
assert(
  exists("write_cutrun_peak_summary_xlsx", envir = app_env, inherits = FALSE) &&
    grepl("download_cutrun_seacr_peak_summary_xlsx", server_source, fixed = TRUE) &&
    grepl("Shared Peak Overlap Details", paste(deparse(app_env$write_cutrun_peak_summary_xlsx), collapse = "\n"), fixed = TRUE),
  "CUT&RUN peak QC includes an Excel summary with overlap details"
)

atac_project <- chip_project
atac_project$analysis_key <- "atac"
atac_project$analysis <- "ATAC-seq"
initial_progress <- app_env$sample_progress(atac_project, jobs = data.frame())$table
initial_a1_bowtie <- initial_progress$status[initial_progress$sample == "A1" & initial_progress$step == "Bowtie2"]
assert(identical(initial_a1_bowtie, "Not started"), "untouched samples start as Not started")
cutadapt_dir <- file.path(root, "cutadapt")
dir.create(cutadapt_dir, recursive = TRUE, showWarnings = FALSE)
partial_cutadapt_output <- file.path(cutadapt_dir, "A1.fastq.gz")
writeBin(as.raw(rep(seq_len(100), 2)), partial_cutadapt_output)
assert(
  isTRUE(app_env$trimmed_checkbox_default(atac_project, FALSE)),
  "ATAC trimmed-read controls default on when Cutadapt output is available"
)
partial_cutadapt_progress <- app_env$sample_progress(atac_project, jobs = data.frame())$table
assert(
  identical(partial_cutadapt_progress$status[partial_cutadapt_progress$sample == "A1" & partial_cutadapt_progress$step == "Cutadapt"], "Completed"),
  "existing validated Cutadapt output is complete at sample level"
)
partial_cutadapt_status <- app_env$project_status(atac_project, jobs = data.frame(), progress = partial_cutadapt_progress)
assert(
  identical(partial_cutadapt_status$status[partial_cutadapt_status$step == "Cutadapt"], "Partial"),
  "mixed completed and untouched Cutadapt samples report Partial rather than Not started"
)
assert(identical(partial_cutadapt_status$detail[partial_cutadapt_status$step == "Cutadapt"], "1/8 samples complete"), "partial Cutadapt status reports the completed sample count")
assert(
  lengths(regmatches(app_text, gregexpr('updateCheckboxInput(session, "atac_bowtie2_use_trimmed", value = TRUE)', app_text, fixed = TRUE))) >= 2L,
  "project switching and Cutadapt submission activate the ATAC trimmed-read default"
)
assert(identical(app_env$status_css_key("Partial"), "partial"), "partial pipeline status has a dedicated visual state")
unlink(partial_cutadapt_output)
partial_targets <- app_env$sample_step_targets(atac_project, "A1", "Bowtie2")
dir.create(dirname(partial_targets[[1]]), recursive = TRUE, showWarnings = FALSE)
writeLines("partial", partial_targets[[1]])
partial_progress <- app_env$sample_progress(atac_project, jobs = data.frame())$table
partial_a1_bowtie <- partial_progress$status[partial_progress$sample == "A1" & partial_progress$step == "Bowtie2"]
assert(identical(partial_a1_bowtie, "Not started"), "partial files do not imply failure before a terminal job state")
unlink(partial_targets[[1]])
running_job <- data.frame(step = "Bowtie2", sample = "A1", slurm_state = "RUNNING", elapsed = "00:00:05", stderr = "", stringsAsFactors = FALSE)
running_progress <- app_env$sample_progress(atac_project, jobs = running_job)$table
assert(identical(running_progress$status[running_progress$sample == "A1" & running_progress$step == "Bowtie2"], "Running"), "active jobs are not marked failed while outputs are incomplete")
finished_job <- running_job
finished_job$slurm_state <- "COMPLETED"
finished_progress <- app_env$sample_progress(atac_project, jobs = finished_job)$table
assert(identical(finished_progress$status[finished_progress$sample == "A1" & finished_progress$step == "Bowtie2"], "Likely failed"), "missing outputs are classified only after a job completes")
retry_ui <- as.character(app_env$sample_retry_ui(atac_project, finished_progress, "Bowtie2"))
assert(grepl("atac_bowtie2_samples", retry_ui, fixed = TRUE) && grepl('wanted=[&quot;A1&quot;]', retry_ui, fixed = TRUE), "retry action selects only terminally incomplete samples before submission")
sample_dir <- file.path(root, "macs2", "A1")
dir.create(sample_dir, recursive = TRUE)
legacy_peak <- file.path(sample_dir, "A1_peaks.narrowPeak")
run_log <- file.path(sample_dir, "A1_macs2.log")
marker <- file.path(sample_dir, "A1_macs2_complete.txt")
writeLines("chr1\t1\t2", legacy_peak)
assert(identical(app_env$atac_macs2_completion_target(atac_project, "A1"), legacy_peak), "legacy ATAC peaks remain recognized")
writeLines(c(
  "chr1\t1\t200\tpeak1\t50\t.\t8\t6\t5\t75",
  "chr1\t300\t500\tpeak2\t40\t.\t7\t5\t4\t60"
), legacy_peak)
atac_file_catalog <- app_env$atac_result_file_catalog(atac_project)
catalog_peak <- atac_file_catalog[atac_file_catalog[["Absolute path"]] == normalizePath(legacy_peak), , drop = FALSE]
assert(
  NROW(catalog_peak) == 1L && identical(catalog_peak$Tool[[1]], "MACS2") &&
    identical(catalog_peak$Sample[[1]], "A1") && identical(catalog_peak$File[[1]], basename(legacy_peak)),
  "ATAC file catalog reports RNA-style tool, sample, filename, size, modified time, and protected absolute path metadata"
)
atac_peak_catalog <- app_env$genome_browser_track_catalog(atac_project)
assert(
  NROW(atac_peak_catalog[atac_peak_catalog$kind == "peaks" & atac_peak_catalog$sample == "A1", , drop = FALSE]) >= 1L,
  "ATAC genome browser discovers per-sample MACS2 peak tracks"
)
atac_peak_navigation <- app_env$cutrun_individual_peak_navigation(legacy_peak)
assert(
  length(atac_peak_navigation$peaks) == 2L &&
    identical(unname(atac_peak_navigation$peaks)[[1]], "chr1:2-200") &&
    grepl("MACS2 -log10(q)", names(atac_peak_navigation$peaks)[[1]], fixed = TRUE) &&
    !any(grepl("SEACR", names(atac_peak_navigation$peaks), fixed = TRUE)),
  "ATAC sample-peak navigation uses MACS2 evidence labels without CUT&RUN terminology"
)
atac_signal_navigation <- app_env$cutrun_individual_peak_navigation(
  legacy_peak,
  score_preference = c("signalValue", "qValue", "pValue", "score"),
  locus_flank = 500L
)
assert(
  grepl("signal 8", names(atac_signal_navigation$peaks)[[1]], fixed = TRUE) &&
    identical(unname(atac_signal_navigation$peaks)[[1]], "chr1:1-700") &&
    identical(atac_signal_navigation$ranking, "MACS2 signalValue"),
  "ATAC peak navigation ranks by MACS2 signalValue and adds flanking baseline around the called interval"
)
missing_q_peak <- file.path(sample_dir, "A1_missing_q_peaks.narrowPeak")
writeLines(c(
  "chr1\t10\t110\tpeak_missing_q\t25\t.\t12.5\t-1\t-1\t50",
  "chr1\t210\t310\tpeak_partial\t20\t.\t.\t-1\t-1\t45"
), missing_q_peak)
missing_q_navigation <- app_env$cutrun_individual_peak_navigation(missing_q_peak)
assert(
  grepl("signal 12.5", names(missing_q_navigation$peaks)[[1]], fixed = TRUE) &&
    !any(grepl(" NA", names(missing_q_navigation$peaks), fixed = TRUE)),
  "ATAC peak navigation falls back from unavailable MACS2 q/p values to finite signal and never prints signal NA"
)
atac_display_bed <- app_env$genome_browser_atac_peak_bed6(legacy_peak)
atac_display_rows <- utils::read.table(atac_display_bed, sep = "\t", header = FALSE, quote = "", comment.char = "")
assert(
  file.exists(atac_display_bed) && NCOL(atac_display_rows) == 6L && NROW(atac_display_rows) == 2L &&
    identical(as.character(atac_display_rows[[1]][[1]]), "chr1") && identical(as.numeric(atac_display_rows[[5]][[1]]), 50),
  "ATAC individual peaks are converted to a visible BED6 annotation while quantitative signal remains in the bigWig"
)
summits_file <- file.path(sample_dir, "A1_summits.bed")
writeLines("chr1\t75\t76\tpeak1\t8", summits_file)
atac_catalog_with_summits <- app_env$genome_browser_track_catalog(atac_project)
atac_called_peak_rows <- app_env$atac_individual_peak_rows(atac_catalog_with_summits, "A1")
assert(
  NROW(atac_called_peak_rows) >= 1L &&
    all(grepl("_peaks\\.(narrowPeak|broadPeak)$", basename(atac_called_peak_rows$path), ignore.case = TRUE)) &&
    !any(grepl("summits", basename(atac_called_peak_rows$path), ignore.case = TRUE)),
  "ATAC individual peak mode excludes coordinate-only MACS2 summit BEDs when scored called-peak files are available"
)
completed_selector_samples <- app_env$completed_samples_for_step(atac_project, "MACS2 Peaks", c("A1", "A2"))
assert("A1" %in% completed_selector_samples && !"A2" %in% completed_selector_samples, "sample selectors identify completed samples without unchecking unfinished samples")
failed_macs_job <- data.frame(step = "MACS2 Peaks", sample = "A1", slurm_state = "FAILED", elapsed = "00:01:00", stderr = "", stringsAsFactors = FALSE)
legacy_peak_progress <- app_env$sample_progress(atac_project, jobs = failed_macs_job)$table
legacy_peak_status <- legacy_peak_progress$status[legacy_peak_progress$sample == "A1" & legacy_peak_progress$step == "MACS2 Peaks"]
assert(identical(legacy_peak_status, "Completed"), "validated legacy ATAC peaks are not hidden by a later failed retry")
writeLines(rep("validated alignment summary", 8), partial_targets[[1]])
failed_bowtie_job <- data.frame(step = "Bowtie2", sample = "A1", slurm_state = "FAILED", elapsed = "00:01:00", stderr = "", stringsAsFactors = FALSE)
failed_bowtie_progress <- app_env$sample_progress(atac_project, jobs = failed_bowtie_job)$table
assert(identical(failed_bowtie_progress$status[failed_bowtie_progress$sample == "A1" & failed_bowtie_progress$step == "Bowtie2"], "Likely failed"), "legacy-output exception is limited to ATAC MACS2 peaks")
unlink(partial_targets[[1]])
writeLines("Traceback (most recent call last):\nOSError: No space left on device", run_log)
assert(identical(app_env$atac_macs2_completion_target(atac_project, "A1"), marker), "new ATAC runs require a completion marker")
assert(app_env$cutrun_macs_fatal_error_signal(atac_project, data.frame(), "MACS2 Peaks", "A1"), "ATAC internal MACS2 exception detection")
writeLines("status\tcomplete", marker)
assert(identical(app_env$atac_macs2_completion_target(atac_project, "A1"), marker), "completed ATAC marker selected")
assert(identical(app_env$atac_macs2_peak_file(atac_project, "A1"), legacy_peak), "validated ATAC peak selected for DiffBind")

unlink(marker)
assert(identical(app_env$chip_macs2_peak_file(chip_project, "A1"), ""), "partial ChIP MACS2 output rejected")
writeLines("status\tcomplete", marker)
assert(identical(app_env$chip_macs2_peak_file(chip_project, "A1"), legacy_peak), "completed ChIP MACS2 peak accepted")
chip_peaks <- app_env$chip_peak_summary_table(chip_project)
assert(NROW(chip_peaks) == 4L && chip_peaks$status[chip_peaks$sample == "A1"] == "Completed", "ChIP matched-input peak summary reports completion")
alignment_dir <- file.path(root, "bowtie2", "A1")
dir.create(alignment_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(c("sample\tA1", "mapped_reads\t100", "deduplicated_reads\t80", "bigwig_normalization\tCPM"), file.path(alignment_dir, "A1_alignment_summary.txt"))
signal_file <- file.path(alignment_dir, "A1Aligned.sortedByCoord_removeDup.out.bw")
writeBin(as.raw(seq_len(64)), signal_file)
fragment_pdf <- file.path(alignment_dir, "A1_insert_size_histogram.pdf")
grDevices::pdf(fragment_pdf, width = 8, height = 5)
graphics::plot(1:10, type = "h", main = "Synthetic insert sizes")
grDevices::dev.off()
fragment_html <- as.character(app_env$fragment_plot_ui(fragment_pdf))
assert(grepl("fragment-plot-frame", fragment_html, fixed = TRUE) && grepl("data:image/png;base64", fragment_html, fixed = TRUE) && !grepl("iframe", fragment_html, fixed = TRUE), "fragment PDFs render directly as standardized images")
assert(identical(app_env$pdf_first_page_data_uri(fragment_pdf), app_env$pdf_first_page_data_uri(fragment_pdf)), "fragment PDF rendering is cached in memory")
chip_alignment <- app_env$chip_alignment_summary_table(chip_project)
assert(NROW(chip_alignment) == 1L && all(c("role", "condition", "matched_input") %in% names(chip_alignment)) && chip_alignment$matched_input[[1]] == "I1", "ChIP alignment summary includes experimental roles")
chip_signal <- app_env$peak_signal_track_table(chip_project)
assert(NROW(chip_signal) == 1L && chip_signal$role[[1]] == "chip" && chip_signal$normalization[[1]] == "CPM", "ChIP signal table reports role and saved normalization")
igv_catalog <- app_env$genome_browser_track_catalog(chip_project)
assert(NROW(igv_catalog) >= 2L && all(c("signal", "peaks") %in% igv_catalog$kind), "embedded genome browser catalogs signal and peak tracks")
assert(identical(app_env$genome_browser_reference(chip_project), "mm39") && identical(app_env$genome_browser_reference(human_chip_project), "hg38"), "embedded genome browser follows the project reference")
shared_scale_config <- app_env$genome_browser_signal_display_config(TRUE, TRUE)
independent_scale_config <- app_env$genome_browser_signal_display_config(FALSE, TRUE)
assert(identical(shared_scale_config$autoscaleGroup, "codespring_comparison_signal"), "comparison bigWigs share an IGV autoscale group")
assert(is.null(independent_scale_config$autoscaleGroup), "manual genome-browser tracks retain independent autoscaling")
range_response <- app_env$genome_browser_range_response(
  list(path = signal_file, content_type = "application/octet-stream"),
  list(REQUEST_METHOD = "GET", HTTP_RANGE = "bytes=10-19")
)
assert(identical(range_response$status, 206) && length(range_response$content) == 10L, "genome browser serves bounded byte ranges for large tracks")
invalid_range <- app_env$genome_browser_range_response(
  list(path = signal_file, content_type = "application/octet-stream"),
  list(REQUEST_METHOD = "GET", HTTP_RANGE = "bytes=999-1000")
)
assert(identical(invalid_range$status, 416), "genome browser rejects out-of-range project track requests")
fake_file_choices <- stats::setNames(c(signal_file, file.path(root, "bowtie2", "B1", "B1_signal.bw")), c("A1 signal", "B1 signal"))
assert(identical(app_env$result_file_sample(chip_project, signal_file), "A1"), "result files resolve to their design sample")
assert(identical(unname(app_env$filter_result_files_by_sample(chip_project, fake_file_choices, "A1")), signal_file), "sample file filter excludes other samples")
assert(identical(app_env$validated_project_result_path(chip_project, signal_file), normalizePath(signal_file)), "current-project result path accepted")
outside_file <- tempfile("outside-result-")
writeLines("outside", outside_file)
assert(identical(app_env$validated_project_result_path(chip_project, outside_file), ""), "result path outside the current project is rejected")
assert(inherits(app_env$atac_summary_cards_ui(chip_project), "shiny.tag"), "ChIP summary cards render with fake results")
assert(inherits(app_env$chip_results_explorer_ui(), "shiny.tag"), "ChIP Results Explorer UI renders locally")

atac_ui_text <- as.character(app_env$atac_results_explorer_ui())
chip_ui_text <- as.character(app_env$chip_results_explorer_ui())
cutrun_ui_text <- as.character(app_env$cutrun_results_explorer_ui())
atac_summary_text <- as.character(app_env$atac_summary_cards_ui(atac_project))
assert(
  grepl("Disk space", atac_summary_text, fixed = TRUE) &&
    grepl("Peak files", atac_summary_text, fixed = TRUE) &&
    grepl("Result files", atac_summary_text, fixed = TRUE) &&
    grepl("atac_overview_sample_progress_ui", atac_ui_text, fixed = TRUE) &&
    !grepl("Completion state for every ATAC-seq analysis stage", atac_ui_text, fixed = TRUE) &&
    !grepl("save_results_design", atac_ui_text, fixed = TRUE),
  "ATAC overview mirrors the compact RNA summary and live sample-progress layout"
)
assert(
  grepl("atac_file_tool_ui", atac_ui_text, fixed = TRUE) &&
    grepl("atac_file_sample_ui", atac_ui_text, fixed = TRUE) &&
    grepl("atac_files_catalog_table", atac_ui_text, fixed = TRUE) &&
    grepl("download_atac_file_catalog", atac_ui_text, fixed = TRUE) &&
    grepl("download_selected_atac_file", atac_ui_text, fixed = TRUE) &&
    grepl("confirm_delete_atac_file", server_source, fixed = TRUE) &&
    grepl("validated_project_result_path", server_source, fixed = TRUE),
  "ATAC Files matches the RNA file catalog with tool/sample filters and protected copy, download, and confirmed-delete actions"
)
for (ui_check in list(
  ATAC = atac_ui_text,
  ChIP = chip_ui_text,
  CUTRUN = cutrun_ui_text
)) {
  assert(grepl("Developed by CSHL's Bioinformatics Shared Resource", ui_check, fixed = TRUE), paste(names(ui_check), "Results Explorer uses the shared branded header"))
  assert(grepl("Overview", ui_check, fixed = TRUE) && grepl("QC", ui_check, fixed = TRUE) && grepl("Files", ui_check, fixed = TRUE), "custom Results Explorer exposes the standard navigation")
  assert(grepl("Signal &amp; Peaks", ui_check, fixed = TRUE) || grepl("Signal & Peaks", ui_check, fixed = TRUE), "custom Results Explorer exposes standardized signal and peak navigation")
}
assert(all(vapply(list(atac_ui_text, chip_ui_text), function(x) grepl('col-sm-3', x, fixed = TRUE) && grepl('col-sm-9', x, fixed = TRUE), logical(1))), "ATAC and ChIP use the same 3:9 control/content layout as CUT&RUN")
assert(grepl("Initial QC", chip_ui_text, fixed = TRUE) && grepl("Fragment Size", chip_ui_text, fixed = TRUE), "ChIP Results Explorer includes RNA-style QC navigation")
assert(grepl("Signal Tracks", atac_ui_text, fixed = TRUE) && grepl("Signal Tracks", chip_ui_text, fixed = TRUE), "ATAC and ChIP Results Explorers expose signal-track navigation")
assert(all(vapply(list(atac_ui_text, chip_ui_text, cutrun_ui_text), grepl, logical(1), pattern = "Genome Browser", fixed = TRUE)), "ATAC, ChIP, and CUT&RUN Results Explorers expose the embedded genome browser")
assert(
  grepl("Available controls are detected from the selected project's saved signal and peak files", atac_ui_text, fixed = TRUE) &&
    !grepl("For CUT&RUN, Individual sample peak mode", atac_ui_text, fixed = TRUE) &&
    !grepl('output$genome_browser_controls_ui <- renderUI({\n    req(identical(input$web_main_tabs', server_source, fixed = TRUE),
  "shared genome browser uses assay-neutral guidance and renders controls in completed-results/read-only layouts"
)
assert(
  grepl("genome_browser_catalog <- reactive", server_source, fixed = TRUE) &&
    grepl('uiOutput("genome_browser_sample_peak_navigation_ui")', server_source, fixed = TRUE) &&
    grepl("output$genome_browser_sample_peak_navigation_ui <- renderUI", server_source, fixed = TRUE) &&
    grepl("is.null(input$genome_browser_sample_peak_sample)", server_source, fixed = TRUE),
  "ATAC genome-browser track discovery is cached and lightweight sample/file controls render before peak ranking and browser loading"
)
assert(
  grepl("peak_loci <- navigation$peaks", server_source, fixed = TRUE) &&
    grepl("selected_locus <- unname(navigation$peaks)[[1]]", server_source, fixed = TRUE),
  "ATAC peak navigation has no blank placeholder and selects the first ranked peak by default"
)
atac_signal_position <- regexpr("Signal &amp; Peaks", atac_ui_text, fixed = TRUE)[[1]]
atac_browser_position <- regexpr("Genome Browser", atac_ui_text, fixed = TRUE)[[1]]
atac_diff_position <- regexpr("Differential Accessibility", atac_ui_text, fixed = TRUE)[[1]]
assert(atac_signal_position > 0L && atac_browser_position > atac_signal_position && atac_diff_position > atac_browser_position, "ATAC Genome Browser is a main tab between Signal & Peaks and Differential Accessibility")
cutrun_signal_position <- regexpr("Signal &amp; Peaks", cutrun_ui_text, fixed = TRUE)[[1]]
cutrun_browser_position <- regexpr("Genome Browser", cutrun_ui_text, fixed = TRUE)[[1]]
cutrun_diff_position <- regexpr("Differential Binding", cutrun_ui_text, fixed = TRUE)[[1]]
assert(cutrun_signal_position > 0L && cutrun_browser_position > cutrun_signal_position && cutrun_diff_position > cutrun_browser_position, "CUT&RUN Genome Browser is a main tab between Signal & Peaks and Differential Binding")
assert(grepl("Open comparison in Genome Browser", cutrun_ui_text, fixed = TRUE), "CUT&RUN Differential Binding links directly to its comparison browser")
assert(grepl("codespring-genome-browser-controls", atac_ui_text, fixed = TRUE), "genome-browser controls use an overflow-safe dropdown container")
assert(all(vapply(list(atac_ui_text, chip_ui_text, cutrun_ui_text), grepl, logical(1), pattern = "Gene Annotation", fixed = TRUE)), "all peak Results Explorers expose gene annotations")
assert(grepl("cutrun_file_sample_ui", cutrun_ui_text, fixed = TRUE), "CUT&RUN file explorer exposes a sample selector")
assert(grepl("height:680px", app_env$app_css, fixed = TRUE), "fragment plots share a fixed display height")

rna_project <- chip_project
rna_project$id <- "fake-rna"
rna_project$name <- "fake-rna"
rna_project$analysis_key <- "rna"
rna_project$analysis <- "RNA-seq"
rna_project$genome <- "mouse"
rna_project$genome_version <- "mouse_gencodeM39"
old_app_home <- app_env$APP_HOME
app_env$APP_HOME <- file.path(root, "fake-app-home")
rna_config <- app_env$write_native_shiny_config(rna_project)
rna_config_text <- paste(readLines(rna_config, warn = FALSE), collapse = "\n")
assert(grepl('genome_species <- "mouse"', rna_config_text, fixed = TRUE), "RNA Results Explorer config records the analysis species")
assert(grepl('genome_version <- "mouse_gencodeM39"', rna_config_text, fixed = TRUE), "RNA Results Explorer config records the analysis reference")
assert(grepl("gencode.vM39.primary_assembly.annotation.gtf", rna_config_text, fixed = TRUE), "RNA Results Explorer config receives the analysis GTF")
rna_viewer <- app_env$load_native_rnaseq_viewer(rna_project)
app_env$APP_HOME <- old_app_home
assert(inherits(rna_viewer$ui, "shiny.tag") && is.function(rna_viewer$server), "RNA Results Explorer loads against a synthetic project")
rna_ui_text <- as.character(rna_viewer$ui)
assert(grepl("RNA-seq Results Explorer", rna_ui_text, fixed = TRUE), "RNA Results Explorer uses the canonical analysis label")
assert(grepl("Overview", rna_ui_text, fixed = TRUE) && grepl("QC", rna_ui_text, fixed = TRUE) && grepl("Files", rna_ui_text, fixed = TRUE), "RNA Results Explorer now follows the shared navigation contract")
assert(
  grepl("Sample progress", rna_ui_text, fixed = TRUE) &&
    grepl("rna_overview_sample_progress_ui", rna_ui_text, fixed = TRUE) &&
    !grepl("Pipeline status", rna_ui_text, fixed = TRUE) &&
    !grepl("Design matrix", rna_ui_text, fixed = TRUE),
  "RNA overview exposes the live sample progress matrix without redundant pipeline or design summaries"
)
assert(
  grepl("rna_file_tool", rna_ui_text, fixed = TRUE) &&
    grepl("rna_file_sample_ui", rna_ui_text, fixed = TRUE) &&
    grepl("download_selected_rna_file", rna_ui_text, fixed = TRUE) &&
    grepl("rna-trigger-file-download", rna_ui_text, fixed = TRUE) &&
    !grepl("rna_file_category", rna_ui_text, fixed = TRUE),
  "RNA Files tab filters the readable filename catalog and exposes secure file actions"
)

fake_jobs <- data.frame(
  step = c("Bowtie2", "Bowtie2", "Bowtie2", "Bowtie2", "FastQC"),
  sample = c("A1", "A2", "A3", "A4", "A1"),
  job_id = c("101", "102", "103", "104", "105"),
  slurm_state = c("RUNNING", "PENDING", "COMPLETED", "CANCELLED", "RUNNING"),
  stringsAsFactors = FALSE
)
active_bowtie <- app_env$active_step_jobs_from_jobs(fake_jobs, "Bowtie2")
assert(identical(sort(active_bowtie$job_id), c("101", "102")), "active-job filtering excludes completed, cancelled, and other-step jobs")
assert(identical(unname(app_env$active_step_sample_choices(fake_jobs, "Bowtie2")), c("A1", "A2")), "cancellation choices contain only samples with active jobs")
assert(identical(app_env$filter_active_jobs_by_samples(active_bowtie, "A2")$job_id, "102"), "selected-sample cancellation resolves only the requested active job")
assert(inherits(app_env$active_jobs_modal_table(active_bowtie), "shiny.tag"), "active-job cancellation summary renders locally")

assay_jobs <- data.frame(
  step = c("STAR", "Bowtie2", "SEACR", "MACS2 Peaks"),
  sample = c("rna_sample", "atac_sample", "cutrun_sample", "chip_sample"),
  job_id = c("201", "202", "203", "204"),
  slurm_state = rep("RUNNING", 4),
  stringsAsFactors = FALSE
)
for (step in assay_jobs$step) {
  expected <- assay_jobs$sample[assay_jobs$step == step]
  assert(identical(unname(app_env$active_step_sample_choices(assay_jobs, step)), expected), paste(step, "supports active sample cancellation"))
}

sample_aware_submitters <- c(
  "submit_cutadapt_jobs", "submit_fastqc_jobs", "submit_star_jobs", "submit_featurecounts_jobs", "submit_rsem_jobs", "submit_kallisto_jobs",
  "submit_cutrun_bowtie2_jobs", "submit_cutrun_seacr_jobs", "submit_cutrun_macs2_jobs",
  "submit_atac_bowtie2_jobs", "submit_atac_macs2_jobs", "submit_chip_bowtie2_jobs", "submit_chip_macs2_jobs"
)
for (function_name in sample_aware_submitters) {
  assert("samples" %in% names(formals(app_env[[function_name]])), paste(function_name, "accepts explicit sample selection"))
}

runner_test_root <- file.path(root, "rna-runner-submit")
dir.create(file.path(runner_test_root, "fastq"), recursive = TRUE)
runner_design_path <- file.path(runner_test_root, "design_matrix.txt")
runner_design <- data.frame(
  sample = "rna1", treatment = "control",
  filename = "rna1_R1.fastq.gz,rna1_R2.fastq.gz", stringsAsFactors = FALSE
)
write.table(runner_design, runner_design_path, sep = "\t", row.names = FALSE, quote = FALSE)
writeBin(as.raw(rep(1:100, 20)), file.path(runner_test_root, "fastq", "rna1_R1.fastq.gz"))
writeBin(as.raw(rep(1:100, 20)), file.path(runner_test_root, "fastq", "rna1_R2.fastq.gz"))
runner_project <- list(
  id = "runner-rna", name = "runner-rna", analysis_key = "rna", analysis = "RNA-seq",
  design_matrix_path = runner_design_path, data_dir = file.path(runner_test_root, "data"),
  results_root = runner_test_root, fastq_dir = file.path(runner_test_root, "fastq"),
  fastq_dirs = file.path(runner_test_root, "fastq"), paired_end = TRUE,
  genome = "human", genome_version = "human_gencode50"
)
dir.create(runner_project$data_dir, recursive = TRUE)
fake_gtf <- file.path(runner_test_root, "genes.gtf")
fake_strand_bed <- file.path(runner_test_root, "strand.bed")
writeLines("annotation", fake_gtf)
writeLines("chr1\t1\t10\tgene1", fake_strand_bed)
captured_submissions <- list()
original_genome_resources <- app_env$genome_resources
original_submit_sbatch <- app_env$submit_sbatch
original_submit_sbatch_wrap <- app_env$submit_sbatch_wrap
original_submit_featurecounts_matrix_job <- app_env$submit_featurecounts_matrix_job
app_env$genome_resources <- function(project) list(
  star_index = file.path(runner_test_root, "star-index"), label = "test-reference",
  gtf = fake_gtf, strand_bed = fake_strand_bed,
  kallisto_index = file.path(runner_test_root, "transcripts.idx"),
  rsem_index = file.path(runner_test_root, "rsem")
)
app_env$submit_sbatch <- function(project, step, script, args, log_name, input_mode = "", sample = "", target = "", reference = "", dependency_ids = character(0), dependency_condition = "afterok") {
  captured_submissions[[length(captured_submissions) + 1L]] <<- list(step = step, script = script, args = args, target = target)
  paste("Submitted", step, sample)
}
app_env$submit_sbatch_wrap <- function(...) "Matrix build captured"
app_env$submit_featurecounts_matrix_job <- function(project, feature = "gene_name", dependency_ids = character(0)) "Matrix build captured"
on.exit({
  app_env$genome_resources <- original_genome_resources
  app_env$submit_sbatch <- original_submit_sbatch
  app_env$submit_sbatch_wrap <- original_submit_sbatch_wrap
  app_env$submit_featurecounts_matrix_job <- original_submit_featurecounts_matrix_job
}, add = TRUE)

invisible(app_env$submit_star_jobs(runner_project, trimmed = FALSE, samples = "rna1"))
star_submission <- captured_submissions[[length(captured_submissions)]]
expected_star_runner <- file.path(app_env$SCRIPTS_DIR, "STAR", "star_PE.sh")
assert(identical(tail(star_submission$args, 1), expected_star_runner), "STAR submission passes an absolute runner path to SLURM")

star_dir <- file.path(runner_project$data_dir, "star", "rna1")
dir.create(star_dir, recursive = TRUE, showWarnings = FALSE)
writeBin(as.raw(rep(1:100, 30)), file.path(star_dir, "rna1Aligned.sortedByCoord.out.bam"))
invisible(app_env$submit_featurecounts_jobs(runner_project, feature = "gene_name", samples = "rna1"))
feature_submission <- captured_submissions[[length(captured_submissions)]]
expected_feature_runner <- file.path(app_env$SCRIPTS_DIR, "featureCounts", "featurecounts_PE.sh")
assert(identical(tail(feature_submission$args, 1), expected_feature_runner), "featureCounts submission passes an absolute runner path to SLURM")

app_env$genome_resources <- function(project) list(
  star_index = file.path(runner_test_root, "maize-star-index"), label = "maize-test-reference",
  gtf = fake_gtf
)
runner_project$genome <- "maize"
runner_project$genome_version <- "maize_b73_nam5"
invisible(app_env$submit_featurecounts_jobs(runner_project, feature = "gene_name", samples = "rna1"))
maize_feature_submission <- captured_submissions[[length(captured_submissions)]]
assert(
  length(maize_feature_submission$args) == 7L &&
    identical(maize_feature_submission$args[[3]], "gene_id") &&
    identical(maize_feature_submission$args[[5]], "none") &&
    identical(maize_feature_submission$args[[6]], runner_project$name) &&
    identical(maize_feature_submission$args[[7]], expected_feature_runner),
  "maize featureCounts preserves every SLURM argument, forces gene_id, and explicitly requests unstranded fallback without a strand BED"
)
runner_project$genome <- "human"
runner_project$genome_version <- "human_gencode50"
app_env$genome_resources <- function(project) list(
  star_index = file.path(runner_test_root, "star-index"), label = "test-reference",
  gtf = fake_gtf, strand_bed = fake_strand_bed,
  kallisto_index = file.path(runner_test_root, "transcripts.idx"),
  rsem_index = file.path(runner_test_root, "rsem")
)

invisible(app_env$submit_kallisto_jobs(runner_project, trimmed = FALSE, samples = "rna1"))
kallisto_submission <- captured_submissions[[length(captured_submissions)]]
expected_kallisto_runner <- file.path(app_env$SCRIPTS_DIR, "Kallisto", "kallisto_PE.sh")
assert(identical(tail(kallisto_submission$args, 1), expected_kallisto_runner), "paired-end Kallisto submission passes an absolute runner path to SLURM")
assert(identical(kallisto_submission$target, file.path(runner_project$data_dir, "kallisto", "rna1", "abundance.tsv")), "Kallisto completion target matches its output directory")

writeBin(as.raw(rep(1:100, 30)), file.path(star_dir, "rna1Aligned.toTranscriptome.out.bam"))
invisible(app_env$submit_rsem_jobs(runner_project, feature = "gene_id", samples = "rna1"))
rsem_submission <- captured_submissions[[length(captured_submissions)]]
expected_rsem_runner <- file.path(app_env$SCRIPTS_DIR, "RSEM", "RSEM_PE.sh")
assert(identical(tail(rsem_submission$args, 1), expected_rsem_runner), "paired-end RSEM submission passes an absolute runner path to SLURM")

runner_project$paired_end <- FALSE
runner_design$filename <- "rna1_R1.fastq.gz"
write.table(runner_design, runner_design_path, sep = "\t", row.names = FALSE, quote = FALSE)
unlink(star_dir, recursive = TRUE)
invisible(app_env$submit_star_jobs(runner_project, trimmed = FALSE, samples = "rna1"))
star_se_submission <- captured_submissions[[length(captured_submissions)]]
expected_star_se_runner <- file.path(app_env$SCRIPTS_DIR, "STAR", "star_SE.sh")
assert(identical(tail(star_se_submission$args, 1), expected_star_se_runner), "single-end STAR submission passes its absolute runner path to SLURM")
assert(length(star_se_submission$args) == 5L, "single-end STAR submission omits the unused R2 argument")

dir.create(star_dir, recursive = TRUE, showWarnings = FALSE)
writeBin(as.raw(rep(1:100, 30)), file.path(star_dir, "rna1Aligned.sortedByCoord.out.bam"))
invisible(app_env$submit_featurecounts_jobs(runner_project, feature = "gene_id", samples = "rna1"))
feature_se_submission <- captured_submissions[[length(captured_submissions)]]
expected_feature_se_runner <- file.path(app_env$SCRIPTS_DIR, "featureCounts", "featurecounts_SE.sh")
assert(identical(tail(feature_se_submission$args, 1), expected_feature_se_runner), "single-end featureCounts submission passes its absolute runner path to SLURM")

invisible(app_env$submit_kallisto_jobs(runner_project, trimmed = FALSE, samples = "rna1"))
kallisto_se_submission <- captured_submissions[[length(captured_submissions)]]
expected_kallisto_se_runner <- file.path(app_env$SCRIPTS_DIR, "Kallisto", "kallisto_SE.sh")
assert(identical(tail(kallisto_se_submission$args, 1), expected_kallisto_se_runner), "single-end Kallisto submission passes its absolute runner path to SLURM")
assert(length(kallisto_se_submission$args) == 5L, "single-end Kallisto submission omits the unused R2 argument")
assert(identical(kallisto_se_submission$target, file.path(runner_project$data_dir, "kallisto", "rna1", "abundance.tsv")), "single-end Kallisto writes to the App completion target")

writeBin(as.raw(rep(1:100, 30)), file.path(star_dir, "rna1Aligned.toTranscriptome.out.bam"))
invisible(app_env$submit_rsem_jobs(runner_project, feature = "gene_id", samples = "rna1"))
rsem_se_submission <- captured_submissions[[length(captured_submissions)]]
expected_rsem_se_runner <- file.path(app_env$SCRIPTS_DIR, "RSEM", "RSEM_SE.sh")
assert(identical(tail(rsem_se_submission$args, 1), expected_rsem_se_runner), "single-end RSEM submission passes its absolute runner path to SLURM")

app_env$genome_resources <- original_genome_resources
app_env$submit_sbatch <- original_submit_sbatch
app_env$submit_sbatch_wrap <- original_submit_sbatch_wrap
app_env$submit_featurecounts_matrix_job <- original_submit_featurecounts_matrix_job

matrix_root <- file.path(root, "quant-matrix-build")
kallisto_dir <- file.path(matrix_root, "kallisto")
rsem_dir <- file.path(matrix_root, "rsem")
counts_dir <- file.path(matrix_root, "counts")
for (sample in c("sample1", "sample2")) {
  dir.create(file.path(kallisto_dir, sample), recursive = TRUE, showWarnings = FALSE)
  write.table(data.frame(
    target_id = c("tx1", "tx2"), length = c(100, 200), eff_length = c(80, 180),
    est_counts = c(5, 10) + match(sample, c("sample1", "sample2")),
    tpm = c(20, 30) + match(sample, c("sample1", "sample2"))
  ), file.path(kallisto_dir, sample, "abundance.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  dir.create(file.path(rsem_dir, sample), recursive = TRUE, showWarnings = FALSE)
  write.table(data.frame(
    gene_id = c("gene1", "gene2"), expected_count = c(5, 10), TPM = c(20, 30), FPKM = c(8, 12)
  ), file.path(rsem_dir, sample, paste0(sample, ".genes.results")), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(
    transcript_id = c("tx1", "tx2"), gene_id = c("gene1", "gene2"),
    expected_count = c(5, 10), TPM = c(20, 30), FPKM = c(8, 12)
  ), file.path(rsem_dir, sample, paste0(sample, ".isoforms.results")), sep = "\t", row.names = FALSE, quote = FALSE)
}
kallisto_matrix_script <- app_env$write_quant_matrix_script(runner_project, "kallisto")
assert(system2("Rscript", c(kallisto_matrix_script, kallisto_dir, counts_dir)) == 0L, "Kallisto matrix builder executes on synthetic quantifications")
kallisto_matrix <- read.delim(file.path(counts_dir, "kallisto_tpm_matrix.txt"), check.names = FALSE)
assert(identical(names(kallisto_matrix), c("target_id", "sample1", "sample2")) && NROW(kallisto_matrix) == 2L, "Kallisto TPM matrix combines all samples")

rsem_matrix_script <- app_env$write_quant_matrix_script(runner_project, "rsem")
assert(system2("Rscript", c(rsem_matrix_script, rsem_dir, counts_dir)) == 0L, "RSEM matrix builder executes on synthetic quantifications")
rsem_matrix <- read.delim(file.path(counts_dir, "rsem_tpm_matrix.txt"), check.names = FALSE)
assert(identical(names(rsem_matrix), c("gene_id", "sample1", "sample2")) && NROW(rsem_matrix) == 2L, "RSEM TPM matrix combines all samples")
assert(file.exists(file.path(counts_dir, "rsem_isoform_tpm_matrix.txt")), "RSEM isoform TPM matrix is generated")
assert(identical(app_env$requested_sample_subset(atac_project, c("A1", "A2"), "A2", "test step"), "A2"), "unchecked samples are excluded from submission")
cutrun_example_project <- chip_project
cutrun_example_project$analysis_key <- "cutrun"
cutrun_example_project$analysis <- "CUT&RUN"
cutrun_example_project$design_matrix_path <- file.path(app_env$example_dataset_paths("cutrun")$design_dir, "design_matrix.txt")
cutrun_targets <- app_env$pipeline_step_sample_candidates(cutrun_example_project, targets_only = TRUE)
assert(length(cutrun_targets) == 4L && !any(grepl("IgG", cutrun_targets)), "CUT&RUN peak-step selectors contain targets but not controls")

seacr_selector_root <- file.path(root, "cutrun-seacr-selector")
dir.create(seacr_selector_root, recursive = TRUE)
seacr_selector_design <- data.frame(
  sample = c("S1", "S2"), include = c(TRUE, TRUE), cell_type = c("AKP", "AKP"),
  mark = c("Creb", "Creb"), target = c("Creb", "Creb"), target_class = c("tf_or_other", "tf_or_other"),
  seacr_stringency = c("auto", "auto"), condition = c("AA", "AA"), replicate = c(1, 2),
  control_sample = c("", ""), filename = c("S1.fastq.gz", "S2.fastq.gz"), stringsAsFactors = FALSE
)
seacr_selector_design_path <- file.path(seacr_selector_root, "design_matrix.txt")
write.table(seacr_selector_design, seacr_selector_design_path, sep = "\t", row.names = FALSE, quote = FALSE)
seacr_selector_project <- cutrun_example_project
seacr_selector_project$id <- "cutrun-seacr-selector"
seacr_selector_project$name <- "cutrun-seacr-selector"
seacr_selector_project$data_dir <- file.path(seacr_selector_root, "data")
seacr_selector_project$design_matrix_path <- seacr_selector_design_path
for (spec in list(c("S1", "norm"), c("S2", "non"))) {
  sample <- spec[[1]]
  norm <- spec[[2]]
  peak <- app_env$cutrun_seacr_peak_path(seacr_selector_project, sample, norm, "stringent")
  summary <- app_env$cutrun_seacr_summary_path(seacr_selector_project, sample, norm, "stringent")
  dir.create(dirname(peak), recursive = TRUE, showWarnings = FALSE)
  file.create(peak)
  writeLines(c(paste0("normalization\t", norm), "stringency\tstringent", "peak_count\t0"), summary)
}
assert(identical(app_env$completed_cutrun_seacr_samples(seacr_selector_project, c("S1", "S2"), "norm", "stringent"), "S1"), "SEACR selector completion is specific to norm/stringency and accepts a completed zero-peak result")
assert(identical(app_env$completed_cutrun_seacr_samples(seacr_selector_project, c("S1", "S2"), "non", "stringent"), "S2"), "SEACR non completion does not hide samples from the norm selector")
shared_seacr_peak <- app_env$cutrun_seacr_peak_path(seacr_selector_project, "S1", "norm", "relaxed", "raw")
shared_seacr_summary <- app_env$cutrun_seacr_summary_path(seacr_selector_project, "S1", "norm", "relaxed", "raw")
dir.create(dirname(shared_seacr_peak), recursive = TRUE, showWarnings = FALSE)
writeLines(c("chr1\t100\t220\t25", "chr1\t400\t520\t15"), shared_seacr_peak)
writeLines(c("normalization\tnorm", "stringency\trelaxed", "peak_count\t2"), shared_seacr_summary)
shared_macs_dir <- file.path(seacr_selector_project$data_dir, "macs2", "S1")
dir.create(shared_macs_dir, recursive = TRUE, showWarnings = FALSE)
shared_macs_peak <- file.path(shared_macs_dir, "S1_peaks.narrowPeak")
writeLines(c("chr1\t150\t250\tpeak1\t100", "chr1\t600\t700\tpeak2\t100"), shared_macs_peak)
writeLines(c("sample\tS1", "qval\t0.01", "peak_type\tnarrow", paste0("peak_file\t", shared_macs_peak), "peak_count\t2"), file.path(shared_macs_dir, "S1_macs2_summary.txt"))
shared_sources <- app_env$cutrun_peak_source_catalog(seacr_selector_project)
shared_seacr_id <- "seacr_raw_norm_relaxed"
shared_macs_id <- "macs2_narrow_q_0_01"
assert(all(c(shared_seacr_id, shared_macs_id) %in% shared_sources$source_id), "CUT&RUN peak-overlap sources distinguish each caller and concrete setting")
assert(
  identical(app_env$cutrun_peak_source_file(seacr_selector_project, shared_seacr_id, "S1"), normalizePath(shared_seacr_peak)) &&
    identical(app_env$cutrun_peak_source_file(seacr_selector_project, shared_macs_id, "S1"), normalizePath(shared_macs_peak)),
  "peak-overlap source selection resolves the correct per-sample SEACR and MACS2 inputs"
)
assert(
  identical(
    unname(app_env$cutrun_peak_source_files(seacr_selector_project, shared_macs_id, c("S1", "S2"), sources = shared_sources)),
    c(normalizePath(shared_macs_peak), "")
  ),
  "peak-overlap controls resolve selected source paths in one cached pass rather than rescanning all summaries per sample"
)
shared_overlap_bed <- app_env$cutrun_peak_overlap_bed(seacr_selector_project, shared_seacr_id, shared_macs_id, "S1")
shared_overlap_summary <- app_env$cutrun_peak_overlap_summary_path(seacr_selector_project, shared_seacr_id, shared_macs_id, "S1")
dir.create(dirname(shared_overlap_bed), recursive = TRUE, showWarnings = FALSE)
writeLines(c("chr1\t150\t220", "chr1\t400\t450"), shared_overlap_bed)
shared_overlap_ranking <- sub("\\.bed$", "_ranking.tsv", shared_overlap_bed)
writeLines(c(
  "chrom\tstart\tend\tcombined_evidence_rank\tsource_a_rank\tsource_b_rank\tsource_a_evidence_score\tsource_b_evidence_score",
  "chr1\t400\t450\t5\t3\t2\t15\t12",
  "chr1\t150\t220\t2\t1\t1\t25\t20"
), shared_overlap_ranking)
writeLines(c("sample\tS1", paste0("overlap_name\t", app_env$cutrun_peak_overlap_name(shared_seacr_id, shared_macs_id)), "source_a_peaks\t2", "source_b_peaks\t2", "overlap_peaks\t2", "minimum_reciprocal_overlap\t0", paste0("overlap_bed\t", shared_overlap_bed), paste0("ranking_tsv\t", shared_overlap_ranking)), shared_overlap_summary)
shared_overlap_summary_table <- app_env$cutrun_peak_overlap_summary_table(seacr_selector_project)
assert(NROW(shared_overlap_summary_table) == 1L && shared_overlap_summary_table[["Shared overlap peaks"]][[1]] == "2", "shared peak-overlap output produces a compact per-sample summary table")
diffbind_design <- data.frame(
  sample = c("S1", "S2", "S3", "S4"),
  cell_type = "Model",
  mark = "Creb",
  target_class = "tf_or_other",
  condition = c("A", "A", "B", "B"),
  replicate = c(1, 2, 1, 2),
  control_sample = "",
  filename = paste0("S", 1:4, "_R1.fastq.gz,S", 1:4, "_R2.fastq.gz"),
  stringsAsFactors = FALSE
)
write.table(diffbind_design, seacr_selector_design_path, sep = "\t", row.names = FALSE, quote = FALSE)
for (sample in c("S2", "S3", "S4")) {
  macs_dir <- file.path(seacr_selector_project$data_dir, "macs2", sample)
  dir.create(macs_dir, recursive = TRUE, showWarnings = FALSE)
  macs_peak <- file.path(macs_dir, paste0(sample, "_peaks.narrowPeak"))
  writeLines(c("chr1\t100\t220\tpeak1\t100", "chr1\t400\t520\tpeak2\t100"), macs_peak)
  writeLines(c(paste0("sample\t", sample), "qval\t0.01", "peak_type\tnarrow", paste0("peak_file\t", macs_peak), "peak_count\t2"), file.path(macs_dir, paste0(sample, "_macs2_summary.txt")))
  overlap_bed <- app_env$cutrun_peak_overlap_bed(seacr_selector_project, shared_seacr_id, shared_macs_id, sample)
  dir.create(dirname(overlap_bed), recursive = TRUE, showWarnings = FALSE)
  writeLines(c("chr1\t150\t220", "chr1\t430\t500"), overlap_bed)
  writeLines(c(
    paste0("sample\t", sample),
    paste0("overlap_name\t", app_env$cutrun_peak_overlap_name(shared_seacr_id, shared_macs_id)),
    "source_a_peaks\t2", "source_b_peaks\t2", "overlap_peaks\t2",
    paste0("overlap_bed\t", overlap_bed)
  ), sub("\\.bed$", "_summary.txt", overlap_bed))
}
for (sample in paste0("S", 1:4)) {
  bam <- app_env$cutrun_bowtie2_signal_bam(seacr_selector_project, sample)
  dir.create(dirname(bam), recursive = TRUE, showWarnings = FALSE)
  writeBin(as.raw(seq_len(64)), bam)
}
diffbind_sources <- app_env$cutrun_diffbind_peak_source_catalog(seacr_selector_project)
shared_diffbind_id <- paste0("shared_", app_env$cutrun_peak_overlap_name(shared_seacr_id, shared_macs_id))
assert(all(c(shared_macs_id, shared_diffbind_id) %in% diffbind_sources$source_id), "CUT&RUN DiffBind can select MACS2 or a completed shared-overlap peak source")
assert(
  identical(
    app_env$cutrun_peak_source_file(seacr_selector_project, shared_diffbind_id, "S1", sources = diffbind_sources),
    normalizePath(shared_overlap_bed)
  ),
  "shared-overlap DiffBind source resolves the selected per-sample BED"
)
two_rep_plan <- app_env$cutrun_diffbind_comparison_plan(seacr_selector_project, "A", 1L, shared_diffbind_id, 2L)
assert(NROW(two_rep_plan) == 1L && isTRUE(two_rep_plan$eligible[[1]]) && two_rep_plan$comparison_replicates[[1]] == 2L && two_rep_plan$reference_replicates[[1]] == 2L, "CUT&RUN DiffBind accepts exactly two biological replicates per condition")
diffbind_run_slug <- app_env$cutrun_diffbind_run_slug(shared_diffbind_id, "Model", "Creb", "B", "A")
ready_status <- app_env$cutrun_diffbind_comparison_status(seacr_selector_project, two_rep_plan, jobs = data.frame())
assert(NROW(ready_status) == 1L && identical(ready_status$Status[[1]], "Ready"), "eligible CUT&RUN DiffBind comparisons are shown as ready before submission")
running_status <- app_env$cutrun_diffbind_comparison_status(
  seacr_selector_project,
  two_rep_plan,
  jobs = data.frame(
    step = "Differential Peaks", sample = diffbind_run_slug, job_id = "12345",
    slurm_state = "RUNNING", elapsed = "00:03:12", stringsAsFactors = FALSE
  )
)
assert(
  identical(running_status$Status[[1]], "Running") &&
    identical(running_status$`Job ID`[[1]], "12345") &&
    identical(running_status$Elapsed[[1]], "00:03:12"),
  "CUT&RUN DiffBind status table tracks each comparison's SLURM job independently"
)
diffbind_complete_marker <- file.path(seacr_selector_project$data_dir, "cutrun_diffbind", diffbind_run_slug, "_COMPLETE")
dir.create(dirname(diffbind_complete_marker), recursive = TRUE, showWarnings = FALSE)
writeLines(as.character(Sys.time()), diffbind_complete_marker)
write.table(
  data.frame(
    cell_type = "Model", mark = "Creb", comparison = "B", reference = "A",
    peak_source = shared_diffbind_id, stringsAsFactors = FALSE
  ),
  file.path(dirname(diffbind_complete_marker), "cutrun_diffbind_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
diffbind_result_label <- app_env$cutrun_diffbind_result_label(seacr_selector_project, dirname(diffbind_complete_marker))
assert(
  grepl("^Model — Creb — B vs A — peaks: Shared overlap", diffbind_result_label),
  "CUT&RUN DiffBind result labels lead with the biological comparison and show the peak source last"
)
complete_status <- app_env$cutrun_diffbind_comparison_status(seacr_selector_project, two_rep_plan, jobs = data.frame())
assert(identical(complete_status$Status[[1]], "Complete"), "CUT&RUN DiffBind status uses each comparison's own completion marker")
strict_peak_plan <- app_env$cutrun_diffbind_comparison_plan(seacr_selector_project, "A", 1L, shared_diffbind_id, 3L)
assert(NROW(strict_peak_plan) == 1L && !isTRUE(strict_peak_plan$eligible[[1]]) && grepl("Below 3 peaks", strict_peak_plan$reason[[1]], fixed = TRUE), "CUT&RUN DiffBind excludes comparisons when any selected source has fewer than the requested peaks")
assert(
  !NROW(app_env$cutrun_diffbind_comparison_status(seacr_selector_project, strict_peak_plan, jobs = data.frame())),
  "ineligible CUT&RUN comparisons are not shown as unavailable status rows"
)
resolved_shared_sheet <- app_env$cutrun_diffbind_sample_sheet(seacr_selector_project, "A", 1L, "Model", "Creb", "B", shared_diffbind_id, 2L, "cpm")
resolved_shared <- read.delim(resolved_shared_sheet, check.names = FALSE, stringsAsFactors = FALSE)
assert(
  NROW(resolved_shared) == 4L && all(resolved_shared$PeakCount == 2L) &&
    all(grepl("/peak_overlap/", resolved_shared$Peaks, fixed = TRUE)),
  "resolved CUT&RUN DiffBind sheet records the chosen shared-overlap paths and per-sample peak counts"
)
inferred_project <- seacr_selector_project
inferred_project$data_dir <- file.path(root, "cutrun_inferred_results")
inferred_project$design_matrix_path <- file.path(inferred_project$data_dir, "manifest", "design_matrix.txt")
for (sample in c("AKPS_Creb-AA1", "AKPS_Creb-AA2", "AKPS_Creb-Veh1", "AKPS_Creb-Veh2", "AKPS_IgG_AA", "AKPS_IgG_Veh")) {
  dir.create(file.path(inferred_project$data_dir, "bowtie2", sample), recursive = TRUE, showWarnings = FALSE)
}
inferred_design <- app_env$project_design_df(inferred_project)
assert(
  NROW(inferred_design) == 6L &&
    all(c("AA", "Veh") %in% inferred_design$condition) &&
    all(c("1", "2") %in% inferred_design$replicate) &&
    sum(tolower(inferred_design$mark) == "igg") == 2L,
  "completed CUT&RUN outputs without a saved manifest infer sample, condition, replicate, mark, and control metadata from sample folders"
)
shared_overlap_navigation <- app_env$cutrun_individual_peak_navigation(shared_overlap_bed, max_peaks = 2L)
assert(
  identical(unname(shared_overlap_navigation$peaks[[1]]), "chr1:151-220") && grepl("combined rank 2", names(shared_overlap_navigation$peaks)[[1]], fixed = TRUE),
  "shared peak-overlap browser navigation uses the precomputed lowest combined caller rank first"
)
shared_project_summary <- app_env$cutrun_seacr_peak_summary_table(seacr_selector_project)
assert(any(grepl("Shared Peaks:", names(shared_project_summary), fixed = TRUE)), "project peak summary automatically adds each shared-overlap peak count column")
app_source_text <- paste(readLines(file.path(repo_root, "app.R"), warn = FALSE), collapse = "\n")
assert(grepl("Select all samples", app_source_text, fixed = TRUE) && grepl("Clear selection", app_source_text, fixed = TRUE), "sample-level step selectors expose select-all and clear controls")
assert(grepl("cslRestoreToolPanels", app_source_text, fixed = TRUE) && grepl("server = TRUE", app_source_text, fixed = TRUE), "pipeline panels preserve open state and tables use server-side rendering")
assert(grepl("Track normalization to generate", app_source_text, fixed = TRUE) && grepl("Generate selected signal tracks", app_source_text, fixed = TRUE), "CUT&RUN app exposes normalization regeneration from existing aligned BAMs")
assert(grepl("run_cutrun_peak_overlap", app_source_text, fixed = TRUE) && grepl("Shared Peaks", app_source_text, fixed = TRUE), "CUT&RUN app exposes selectable shared peak-overlap runs and results")

assert(system2("bash", c("-n", shQuote(file.path(repo_root, "run_codespringweb.sh")))) == 0L, "CodeSpringApp launcher shell syntax is valid")
assert(system2("bash", c("-n", shQuote(file.path(lab_root, "scripts_DoNotTouch", "CUTRUN", "cutrun_peak_overlap.sh")))) == 0L, "CUT&RUN overlap runner shell syntax is valid")
assert(system2("bash", c("-n", shQuote(file.path(lab_root, "scripts_DoNotTouch", "CUTRUN", "qsub_cutrun_peak_overlap.sh")))) == 0L, "CUT&RUN overlap submission wrapper shell syntax is valid")

bad_q <- app_env$submit_atac_macs2_jobs(atac_project, "not-a-number", "A1")
assert(grepl("q-value must be", bad_q), "invalid ATAC MACS2 q-value rejected before submission")
assert(grepl("two different", app_env$submit_atac_diffbind_job(atac_project, "condition", "A", "A")), "identical ATAC DiffBind conditions rejected")
assert(grepl("two different", app_env$submit_chip_diffbind_job(chip_project, "condition", "A", "A")), "identical ChIP DiffBind conditions rejected")

comparison_dir <- file.path(root, "diffbind", "B_vs_A")
dir.create(comparison_dir, recursive = TRUE)
legacy_result <- file.path(comparison_dir, "DifferentialPeaks_B_vs_A_ref.txt")
writeLines("Fold\tFDR\n1\t0.01", legacy_result)
assert(identical(app_env$peak_diffbind_status(atac_project), "Complete"), "legacy DiffBind comparison remains recognized")
writeLines(character(0), legacy_result)
assert(!app_env$diffbind_comparison_complete(comparison_dir), "empty legacy result is not accepted")
writeLines("Fold\tFDR\n1\t0.01", legacy_result)
writeLines("status\trunning", file.path(comparison_dir, "_RUN_STARTED"))
assert(identical(app_env$peak_diffbind_status(atac_project), "Likely failed"), "partial DiffBind output is not accepted")
assert(!app_env$diffbind_comparison_complete(comparison_dir), "started DiffBind comparison hidden from Results Explorer")
active_jobs <- data.frame(
  step = "Differential Peaks", slurm_state = "RUNNING", sample = basename(comparison_dir),
  target = file.path(comparison_dir, "_COMPLETE"), stringsAsFactors = FALSE
)
assert(app_env$diffbind_comparison_active(atac_project, comparison_dir, jobs = active_jobs), "active DiffBind comparison recognized")
unlink(file.path(comparison_dir, "_RUN_STARTED"))
writeLines("status\tcomplete", file.path(comparison_dir, "_COMPLETE"))
assert(app_env$diffbind_comparison_complete(comparison_dir), "final DiffBind marker accepted")
writeLines(c(
  "seqnames\tstart\tend\tFold\tp.value\tFDR",
  "chr1\t101\t220\t2.5\t0.0002\t0.004",
  "chr1\t501\t650\t-1.8\t0.003\t0.02"
), legacy_result)

differential_bed <- file.path(comparison_dir, "DifferentialPeaks_B_vs_A_ref.with_stats.bed")
writeLines(c(
  "chr1\t100\t220\tpeak_1|Fold=2.5|p.value=0.0002|FDR=0.004\t2.5",
  "chr1\t500\t650\tpeak_2|Fold=-1.8|p.value=0.003|FDR=0.02\t-1.8"
), differential_bed)
expanded_bed <- app_env$safe_read_result_table(differential_bed)
assert(all(c("Fold", "p.value", "FDR") %in% names(expanded_bed)), "ATAC with-stats BED exposes p-value and FDR columns in the Results Explorer")
comparison_annotation <- file.path(comparison_dir, "DifferentialPeaks_B_vs_A_ref_annotated_with_stats.txt")
writeLines(c(
  "PeakID (cmd=annotatePeaks.pl synthetic.with_stats.bed mm39)\tGene Name\tAnnotation\tDetailed Annotation\tDistance to TSS\tGene Description",
  "chr1:101-220|Fold=2.5|p.value=0.0002|FDR=0.004\tGeneA\tPromoter\tpromoter-TSS\t15\tSynthetic gene A",
  "chr1:501-650|Fold=-1.8|p.value=0.003|FDR=0.02\tGeneB\tIntron\tintron (GeneB)\t850\tSynthetic gene B"
), comparison_annotation)
navigation <- app_env$genome_browser_comparison_navigation(comparison_dir)
assert(all(c("GeneA", "GeneB") %in% unname(navigation$genes)), "genome browser offers annotated differential-peak genes")
assert(identical(unname(navigation$genes), c("GeneA", "GeneB")), "genome-browser differential genes are sorted alphabetically")
assert(length(navigation$peaks) == 2L && grepl("1. chr1:101-220", names(navigation$peaks)[[1]], fixed = TRUE), "genome browser offers ranked searchable differential-peak intervals")
assert(grepl("GeneA", names(navigation$peaks)[[1]], fixed = TRUE), "top differential-peak selector includes HOMER gene-name context")
assert(identical(unname(navigation$peaks)[[1]], "chr1:101-220"), "selected differential peak inserts its exact genomic interval into the locus search")
comparison_direction_navigation <- app_env$genome_browser_comparison_navigation(comparison_dir, direction = "comparison")
reference_direction_navigation <- app_env$genome_browser_comparison_navigation(comparison_dir, direction = "reference")
assert(
  identical(unname(comparison_direction_navigation$peaks), "chr1:101-220") &&
    identical(unname(comparison_direction_navigation$genes), "GeneA"),
  "positive DiffBind Fold filtering retains only peaks higher in the comparison condition"
)
assert(
  identical(unname(reference_direction_navigation$peaks), "chr1:501-650") &&
    identical(unname(reference_direction_navigation$genes), "GeneB"),
  "negative DiffBind Fold filtering retains only peaks higher in the reference condition"
)
comparison_direction_bed <- app_env$genome_browser_filtered_differential_bed(differential_bed, "comparison")
reference_direction_bed <- app_env$genome_browser_filtered_differential_bed(differential_bed, "reference")
assert(NROW(read.table(comparison_direction_bed, sep = "\t")) == 1L && read.table(comparison_direction_bed, sep = "\t")[[5]][[1]] > 0, "comparison-direction browser BED contains only positive-Fold peaks")
assert(NROW(read.table(reference_direction_bed, sep = "\t")) == 1L && read.table(reference_direction_bed, sep = "\t")[[5]][[1]] < 0, "reference-direction browser BED contains only negative-Fold peaks")
ranking_dir <- file.path(root, "diffbind", "ranking_test")
dir.create(ranking_dir, recursive = TRUE)
ranking_rows <- vapply(seq_len(205L), function(i) {
  paste("chr2", i * 100L, i * 100L + 50L, paste0("peak_", i, "|Fold=", i / 10, "|p.value=", (206L - i) / 10000, "|FDR=", (206L - i) / 1000), i / 10, sep = "\t")
}, character(1))
writeLines(ranking_rows, file.path(ranking_dir, "DifferentialPeaks_ranked.with_stats.bed"))
ranking_annotations <- c("PeakID\tGene Name\tFDR\tp.value\tFold", vapply(seq_len(205L), function(i) {
  paste0("peak_", i, "\tGene", i, "\t", (206L - i) / 1000, "\t", (206L - i) / 10000, "\t", i / 10)
}, character(1)))
writeLines(ranking_annotations, file.path(ranking_dir, "DifferentialPeaks_ranked_annotated_with_stats.txt"))
ranked_navigation <- app_env$genome_browser_comparison_navigation(ranking_dir, max_peaks = 200L)
assert(length(ranked_navigation$peaks) == 200L, "genome browser caps the differential-peak selector at 200 entries")
assert(grepl("1. chr2:20501-20550", names(ranked_navigation$peaks)[[1]], fixed = TRUE), "strongest differential peak is first in the selector")
assert(length(ranked_navigation$genes) <= 200L && "Gene205" %in% unname(ranked_navigation$genes), "gene selector uses only top annotated DiffBind genes")
writeLines(c(
  "PeakID\tGene Name",
  "chr9:1-10|Fold=1\tLeanGene"
), file.path(comparison_dir, "DifferentialPeaks_Z_vs_A_ref_annotated_with_stats.txt"))
differential_table <- app_env$differential_accessibility_result_table(comparison_dir)
assert(identical(differential_table[["Genomic interval"]], c("chr1:101-220", "chr1:501-650")), "differential accessibility table includes complete genomic intervals")
assert(identical(differential_table[["Gene Name"]], c("GeneA", "GeneB")), "differential accessibility table prefers the HOMER-annotated result and displays Gene Name")
assert(all(c("Fold", "p.value", "FDR") %in% names(differential_table)), "annotated differential accessibility table retains expanded DiffBind statistics")
assert(all(c("Detailed Annotation", "Distance to TSS", "Gene Description") %in% names(differential_table)), "differential accessibility table chooses the most detailed annotated result available")
chip_sheet_dir <- file.path(root, "manifest", "chip_diffbind", basename(comparison_dir))
dir.create(chip_sheet_dir, recursive = TRUE, showWarnings = FALSE)
comparison_samples <- data.frame(
  SampleID = c("A1", "A2", "B1", "B2"), Condition = c("A", "A", "B", "B"), Replicate = c(1, 2, 1, 2),
  bamReads = "synthetic.bam", Peaks = "synthetic.narrowPeak", PeakCaller = "narrowpeak", stringsAsFactors = FALSE
)
write.table(comparison_samples, file.path(chip_sheet_dir, "chip_diffbind_samples.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
for (sample in c("A2", "B1", "B2")) {
  sample_signal_dir <- file.path(root, "bowtie2", sample)
  dir.create(sample_signal_dir, recursive = TRUE, showWarnings = FALSE)
  writeBin(as.raw(seq_len(64)), file.path(sample_signal_dir, paste0(sample, "Aligned.sortedByCoord_removeDup.out.bw")))
}
comparison_catalog <- app_env$genome_browser_comparison_catalog(chip_project)
assert(NROW(comparison_catalog) == 1L, "genome browser discovers a completed differential comparison")
assert(identical(comparison_catalog$samples[[1]], c("B1", "B2", "A1", "A2")), "comparison browser places the experimental condition above the reference condition")
assert(identical(comparison_catalog$reference_condition[[1]], "A"), "comparison browser records the DiffBind reference condition")
assert(
  identical(unname(app_env$genome_browser_comparison_condition_labels(comparison_catalog)), c("B", "A")),
  "differential direction controls use the actual comparison and reference condition labels"
)
assert(identical(comparison_catalog$differential_bed[[1]], normalizePath(differential_bed)), "comparison browser selects the differential BED annotation")
comparison_tracks <- app_env$genome_browser_preferred_signal_rows(
  chip_project, app_env$genome_browser_track_catalog(chip_project),
  comparison_catalog$samples[[1]], comparison_catalog$sample_metadata[[1]]
)
assert(NROW(comparison_tracks) == 4L && identical(comparison_tracks$sample, c("B1", "B2", "A1", "A2")), "comparison browser loads reference signal tracks below experimental tracks")

cutrun_comparison_dir <- file.path(root, "cutrun_diffbind", "Creb", "B_vs_A")
dir.create(cutrun_comparison_dir, recursive = TRUE, showWarnings = FALSE)
writeLines("chr1\t100\t220\tpeak_1\t2.5", file.path(cutrun_comparison_dir, "significant_differential_peaks.bed"))
write.table(
  data.frame(
    seqnames = c("chr1", "chr2"), start = c(101, 501), end = c(220, 650),
    Fold = c(2.5, -1.8), p.value = c(0.0002, 0.003), FDR = c(0.004, 0.02),
    stringsAsFactors = FALSE
  ),
  file.path(cutrun_comparison_dir, "all_differential_peaks.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
writeLines(c(
  "PeakID\tGene Name\tFold\tp.value\tFDR",
  "chr1:101-220|Fold=2.5|p.value=0.0002|FDR=0.004\tCutGeneA\t2.5\t0.0002\t0.004",
  "chr2:501-650|Fold=-1.8|p.value=0.003|FDR=0.02\tCutGeneB\t-1.8\t0.003\t0.02"
), file.path(cutrun_comparison_dir, "all_differential_peaks_annotated_with_stats.txt"))
comparison_samples$normalization_mode <- "spikein"
write.table(comparison_samples, file.path(cutrun_comparison_dir, "diffbind_sample_sheet.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cutrun_browser_project <- chip_project
cutrun_browser_project$analysis_key <- "cutrun"
cutrun_browser_project$analysis <- "CUT&RUN"
cutrun_comparisons <- app_env$genome_browser_comparison_catalog(cutrun_browser_project)
assert(NROW(cutrun_comparisons) == 1L && basename(cutrun_comparisons$id[[1]]) == "B_vs_A", "genome browser discovers nested CUT&RUN differential comparisons")
assert(identical(cutrun_comparisons$samples[[1]], c("B1", "B2", "A1", "A2")), "CUT&RUN genome browser places the reference condition below the comparison")
cutrun_navigation <- app_env$genome_browser_comparison_navigation(cutrun_comparison_dir, project = cutrun_browser_project)
assert(
  length(cutrun_navigation$peaks) == 2L &&
    identical(unname(cutrun_navigation$peaks)[[1]], "chr1:101-220") &&
    "CutGeneA" %in% unname(cutrun_navigation$genes),
  "CUT&RUN comparison browser ranks all differential regions and exposes annotated gene navigation"
)

atac_browser_root <- file.path(root, "atac_browser_case")
atac_browser_samples <- paste0(rep(c("Control", "Treated"), each = 3), rep(1:3, 2))
atac_browser_design <- data.frame(
  sample = atac_browser_samples, condition = rep(c("Control", "Treated"), each = 3), replicate = rep(1:3, 2),
  filename = paste0(atac_browser_samples, "_R1.fastq.gz"), stringsAsFactors = FALSE
)
dir.create(atac_browser_root, recursive = TRUE, showWarnings = FALSE)
atac_browser_design_path <- file.path(atac_browser_root, "design_matrix.txt")
write.table(atac_browser_design, atac_browser_design_path, sep = "\t", row.names = FALSE, quote = FALSE)
atac_manifest_dir <- file.path(atac_browser_root, "manifest", "atac_diffbind", "all_samples", "condition")
dir.create(atac_manifest_dir, recursive = TRUE, showWarnings = FALSE)
write.table(atac_browser_design, file.path(atac_manifest_dir, "design_matrix.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
atac_comparison_dir <- file.path(atac_browser_root, "diffbind", "Treated_vs_Control")
dir.create(atac_comparison_dir, recursive = TRUE, showWarnings = FALSE)
writeLines("chr1\t200\t320\tpeak_1\t3", file.path(atac_comparison_dir, "DifferentialPeaks_Treated_vs_Control_ref.with_stats.bed"))
for (sample in atac_browser_samples) {
  sample_signal_dir <- file.path(atac_browser_root, "bowtie2", sample)
  dir.create(sample_signal_dir, recursive = TRUE, showWarnings = FALSE)
  writeBin(as.raw(seq_len(64)), file.path(sample_signal_dir, paste0(sample, "Aligned.sortedByCoord_removeDup.out.bw")))
}
atac_browser_project <- list(
  id = "atac-browser", name = "atac-browser", analysis_key = "atac", analysis = "ATAC-seq",
  design_matrix_path = atac_browser_design_path, data_dir = atac_browser_root, results_root = root,
  fastq_dir = atac_browser_root, fastq_dirs = atac_browser_root, paired_end = TRUE, genome = "mouse"
)
atac_comparisons <- app_env$genome_browser_comparison_catalog(atac_browser_project)
expected_atac_track_order <- c("Treated1", "Treated2", "Treated3", "Control1", "Control2", "Control3")
assert(NROW(atac_comparisons) == 1L && identical(atac_comparisons$samples[[1]], expected_atac_track_order), "existing ATAC comparisons place vehicle/control tracks below treatment tracks")

for (scope in c("AKP", "AKPS")) {
  scoped_manifest_dir <- file.path(atac_browser_root, "manifest", "atac_diffbind", paste0("cell_type_", scope), "condition")
  dir.create(scoped_manifest_dir, recursive = TRUE, showWarnings = FALSE)
  scoped_samples <- paste0(scope, "_", rep(c("AA", "Veh"), each = 2L), rep(1:2, 2L))
  scoped_design <- data.frame(
    sample = scoped_samples,
    condition = rep(c("AA", "Veh"), each = 2L),
    filename = paste0(scoped_samples, "_R1.fastq.gz"),
    stringsAsFactors = FALSE
  )
  write.table(scoped_design, file.path(scoped_manifest_dir, "design_matrix.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  for (sample in scoped_samples) {
    sample_signal_dir <- file.path(atac_browser_root, "bowtie2", sample)
    dir.create(sample_signal_dir, recursive = TRUE, showWarnings = FALSE)
    writeBin(as.raw(seq_len(64)), file.path(sample_signal_dir, paste0(sample, "Aligned.sortedByCoord_removeDup.out.bw")))
  }
}
akps_comparison_dir <- file.path(atac_browser_root, "diffbind", "AKPS_AA_vs_AKPS_Veh")
dir.create(akps_comparison_dir, recursive = TRUE, showWarnings = FALSE)
writeLines("chr1\t300\t420\tpeak_akps\t4", file.path(akps_comparison_dir, "DifferentialPeaks_AA_vs_Veh_ref.with_stats.bed"))
akps_manifest_sheet <- app_env$genome_browser_atac_manifest_sheet(atac_browser_project, akps_comparison_dir)
assert(NROW(akps_manifest_sheet) == 4L && all(startsWith(akps_manifest_sheet$sample, "AKPS_")), "ATAC comparison scope matches AKPS exactly and never falls back to the AKP prefix")
external_atac_project <- atac_browser_project
external_atac_project$design_matrix_path <- file.path(atac_browser_root, "manifest", "missing_design_matrix.txt")
external_atac_project$fastq_dir <- ""
external_atac_project$fastq_dirs <- character(0)
external_atac_sheet <- app_env$genome_browser_infer_atac_comparison_sheet(external_atac_project, akps_comparison_dir)
assert(NROW(external_atac_sheet) == 4L && all(startsWith(external_atac_sheet$sample, "AKPS_")), "completed ATAC import infers the exact AKPS samples from its comparison directory")
assert(identical(external_atac_sheet$condition, c("AKPS_AA", "AKPS_AA", "AKPS_Veh", "AKPS_Veh")), "completed ATAC import keeps AA signal tracks above Veh tracks")
assert(identical(app_env$validate_completed_atac_results(atac_browser_root), normalizePath(atac_browser_root, winslash = "/", mustWork = TRUE)), "completed ATAC import validates BigWigs and DiffBind BED output before registration")

annotation_inputs <- app_env$peak_annotation_input_files(atac_project)
assert(legacy_peak %in% annotation_inputs, "peak annotation discovers completed per-sample MACS2 peaks")
annotation_root <- file.path(root, "peak_annotation")
dir.create(annotation_root, recursive = TRUE)
annotation_jobs <- data.frame(
  step = c("Peak Annotation", "Peak Annotation"),
  slurm_state = c("RUNNING", "FAILED"), stringsAsFactors = FALSE
)
assert(identical(app_env$peak_annotation_status(atac_project, annotation_jobs), "Active"), "an active annotation job is not hidden by a newer stale job record")
writeLines("status\trunning", file.path(annotation_root, "_RUN_STARTED"))
assert(identical(app_env$peak_annotation_status(atac_project, data.frame()), "Likely failed"), "orphaned annotation run marker reports an incomplete job")
unlink(file.path(annotation_root, "_RUN_STARTED"))
annotated_peak <- file.path(sample_dir, "A1_peaks_annotated.txt")
writeLines("PeakID\tAnnotation\npeak1|Fold=2.1|p.value=0.0004|FDR=0.008\tPromoter", annotated_peak)
expanded_annotation <- app_env$safe_read_result_table(annotated_peak)
assert(all(c("Fold", "p.value", "FDR") %in% names(expanded_annotation)), "embedded peak statistics expand into explicit result columns")
assert(is.numeric(expanded_annotation$p.value) && identical(expanded_annotation$p.value[[1]], 0.0004), "expanded ATAC differential p-value remains numeric")
write.table(data.frame(
  result_type = "MACS2", sample_or_comparison = "A1", peak_count = 2,
  source_peak_file = legacy_peak, annotated_file = annotated_peak, status = "complete",
  stringsAsFactors = FALSE
), file.path(annotation_root, "peak_annotation_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
Sys.sleep(0.01)
writeLines(c("status\tcomplete", "annotated_files\t1"), file.path(annotation_root, "_COMPLETE"))
assert(app_env$peak_annotation_is_current(atac_project, legacy_peak), "current annotations are not resubmitted")
assert(identical(app_env$peak_annotation_status(atac_project, data.frame()), "Complete"), "annotation completion marker drives pipeline status")
assert(NROW(app_env$peak_annotation_summary_table(atac_project)) == 1L, "annotation summary is available to the Results Explorer")
assert(annotated_peak %in% unname(app_env$peak_annotation_result_files(atac_project)), "annotated peak tables are discoverable in Results Explorer")

counts_upload <- file.path(root, "uploaded_counts.csv")
write.csv(data.frame(
  Gene = c("GeneA", "GeneA", "GeneB"),
  Chr = c("chr1", "chr1", "chr2"),
  Length = c(100, 100, 200),
  S1 = c(1.5, 2.4, 0),
  S2 = c(3.5, 1.4, 5),
  S3 = c(2, 2, 7),
  check.names = FALSE
), counts_upload, row.names = FALSE, quote = FALSE)
counts_only_root <- file.path(root, "counts-only", "data")
count_result <- app_env$standardize_uploaded_count_matrix(
  counts_upload,
  file.path(counts_only_root, "counts", "count_matrix.txt")
)
standardized_counts <- read.delim(count_result$path, check.names = FALSE)
assert(identical(count_result$samples, c("S1", "S2", "S3")), "counts upload preserves numeric sample column names")
assert(NROW(standardized_counts) == 2L && standardized_counts$S1[standardized_counts$Geneid == "GeneA"] == 4, "counts upload rounds half-up and sums duplicate genes")
counts_design <- app_env$write_counts_only_design(
  count_result$samples,
  file.path(counts_only_root, "manifest", "design_matrix.txt"),
  metadata_cols = "treatment, batch"
)
counts_design_df <- read.delim(counts_design, check.names = FALSE)
counts_design_df$treatment <- c("Control", "Control", "Treated")
counts_design_df$batch <- c("B1", "B2", "B1")
write.table(counts_design_df, counts_design, sep = "\t", row.names = FALSE, quote = FALSE)
counts_project <- list(
  id = "rna/counts-only", name = "counts-only", analysis_key = "rna", analysis = "RNA-seq",
  counts_only = TRUE, design_matrix_path = counts_design, data_dir = counts_only_root,
  results_root = root, fastq_dir = "", fastq_dirs = character(0), paired_end = TRUE, genome = "mouse"
)
assert(identical(app_env$pipeline_order(counts_project), c("Design matrix", "DESeq2", "GSEA")), "counts-only projects expose only design, DESeq2, and GSEA")
modeled_design <- app_env$deseq_design_for_column(counts_project, "treatment", "batch")
modeled_df <- read.delim(modeled_design, check.names = FALSE)
assert(identical(names(modeled_df), c("sample", "batch", "treatment", "filename")), "selected covariates are retained and the comparison variable is modeled last")
dir.create(file.path(counts_only_root, "deseq2"), recursive = TRUE, showWarnings = FALSE)
writeLines(
  "DESCRIPTION\tS1\tS2\tS3",
  file.path(counts_only_root, "deseq2", "normalized_counts_Treated_vs_Control(ref).txt")
)
writeLines(
  "gene\tbaseMean\tlog2FoldChange\tpvalue\tpadj",
  file.path(counts_only_root, "deseq2", "DEG_Treated_vs_Control(ref).txt")
)
completed_deseq <- app_env$completed_deseq_comparison_catalog(counts_project)
assert(
  NROW(completed_deseq) == 1L &&
    identical(completed_deseq$compare_col[[1]], "treatment") &&
    identical(completed_deseq$comparison[[1]], "Treated") &&
    identical(completed_deseq$reference[[1]], "Control"),
  "GSEA discovers completed DESeq2 comparisons and maps them to the correct design column"
)
assert(
  grepl("Comparisons (select one or more)", app_text, fixed = TRUE) &&
    grepl("vapply(comparisons, function(comparison)", app_text, fixed = TRUE),
  "DESeq2 supports submitting multiple selected comparisons as independent jobs"
)
assert(
  grepl("Completed DESeq2 comparison", app_text, fixed = TRUE) &&
    grepl("completed_deseq_comparison_catalog(current_project())", app_text, fixed = TRUE),
  "GSEA has an independent dropdown populated from completed DESeq2 runs"
)
assert(
  app_env$PROGRESS_REFRESH_MS >= 20000 &&
    app_env$SLURM_QUERY_TIMEOUT_SECONDS <= 2 &&
    app_env$MAX_SLURM_JOB_IDS_PER_REFRESH <= 100,
  "progress polling is bounded so scheduler latency cannot repeatedly freeze the Shiny event loop"
)
assert(
  isTRUE(app_env$recent_submission(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), now = Sys.time())) &&
    !isTRUE(app_env$recent_submission("2020-01-01 00:00:00", now = Sys.time())),
  "only newly submitted unresolved jobs retain the temporary active status"
)
assert(
  grepl("project_status_state(data.frame())", server_source, fixed = TRUE) &&
    grepl("submission_holds(list())", server_source, fixed = TRUE),
  "switching projects clears status colors and submission handoff state before refreshing"
)
assert(
  grepl("previous_cache$value", app_text, fixed = TRUE) &&
    grepl("sacct_ids <- setdiff(ids, ids_seen)", app_text, fixed = TRUE),
  "terminal SLURM states are cached and active queue jobs are not redundantly queried through sacct"
)
assert(
  grepl("completed_deseq_result_files(current_project())", app_text, fixed = TRUE) &&
    grepl("completed_gsea_result_files(current_project())", app_text, fixed = TRUE),
  "Results Explorer selectors only list files from completed DESeq2 and GSEA runs"
)

scrna_source <- file.path(root, "single_input.rds")
file.create(scrna_source)
scrna_direct <- app_env$scrna_manifest_from_setup("", "donor_01", scrna_source)
assert(
  NROW(scrna_direct) == 1L && identical(scrna_direct$sample_id[[1]], "donor_01") && identical(scrna_direct$input_path[[1]], scrna_source),
  "scRNA projects can start directly from a single server input without a manifest"
)
scrna_manifest_path <- file.path(root, "provided_scRNA_manifest.tsv")
write.table(data.frame(sample_id = "donor_02", input_path = scrna_source), scrna_manifest_path, sep = "\t", row.names = FALSE, quote = FALSE)
assert(
  identical(app_env$scrna_manifest_from_setup(scrna_manifest_path)$sample_id[[1]], "donor_02"),
  "scRNA projects still accept an optional supplied manifest"
)
scrna_project <- list(
  id = "scrna/test", name = "test", analysis_key = "scrna", analysis = "scRNA-seq",
  data_dir = file.path(root, "scRNA-results", "data"), results_root = file.path(root, "scRNA-results"),
  design_matrix_path = scrna_manifest_path, scrna_input_manifest = scrna_manifest_path, scrna_engine = "auto"
)
failed_scrna_job <- data.frame(step = "Input inspection", slurm_state = "FAILED", stringsAsFactors = FALSE)
assert(
  identical(app_env$project_status(scrna_project, jobs = failed_scrna_job)$status[[1]], "Likely failed"),
  "scRNA terminal failures are surfaced as failed rather than not started"
)
processed_rds <- file.path(root, "existing_processed.rds")
writeLines("fixture", processed_rds)
processed_manifest <- file.path(root, "existing_processed_manifest.tsv")
write.table(data.frame(sample_id = "existing", input_path = processed_rds), processed_manifest, sep = "\t", row.names = FALSE, quote = FALSE)
processed_project <- list(
  id = "scrna/processed", name = "processed", analysis_key = "scrna", analysis = "scRNA-seq",
  data_dir = file.path(root, "processed-results", "data"), results_root = file.path(root, "processed-results"),
  design_matrix_path = processed_manifest, scrna_input_manifest = processed_manifest, scrna_engine = "seurat", scrna_input_mode = "single"
)
dir.create(file.path(processed_project$data_dir, "scrna", "tables"), recursive = TRUE)
write.table(data.frame(sample_id = "existing", input_kind = "seurat_rds", pca_detected = TRUE, umap_detected = TRUE, clusters_detected = TRUE, annotation_columns_detected = "cell_type"), file.path(processed_project$data_dir, "scrna", "tables", "input_processing_detected.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
dir.create(file.path(processed_project$data_dir, "scrna", "objects"), recursive = TRUE)
writeLines("project-local processed object", file.path(processed_project$data_dir, "scrna", "objects", "processed_seurat.rds"))
assert(
  app_env$scrna_source_is_processed_object(processed_project) && app_env$scrna_existing_processed_ready(processed_project) &&
    all(app_env$project_status(processed_project)$status[app_env$project_status(processed_project)$step %in% c("Normalize & PCA", "UMAP & clustering", "Annotate & markers")] == "Complete"),
  "a processed object with PCA, UMAP, clusters, and annotations is adopted for continuation without mandatory reconstruction"
)
assert(
  grepl("scanpy_container_check", app_text, fixed = TRUE) && grepl("scrna_manifest_from_setup", app_text, fixed = TRUE),
  "scRNA setup supports an optional manifest and a shared Scanpy container"
)

fastq_dir <- file.path(root, "donor_fastqs")
dir.create(fastq_dir)
for (read in c("R1", "R2")) {
  con <- gzfile(file.path(fastq_dir, paste0("donor03_S1_L001_", read, "_001.fastq.gz")), "wt")
  writeLines(c("@read", "ACGT", "+", "####"), con)
  close(con)
}
fastq_manifest_path <- file.path(root, "fastq_scRNA_design.tsv")
write.table(data.frame(sample_id = "donor03", input_path = fastq_dir, input_type = "fastq_folder", fastq_sample = "donor03"), fastq_manifest_path, sep = "\t", row.names = FALSE, quote = FALSE)
fastq_project <- list(
  id = "scrna/fastq", name = "fastq", analysis_key = "scrna", analysis = "scRNA-seq",
  data_dir = file.path(root, "fastq-scRNA-results", "data"), results_root = file.path(root, "fastq-scRNA-results"),
  design_matrix_path = fastq_manifest_path, scrna_input_manifest = fastq_manifest_path, scrna_engine = "seurat", scrna_input_mode = "multiple"
)
validated_fastq <- app_env$validate_scrna_manifest(app_env$scrna_manifest(fastq_project))
assert(identical(validated_fastq$input_type[[1]], "fastq_folder"), "scRNA sample designs detect matched 10x FASTQ folders")
assert(identical(app_env$scrna_fastq_sample_names(fastq_dir), "donor03"), "10x FASTQ filename prefixes are inferred for Cell Ranger sample selection")
alignment_not_started <- app_env$scrna_alignment_sample_progress(fastq_project, jobs = data.frame())$table
alignment_waiting <- app_env$scrna_alignment_sample_progress(fastq_project, jobs = data.frame(
  step = "Alignment & counting", sample = "donor03", slurm_state = "Submitted", elapsed = "00:00:00", stringsAsFactors = FALSE
))$table
alignment_running <- app_env$scrna_alignment_sample_progress(fastq_project, jobs = data.frame(
  step = "Alignment & counting", sample = "donor03", slurm_state = "RUNNING", elapsed = "00:04:12", stringsAsFactors = FALSE
))$table
alignment_failed <- app_env$scrna_alignment_sample_progress(fastq_project, jobs = data.frame(
  step = "Alignment & counting", sample = "donor03", slurm_state = "FAILED", elapsed = "00:01:05", stringsAsFactors = FALSE
))$table
assert(
  identical(alignment_not_started$status[[1]], "Not started") &&
    identical(alignment_waiting$status[[1]], "Waiting") &&
    identical(alignment_running$status[[1]], "Running") && identical(alignment_running$time_running[[1]], "00:04:12") &&
    identical(alignment_failed$status[[1]], "Likely failed"),
  "scRNA alignment progress distinguishes not-started, queued, running, and failed Cell Ranger samples"
)
alignment_optimistic <- app_env$optimistic_step_progress(fastq_project, "Alignment & counting", samples = "donor03")
assert(
  NROW(alignment_optimistic) == 1L && identical(alignment_optimistic$status[[1]], "Waiting") &&
    grepl("cellranger/donor03$", alignment_optimistic$target[[1]]),
  "new Cell Ranger submissions immediately render as waiting for the selected sample"
)
pooled_fastq_dir <- file.path(root, "pooled_fastqs")
dir.create(pooled_fastq_dir)
pooled_prefixes <- c("Amor_AH07_Y_UT_F", "Amor_AH07_Y_UT_M", "Amor_AH07_O_uPAR_F", "Amor_AH07_O_uPAR_M")
for (prefix in pooled_prefixes) for (lane in c("L001", "L002")) for (read in c("R1", "R2")) {
  con <- gzfile(file.path(pooled_fastq_dir, paste0(prefix, "_S1_", lane, "_", read, "_001.fastq.gz")), "wt")
  writeLines(c("@read", "ACGT", "+", "####"), con)
  close(con)
}
pooled_rows <- app_env$scrna_fastq_input_rows(pooled_fastq_dir)
assert(
  NROW(pooled_rows) == length(pooled_prefixes) &&
    identical(sort(pooled_rows$fastq_sample), sort(pooled_prefixes)) &&
    length(unique(pooled_rows$input_path)) == 1L,
  "one pooled FASTQ folder expands into one Cell Ranger design row per sample prefix while lanes remain grouped"
)
matrix_pool <- file.path(root, "matrix_pool")
for (sample in c("sample_a", "sample_b")) {
  sample_matrix <- file.path(matrix_pool, sample, "outs", "filtered_feature_bc_matrix")
  dir.create(sample_matrix, recursive = TRUE)
  for (name in c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")) {
    con <- gzfile(file.path(sample_matrix, name), "wt"); writeLines("fixture", con); close(con)
  }
}
matrix_rows <- app_env$scrna_matrix_input_rows(matrix_pool)
assert(
  NROW(matrix_rows) == 2L && identical(sort(matrix_rows$sample_id), c("sample_a", "sample_b")),
  "one parent folder discovers multiple filtered feature-barcode matrix subfolders"
)
browser_fixture <- file.path(root, "absolute_path_fixture", "nested")
dir.create(browser_fixture, recursive = TRUE)
browser_file_fixture <- file.path(browser_fixture, "samples.tsv")
writeLines("sample_id\tinput_path", browser_file_fixture)
assert(identical(normalizePath(app_env$browser_start_path(browser_fixture, "dir"), winslash = "/"), normalizePath(browser_fixture, winslash = "/")), "a pasted absolute folder opens at that exact folder")
assert(identical(normalizePath(app_env$browser_start_path(browser_file_fixture, "file"), winslash = "/"), normalizePath(browser_fixture, winslash = "/")), "a pasted absolute file opens at its containing folder")
assert(
  grepl('tabPanel("Samples & Design"', app_text, fixed = TRUE) &&
    grepl('updateTabsetPanel(session, "web_main_tabs", selected = "Samples & Design")', server_source, fixed = TRUE) &&
    grepl('output$new_scrna_inputs_editor_ui <- renderUI', server_source, fixed = TRUE) &&
    grepl('design_form_table_ui(samples, prefix = "new_scrna_form"', server_source, fixed = TRUE),
  "detected scRNA samples open in the top-level RNA-style form editor"
)
assert(
  grepl('typed_value <- path.expand', server_source, fixed = TRUE) &&
    grepl('value <- if (nzchar(typed_value))', server_source, fixed = TRUE) &&
    grepl('path_browser$path', server_source, fixed = TRUE),
  "the server browser uses a pasted absolute path instead of silently falling back to its previous folder"
)
object_project_config <- app_env$new_project_from_inputs(list(
  new_project_analysis = "scRNA-seq", new_project_name = "object_fixture", new_project_mode = "new",
  new_results_root = root, new_scrna_start_mode = "object", new_scrna_folder_type = "filtered_10x_matrix"
))
fastq_project_config <- app_env$new_project_from_inputs(list(
  new_project_analysis = "scRNA-seq", new_project_name = "fastq_fixture", new_project_mode = "new",
  new_results_root = root, new_scrna_start_mode = "new", new_scrna_folder_type = "fastq_folder",
  new_species = "human", new_genome_version = "refdata-gex-GRCh38-2024-A"
))
assert(identical(object_project_config$genome, "auto") && !nzchar(object_project_config$genome_version), "processed scRNA objects do not inherit hidden genome/reference selections")
assert(identical(fastq_project_config$genome, "human") && identical(fastq_project_config$genome_version, "refdata-gex-GRCh38-2024-A"), "FASTQ scRNA projects retain the selected Cell Ranger transcriptome version")
assert("Alignment & counting" %in% app_env$scrna_pipeline_order(fastq_project), "FASTQ-backed scRNA projects expose alignment and counting before input inspection")
assert(
  identical(app_env$canonical_job_step("Cell Ranger count"), "Alignment & counting") && identical(app_env$canonical_job_step("Alignment & counting"), "Alignment & counting"),
  "legacy Cell Ranger job records remain attached to the renamed alignment and counting stage"
)
assert(
  grepl("run_scrna_cellranger", app_text, fixed = TRUE) && any(grepl("module load CellRanger/9.0.1", runtime_text, fixed = TRUE)),
  "the app submits FASTQ samples through the maintained CellRanger/9.0.1 cluster module"
)
assert(
  grepl("Reference inspection is running", app_text, fixed = TRUE) &&
    grepl("Inspecting reference labels", app_text, fixed = TRUE) &&
    grepl("scrna_reference_inspect_button_ui", app_text, fixed = TRUE),
  "reference-label inspection has a persistent running state and disabled progress button"
)
assert(
  grepl("scrna_annotation_ui_settings", app_text, fixed = TRUE) &&
    grepl("Uploaded reference retained for this project", app_text, fixed = TRUE) &&
    grepl("reference-label-preview", app_text, fixed = TRUE),
  "reference annotation retains the selected method/upload and shows readable label previews"
)
assert(
  grepl('"reference_ortholog_file"', app_text, fixed = TRUE) &&
    grepl("Reference and query species are detected", app_text, fixed = TRUE) &&
    grepl("bundled MGI ortholog table", app_text, fixed = TRUE),
  "reference annotation submits the bundled ortholog table and explains automatic mouse-human conversion"
)
assert(
  grepl('checkboxInput("scrna_find_cluster_markers"', app_text, fixed = TRUE) &&
    grepl('find_cluster_markers = isTRUE(input$scrna_find_cluster_markers)', app_text, fixed = TRUE) &&
    grepl('"find_cluster_markers"', app_text, fixed = TRUE) &&
    grepl("This metadata-only step adds the transferred label", app_text, fixed = TRUE),
  "cluster-marker discovery is an explicit optional annotation setting rather than an automatic large-cell calculation"
)
assert(
  grepl('button_ui = uiOutput("scrna_annotation_run_button_ui")', app_text, fixed = TRUE) &&
    grepl('output$scrna_annotation_run_button_ui <- renderUI', app_text, fixed = TRUE) &&
    grepl('if (!ready || is.null(input$scrna_reference_label_column)) return(NULL)', app_text, fixed = TRUE) &&
    grepl('"Run reference annotation"', app_text, fixed = TRUE),
  "the annotation submit button renders only after its dynamic controls and inspected reference labels are ready"
)
assert(
  grepl("cslReportOpenToolPanels", app_text, fixed = TRUE) &&
    grepl("preserve_reference_annotation", app_text, fixed = TRUE) &&
    grepl('"tool_panel_annotate_markers" %in% isolate(open_tool_panels())', app_text, fixed = TRUE) &&
    grepl("its dynamic upload controls stay", app_text, fixed = TRUE),
  "background status refreshes preserve an open reference-upload annotation panel"
)
assert(
  grepl('scrna_stage_resource_options <- function', app_text, fixed = TRUE) &&
    identical(app_env$scrna_stage_resource_options("inspect", 0, "seurat"), c("--cpus-per-task=6", "--mem=48G", "--time=1-00:00:00")) &&
    identical(app_env$scrna_stage_resource_options("cluster", 3 * 1024^3, "seurat"), c("--cpus-per-task=16", "--mem=160G", "--time=2-00:00:00")) &&
    identical(app_env$scrna_stage_resource_options("cluster", 3 * 1024^3, "scanpy"), c("--cpus-per-task=20", "--mem=160G", "--time=2-00:00:00")) &&
    identical(app_env$scrna_stage_resource_options("annotate", 10 * 1024^3, "seurat"), c("--cpus-per-task=20", "--mem=192G", "--time=3-00:00:00")) &&
    identical(app_env$scrna_stage_resource_options("annotate", 30 * 1024^3, "seurat"), c("--cpus-per-task=24", "--mem=256G", "--time=3-00:00:00")) &&
    grepl('scrna_stage_input_bytes(stage, resolved_engine, manifest, out_dir, reference_file)', app_text, fixed = TRUE) &&
    grepl('scrna_cellranger_resource_options(c(pairs$r1, pairs$r2))', app_text, fixed = TRUE) &&
    grepl('sbatch_options <- scrna_stage_resource_options', app_text, fixed = TRUE),
  "single-cell stages scale SLURM resources from the relevant input/checkpoint size and reserve the highest tier for genuinely large objects"
)
idle_reference_html <- as.character(app_env$scrna_reference_label_selector_content(
  list(status = "idle", project_id = "scrna/test", choices = data.frame()),
  "scrna/test"
))
submitted_reference_html <- as.character(app_env$scrna_reference_label_selector_content(
  list(status = "submitted", project_id = "scrna/test", job_id = "123", choices = data.frame()),
  "scrna/test"
))
ready_reference_html <- as.character(app_env$scrna_reference_label_selector_content(
  list(
    status = "ready", project_id = "scrna/test", message = "Ready",
    choices = data.frame(value = c("", "cell_type"), source = c("Active identities", "cell_type"), label_count = c(32L, 22L), non_missing_cells = c(7497L, 7497L), preview = c("A | B", "C | D"), stringsAsFactors = FALSE)
  ),
  "scrna/test"
))
other_project_reference_html <- as.character(app_env$scrna_reference_label_selector_content(
  list(status = "ready", project_id = "scrna/other", choices = data.frame(value = "", source = "Active identities", label_count = 32L, non_missing_cells = 7497L)),
  "scrna/test"
))
assert(
  !grepl("Active identities", idle_reference_html, fixed = TRUE) &&
    !grepl("scrna_reference_label_column", idle_reference_html, fixed = TRUE) &&
    !grepl("scrna_reference_label_column", submitted_reference_html, fixed = TRUE) &&
    grepl("Active identities (32 labels)", ready_reference_html, fixed = TRUE) &&
    grepl("cell_type (22 labels)", ready_reference_html, fixed = TRUE) &&
    !grepl("Active identities", other_project_reference_html, fixed = TRUE),
  "reference labels render only after a successful inspection of the current project's actual reference"
)
assert(
  grepl("scrna_annotation_method_intent", app_text, fixed = TRUE) &&
    grepl("click.cslAnnotationIntent_", app_text, fixed = TRUE) &&
    grepl('annotation_method <- annotation_ui$method %||% input$scrna_annotation_method', app_text, fixed = TRUE) &&
    !grepl('tags$input(type = "hidden", id = "scrna_annotation_method"', app_text, fixed = TRUE) &&
    !grepl('tags$input(type = "hidden", id = "scrna_marker_source"', app_text, fixed = TRUE) &&
    grepl('selected_method <- if (pbmc_example) "markers"', app_text, fixed = TRUE) &&
    grepl('selected_marker_source <- if (pbmc_example) "server"', app_text, fixed = TRUE) &&
    grepl('scrna_marker_assignment: \'scrna_marker_assignment_intent\'', app_text, fixed = TRUE) &&
    grepl('marker_assignment = "automatic"', app_text, fixed = TRUE) &&
    grepl('controls_ready <- all(vapply(suggestions$cluster', app_text, fixed = TRUE) &&
    grepl('default_signature_source <- if (pbmc_example', app_text, fixed = TRUE),
  "annotation method and source persistence is driven by explicit user intent rather than dynamic-control recreation"
)
assert(
  is.infinite(app_env$scrna_cellranger_max_parallel()) &&
    grepl("dependency_condition = \"afterany\"", app_text, fixed = TRUE) &&
    grepl("BAM generation is disabled", app_text, fixed = TRUE) &&
    any(grepl("--create-bam=false", runtime_text, fixed = TRUE)),
  "Cell Ranger defaults to concurrent resumable low-disk jobs with BAM output disabled"
)
matrix_dir <- app_env$scrna_cellranger_matrix_dir(fastq_project, "donor03")
dir.create(matrix_dir, recursive = TRUE)
for (name in c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")) {
  con <- gzfile(file.path(matrix_dir, name), "wt"); writeLines("fixture", con); close(con)
}
writeLines("complete", file.path(app_env$scrna_cellranger_output_dir(fastq_project, "donor03"), "_CELLRANGER_COMPLETE"))
metrics_path <- file.path(app_env$scrna_cellranger_output_dir(fastq_project, "donor03"), "outs", "metrics_summary.csv")
utils::write.csv(data.frame(`Estimated Number of Cells` = "1,234", `Mean Reads per Cell` = "45,678", check.names = FALSE), metrics_path, row.names = FALSE, quote = TRUE)
alignment_complete <- app_env$scrna_alignment_sample_progress(fastq_project, jobs = data.frame(
  step = "Alignment & counting", sample = "donor03", slurm_state = "COMPLETED", elapsed = "00:31:00", stringsAsFactors = FALSE
))$table
alignment_display <- app_env$sample_progress_step_table(alignment_complete, "Alignment & counting")
alignment_ui <- as.character(app_env$sample_progress_step_ui(alignment_complete, "Alignment & counting"))
assert(
  identical(alignment_complete$status[[1]], "Completed") &&
    all(c("Sample", "Status", "Time running", "Estimated cells", "Mean reads/cell") %in% names(alignment_display)) &&
    identical(alignment_display$Sample[[1]], "donor03") &&
    grepl("tool-progress-table", alignment_ui, fixed = TRUE) && grepl("sample-status completed", alignment_ui, fixed = TRUE),
  "completed Cell Ranger samples use the RNA-seq per-sample status table with elapsed time and available metrics"
)
resolved_fastq <- app_env$scrna_resolved_manifest(fastq_project)
assert(
  identical(resolved_fastq$data$input_type[[1]], "filtered_10x_matrix") && identical(resolved_fastq$data$input_path[[1]], matrix_dir),
  "completed Cell Ranger matrices replace FASTQ folders only in the downstream execution manifest"
)

object_dir <- file.path(app_env$scrna_output_dir(scrna_project), "objects")
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
assert(is.null(app_env$scrna_processed_object_info(scrna_project)), "the object download stays unavailable before a final object exists")
scanpy_object <- file.path(object_dir, "processed_scanpy.h5ad")
seurat_object <- file.path(object_dir, "processed_seurat.rds")
writeBin(as.raw(seq_len(16)), scanpy_object)
writeBin(as.raw(seq_len(8)), seurat_object)
scrna_project$scrna_engine <- "scanpy"
download_info <- app_env$scrna_processed_object_info(scrna_project)
assert(
  identical(download_info$engine, "scanpy") && identical(download_info$filename, "processed_scanpy.h5ad") && identical(download_info$size, "16 B"),
  "the scRNA object download selects the project engine's non-empty final object"
)
assert(
  grepl('uiOutput("scrna_processed_object_download_ui")', app_text, fixed = TRUE) &&
    grepl("output$download_scrna_processed_object <- downloadHandler", server_source, fixed = TRUE),
  "the three-tab scRNA explorer exposes the processed-object download section and handler"
)

cat("CodeSpringApp fake-data helper smoke tests passed.\n")
