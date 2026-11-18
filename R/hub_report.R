# CodeSpringHub post-run report: one self-contained HTML file per run.
#
# The packet only gathers what CodeSpringApp's own steps already wrote (the DESeq2 step's plots and
# tables, GSEApy's reports, the tools' QC numbers, the app's methods text). It computes no new
# statistics. Plots and small data files are embedded so the page works offline, can be emailed or
# downloaded from a chat, and every plot/file has a download link. Anything that would make the
# page too big (count matrices, BAMs, Kallisto/RSEM matrices) is listed with its path on the
# cluster instead, and the page says how to open the full results in CodeSpringApp.
#
# Sourced after hub_bridge.R by bin/codespring-run (`codespring-run report PROJECT_ID`), and called
# by the run controller when a run finishes.

HUB_REPORT_MAX_MB <- local({ n <- suppressWarnings(as.numeric(Sys.getenv("CSL_HUB_REPORT_MAX_MB", unset = "15"))); if (!is.finite(n) || n < 1) 15 else n })
HUB_REPORT_FILE_MB <- local({ n <- suppressWarnings(as.numeric(Sys.getenv("CSL_HUB_REPORT_FILE_MB", unset = "5"))); if (!is.finite(n) || n <= 0) 5 else n })
HUB_REPORT_MAX_PLOTS <- 30L

hub_report_path <- function(p) file.path(project_result_dir(p), "report", paste0(p$name, "_codespring_report.html"))

hub_html <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}

hub_mime <- function(path) {
  switch(tolower(tools::file_ext(path)), png = "image/png", csv = "text/csv", txt = , tsv = "text/tab-separated-values",
         "application/octet-stream")
}

hub_data_uri <- function(path) paste0("data:", hub_mime(path), ";base64,", gsub("[\r\n]", "", jsonlite::base64_enc(readBin(path, "raw", file.size(path)))))

hub_html_table <- function(df, max_rows = 25L, digits = 3L) {
  if (!NROW(df)) return("<p class=\"muted\">Nothing to show.</p>")
  df <- utils::head(df, max_rows)
  fmt <- function(v) {
    if (!is.numeric(v) || is.na(v)) return(if (is.na(v)) "" else as.character(v))
    if (v != 0 && abs(v) < 1e-3) return(formatC(v, format = "e", digits = 2))
    if (v == round(v)) return(formatC(v, format = "d", big.mark = ","))
    formatC(round(v, digits), format = "f", digits = digits, big.mark = ",", drop0trailing = TRUE)
  }
  head_row <- paste0("<tr>", paste0("<th>", hub_html(names(df)), "</th>", collapse = ""), "</tr>")
  body <- vapply(seq_len(NROW(df)), function(i) paste0("<tr>", paste0("<td>", hub_html(vapply(df, function(col) fmt(col[i]), "")), "</td>", collapse = ""), "</tr>"), "")
  paste0("<div class=\"scroll\"><table>", head_row, paste(body, collapse = ""), "</table></div>")
}

hub_read_table <- function(path) {
  # DESeq2 tables have one fewer header field than columns (gene names as row names).
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  df <- tryCatch(utils::read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "", check.names = FALSE,
                                   stringsAsFactors = FALSE, fill = TRUE), error = function(e) NULL)
  if (is.null(df)) return(NULL)
  if (!is.null(rownames(df)) && !identical(rownames(df), as.character(seq_len(NROW(df))))) df <- cbind(gene = rownames(df), df, stringsAsFactors = FALSE)
  rownames(df) <- NULL
  df
}

hub_report_files <- function(p) {
  # Everything the report could show, most useful first. kind: plot / data / big (path only).
  data_dir <- p$data_dir
  deseq <- file.path(data_dir, "deseq2")
  out <- list()
  add <- function(path, kind, section, label, comparison = NA_character_) {
    if (length(path) && file.exists(path) && file.size(path) > 0) out[[length(out) + 1]] <<- list(path = path, kind = kind, section = section, label = label, comparison = comparison)
  }
  degs <- sort(list.files(deseq, pattern = "^DEG_.*\\.txt$", full.names = TRUE))
  full <- degs[!grepl("^DEG_fdr", basename(degs))]
  comparisons <- sub("^DEG_(.*)\\.txt$", "\\1", basename(full))
  pngs <- sort(list.files(deseq, pattern = "\\.png$", full.names = TRUE))
  used <- character(0)
  for (cmp in comparisons) {
    mine <- pngs[endsWith(basename(pngs), paste0(cmp, ".png"))]
    for (f in mine) add(f, "plot", "deseq2", basename(f), cmp)
    used <- c(used, mine)
  }
  for (f in setdiff(pngs, used)) add(f, "plot", "deseq2", basename(f))
  gsea_dir <- file.path(data_dir, "gseapy")
  gsea_reports <- sort(list.files(gsea_dir, pattern = "^report\\.gseapy\\..*\\.csv$", recursive = TRUE, full.names = TRUE))
  for (f in sort(list.files(gsea_dir, pattern = "^DotPlot.*\\.png$", recursive = TRUE, full.names = TRUE))) {
    add(f, "plot", "gsea", paste(basename(dirname(dirname(f))), basename(dirname(f)), sep = " / "))
  }
  for (cmp in comparisons) add(file.path(deseq, paste0("DEG_fdr0.05_", cmp, ".txt")), "data", "deseq2", paste0("DEG_fdr0.05_", cmp, ".txt"), cmp)
  add(file.path(deseq, "design_matrix.txt"), "data", "deseq2", "design_matrix.txt")
  for (f in gsea_reports) add(f, "data", "gsea", paste0(basename(dirname(dirname(f))), "_", basename(f)))
  for (cmp in comparisons) add(file.path(deseq, paste0("DEG_", cmp, ".txt")), "data", "deseq2", paste0("DEG_", cmp, ".txt"), cmp)
  for (f in sort(list.files(deseq, pattern = "^normalized_counts_.*\\.txt$", full.names = TRUE))) add(f, "data", "deseq2", basename(f))
  counts <- file.path(data_dir, "counts")
  add(file.path(counts, "featurecounts_summary.txt"), "data", "counts", "featurecounts_summary.txt")
  add(file.path(counts, "count_matrix.txt"), "data", "counts", "count_matrix.txt")
  for (f in sort(list.files(counts, pattern = "^(kallisto|rsem)_.*\\.txt$", full.names = TRUE))) add(f, "data", "counts", basename(f))
  # ATAC-seq / ChIP-seq: DiffBind's plots and tables per comparison, then the HOMER annotation summary.
  diffbind <- sort(list.dirs(file.path(data_dir, "diffbind"), recursive = FALSE, full.names = TRUE))
  for (d in diffbind) {
    for (f in sort(list.files(d, pattern = "\\.png$", full.names = TRUE))) add(f, "plot", "diffbind", basename(f), basename(d))
    for (f in sort(list.files(d, pattern = "^DifferentialPeaks.*\\.txt$", full.names = TRUE))) add(f, "data", "diffbind", paste0(basename(d), "_", basename(f)), basename(d))
    add(file.path(d, "diffbind_sample_sheet.tsv"), "data", "diffbind", paste0(basename(d), "_diffbind_sample_sheet.tsv"), basename(d))
  }
  add(file.path(data_dir, "peak_annotation", "peak_annotation_summary.tsv"), "data", "annotation", "peak_annotation_summary.tsv")
  list(files = out, comparisons = comparisons, gsea_reports = gsea_reports, diffbind = basename(diffbind))
}

hub_report_qc_table <- function(qc_rows) {
  num <- function(x) suppressWarnings(as.numeric(gsub("[,%]", "", as.character(x %||% NA))))
  rows <- lapply(qc_rows, function(r) {
    data.frame(
      sample = r$sample,
      `reads in` = num(r$cutadapt[["Reads before"]] %||% r$star$input_reads %||% NA),
      `kept after trimming` = num(r$cutadapt[["Reads after"]] %||% NA),
      `FastQC pass/warn/fail` = if (!is.null(r$fastqc)) paste(r$fastqc$pass, r$fastqc$warn, r$fastqc$fail, sep = "/") else "",
      `FastQC FAIL modules` = paste(unlist(r$fastqc$failed_modules), collapse = "; "),
      `STAR unique %` = num(r$star$uniquely_mapped_pct %||% NA),
      `STAR multi %` = num(r$star$multi_mapped_pct %||% NA),
      `featureCounts assigned %` = num(r$featurecounts[["Assigned %"]] %||% NA),
      `Kallisto pseudoaligned %` = num(r$kallisto$pseudoaligned_pct %||% NA),
      role = r$role %||% "",
      `Bowtie2 mapped reads` = num(r$bowtie2$mapped_reads %||% NA),
      `kept after dedup %` = num(r$bowtie2$kept_after_dedup_pct %||% NA),
      `MACS2 peaks` = num(r$macs2$peaks %||% NA),
      check.names = FALSE, stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  if (!NROW(df)) return(df)
  keep <- vapply(df, function(col) !all(is.na(col) | col == ""), logical(1))
  df[, keep, drop = FALSE]
}

hub_report <- function(project_id) {
  p <- hub_project(project_id)
  st <- tryCatch(hub_status(p$id), error = function(e) list(state = "unknown", steps = list()))
  qc <- tryCatch(hub_qc(p$id), error = function(e) list(samples = list(), summary = list()))
  found <- hub_report_files(p)
  budget <- HUB_REPORT_MAX_MB * 1024^2 - 300 * 1024  # leave room for the page itself
  per_file <- HUB_REPORT_FILE_MB * 1024^2
  embedded <- list()
  left_out <- list()
  n_plots <- 0L
  for (f in found$files) {
    size <- file.size(f$path) * 4 / 3  # base64
    why <- if (identical(f$kind, "plot") && n_plots >= HUB_REPORT_MAX_PLOTS) "too many plots for one page"
      else if (!identical(f$kind, "plot") && size > per_file) paste0("bigger than ", HUB_REPORT_FILE_MB, " MB")
      else if (size > budget) paste0("the page would go over ", HUB_REPORT_MAX_MB, " MB")
      else NULL
    if (is.null(why)) {
      f$uri <- hub_data_uri(f$path)
      f$id <- paste0("f", length(embedded) + 1L)
      budget <- budget - nchar(f$uri)
      if (identical(f$kind, "plot")) n_plots <- n_plots + 1L
      embedded[[length(embedded) + 1]] <- f
    } else {
      left_out[[length(left_out) + 1]] <- list(label = f$label, path = f$path, size_mb = round(file.size(f$path) / 1024^2, 1), why = why)
    }
  }
  plots <- Filter(function(f) identical(f$kind, "plot"), embedded)
  datas <- Filter(function(f) identical(f$kind, "data"), embedded)

  figure <- function(f) paste0("<figure><img id=\"", f$id, "\" src=\"", f$uri, "\" alt=\"", hub_html(f$label), "\"><figcaption>",
                               hub_html(f$label), " <a class=\"dl\" data-img=\"", f$id, "\" download=\"", hub_html(basename(f$path)),
                               "\" href=\"#\">Download PNG</a></figcaption></figure>")
  file_link <- function(f) paste0("<li><a download=\"", hub_html(f$label), "\" href=\"", f$uri, "\">", hub_html(f$label), "</a> <span class=\"muted\">(",
                                  format(max(1, round(file.size(f$path) / 1024, 0)), big.mark = ","), " KB)</span></li>")
  links_for <- function(section, comparison = NULL) {
    fs <- Filter(function(f) identical(f$section, section) && (is.null(comparison) || identical(f$comparison, comparison)), datas)
    if (!length(fs)) "" else paste0("<ul class=\"files\">", paste(vapply(fs, file_link, ""), collapse = ""), "</ul>")
  }

  # ---- QC -----------------------------------------------------------------------------------
  qc_df <- hub_report_qc_table(qc$samples)
  qc_csv <- if (NROW(qc_df)) {
    tmp <- tempfile(fileext = ".csv"); utils::write.csv(qc_df, tmp, row.names = FALSE, na = ""); tmp
  } else NULL
  qc_html <- c("<section><h2>Quality control</h2>",
               "<p class=\"muted\">Numbers come straight from each tool's own report. FastQC's PASS/WARN/FAIL are FastQC's own judgements; ",
               "for RNA-seq a FAIL on <em>Per base sequence content</em> is common and usually expected. Whether a sample is good enough is the researcher's call (ask BSR if unsure).</p>",
               hub_html_table(qc_df, max_rows = 500L, digits = 1L),
               if (!is.null(qc_csv)) paste0("<ul class=\"files\"><li><a download=\"", hub_html(p$name), "_qc.csv\" href=\"", hub_data_uri(qc_csv), "\">", hub_html(p$name), "_qc.csv</a></li></ul>") else "",
               links_for("counts"), "</section>")

  # ---- DESeq2 -------------------------------------------------------------------------------
  deg_summary <- list()
  de_html <- character(0)
  for (cmp in found$comparisons) {
    sig <- hub_read_table(file.path(p$data_dir, "deseq2", paste0("DEG_fdr0.05_", cmp, ".txt")))
    lfc <- if (!is.null(sig) && "log2FoldChange" %in% names(sig)) suppressWarnings(as.numeric(sig$log2FoldChange)) else numeric(0)
    deg_summary[[cmp]] <- list(significant = length(lfc), up = sum(lfc > 0, na.rm = TRUE), down = sum(lfc < 0, na.rm = TRUE))
    top <- if (!is.null(sig) && "padj" %in% names(sig)) sig[order(suppressWarnings(as.numeric(sig$padj))), intersect(c("gene", "baseMean", "log2FoldChange", "padj"), names(sig)), drop = FALSE] else NULL
    figs <- Filter(function(f) identical(f$section, "deseq2") && identical(f$comparison, cmp), plots)
    de_html <- c(de_html, paste0("<h3>", hub_html(cmp), "</h3>"),
                 paste0("<p><strong>", deg_summary[[cmp]]$significant, "</strong> genes in DESeq2's FDR &lt; 0.05 list (",
                        deg_summary[[cmp]]$up, " up, ", deg_summary[[cmp]]$down, " down; direction from the log2 fold change).</p>"),
                 if (length(figs)) paste0("<div class=\"grid\">", paste(vapply(figs, figure, ""), collapse = ""), "</div>") else "",
                 if (!is.null(top)) paste0("<h4>Top genes by adjusted p-value</h4>", hub_html_table(top, 15L)) else "",
                 links_for("deseq2", cmp))
  }
  other_figs <- Filter(function(f) identical(f$section, "deseq2") && is.na(f$comparison), plots)
  if (length(other_figs)) de_html <- c(paste0("<h3>All samples</h3><div class=\"grid\">", paste(vapply(other_figs, figure, ""), collapse = ""), "</div>"), de_html)
  shared <- Filter(function(f) identical(f$section, "deseq2") && is.na(f$comparison), datas)
  if (length(shared)) de_html <- c(de_html, paste0("<ul class=\"files\">", paste(vapply(shared, file_link, ""), collapse = ""), "</ul>"))
  pipeline <- tryCatch(hub_read_json(hub_spec_path(p))$pipeline, error = function(e) NULL) %||% "rnaseq"
  peaks_run <- !identical(pipeline, "rnaseq")
  de_html <- if (length(de_html)) c("<section><h2>Differential expression (DESeq2)</h2>", de_html, "</section>")
             else if (peaks_run) character(0)
             else "<section><h2>Differential expression (DESeq2)</h2><p class=\"muted\">No DESeq2 results yet.</p></section>"

  # ---- DiffBind (ATAC-seq / ChIP-seq) ---------------------------------------------------------
  peak_summary <- list()
  db_html <- character(0)
  for (slug in found$diffbind) {
    tab_path <- list.files(file.path(p$data_dir, "diffbind", slug), pattern = "^DifferentialPeaks.*\\.txt$", full.names = TRUE)[1]
    tab <- if (!is.na(tab_path)) hub_read_table(tab_path) else NULL
    fold <- if (!is.null(tab) && "Fold" %in% names(tab)) suppressWarnings(as.numeric(tab$Fold)) else numeric(0)
    peak_summary[[slug]] <- list(differential_peaks = if (!is.null(tab)) NROW(tab) else NA,
                                 gained = sum(fold > 0, na.rm = TRUE), lost = sum(fold < 0, na.rm = TRUE))
    figs <- Filter(function(f) identical(f$section, "diffbind") && identical(f$comparison, slug), plots)
    top <- if (!is.null(tab) && "FDR" %in% names(tab)) tab[order(suppressWarnings(as.numeric(tab$FDR))),
      intersect(c("seqnames", "start", "end", "Fold", "p.value", "FDR"), names(tab)), drop = FALSE] else NULL
    db_html <- c(db_html, paste0("<h3>", hub_html(slug), "</h3>"),
                 if (!is.null(tab)) paste0("<p><strong>", NROW(tab), "</strong> regions in DiffBind's report (", peak_summary[[slug]]$gained,
                                          " with higher signal in the comparison group, ", peak_summary[[slug]]$lost, " lower; direction from DiffBind's Fold).</p>")
                 else "<p class=\"muted\">No DiffBind table yet.</p>",
                 if (length(figs)) paste0("<div class=\"grid\">", paste(vapply(figs, figure, ""), collapse = ""), "</div>") else "",
                 if (!is.null(top)) paste0("<h4>Top regions by FDR</h4>", hub_html_table(top, 15L)) else "",
                 links_for("diffbind", slug))
  }
  if (length(db_html)) db_html <- c("<section><h2>Differential peaks (DiffBind)</h2>", db_html, "</section>")
  else if (peaks_run) db_html <- "<section><h2>Differential peaks (DiffBind)</h2><p class=\"muted\">No DiffBind results yet.</p></section>"
  ann <- file.path(p$data_dir, "peak_annotation", "peak_annotation_summary.tsv")
  if (file.exists(ann)) {
    ann_tab <- hub_read_table(ann)
    db_html <- c(db_html, "<section><h2>Peak annotation (HOMER)</h2>", hub_html_table(ann_tab, 30L), links_for("annotation"), "</section>")
  }

  # ---- GSEA ---------------------------------------------------------------------------------
  gsea_top <- list()
  gsea_html <- character(0)
  for (f in found$gsea_reports) {
    tab <- hub_read_table(f)
    key <- paste(basename(dirname(dirname(f))), basename(dirname(f)), sep = " / ")
    if (!is.null(tab) && all(c("Term", "NES", "FDR q-val") %in% names(tab))) {
      tab <- tab[order(suppressWarnings(as.numeric(tab[["FDR q-val"]]))), c("Term", "NES", "FDR q-val"), drop = FALSE]
      gsea_top[[key]] <- utils::head(tab, 5)
      gsea_html <- c(gsea_html, paste0("<h3>", hub_html(key), "</h3>"), hub_html_table(tab, 10L))
    }
  }
  gfigs <- Filter(function(f) identical(f$section, "gsea"), plots)
  if (length(gfigs)) gsea_html <- c(gsea_html, paste0("<div class=\"grid\">", paste(vapply(gfigs, figure, ""), collapse = ""), "</div>"))
  if (length(gsea_html)) gsea_html <- c("<section><h2>Pathways (GSEA)</h2>", gsea_html, links_for("gsea"), "</section>")

  # ---- steps, methods, files left out ------------------------------------------------------
  step_rows <- do.call(rbind, lapply(st$steps %||% list(), function(s) data.frame(step = s$step, state = s$state, stringsAsFactors = FALSE)))
  methods <- tryCatch(project_methods_text(p), error = function(e) "")
  tools <- tryCatch(tool_reference_summary(p), error = function(e) NULL)
  left_html <- if (length(left_out)) paste0(
    "<section><h2>Kept on the cluster (too big for this page)</h2><div class=\"scroll\"><table><tr><th>file</th><th>size</th><th>why</th><th>where</th></tr>",
    paste(vapply(left_out, function(x) paste0("<tr><td>", hub_html(x$label), "</td><td>", x$size_mb, " MB</td><td>", hub_html(x$why), "</td><td><code>", hub_html(x$path), "</code></td></tr>"), ""), collapse = ""),
    "</table></div></section>") else ""

  page <- c(
    "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    paste0("<title>", hub_html(p$name), " report</title>"),
    "<style>",
    ":root{--bg:#fff;--fg:#1d2329;--muted:#5d6873;--line:#d9dee3;--accent:#1f6feb;--card:#f6f8fa}",
    "@media (prefers-color-scheme:dark){:root{--bg:#14181c;--fg:#e6e9ec;--muted:#9aa5b0;--line:#2c343c;--accent:#6ea8fe;--card:#1b2127}}",
    "body{background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,Segoe UI,sans-serif;margin:0;padding:24px 16px}",
    "main{max-width:1100px;margin:0 auto}h1{margin:0 0 4px}h2{border-bottom:1px solid var(--line);padding-bottom:4px;margin-top:36px}",
    ".muted{color:var(--muted)}a{color:var(--accent)}code{font-size:12px;word-break:break-all}",
    ".scroll{overflow-x:auto}table{border-collapse:collapse;font-size:13px;margin:8px 0}th,td{border:1px solid var(--line);padding:4px 8px;text-align:left;white-space:nowrap}th{background:var(--card)}",
    ".grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(300px,1fr));gap:16px}figure{margin:0;background:var(--card);border:1px solid var(--line);border-radius:8px;padding:8px}",
    "figure img{width:100%;height:auto;background:#fff;border-radius:4px}figcaption{font-size:13px;margin-top:6px;word-break:break-all}",
    ".files{padding-left:20px}.facts{display:flex;flex-wrap:wrap;gap:8px 24px;margin:12px 0}.facts div{background:var(--card);border:1px solid var(--line);border-radius:6px;padding:6px 10px}",
    "pre{white-space:pre-wrap;background:var(--card);border:1px solid var(--line);border-radius:6px;padding:12px;font-size:13px}",
    "</style></head><body><main>",
    paste0("<h1>", hub_html(p$name), "</h1>"),
    paste0("<p class=\"muted\">CodeSpringApp ", hub_html(switch(pipeline, rnaseq = "bulk RNA-seq", atac = "bulk ATAC-seq", chip = "bulk ChIP-seq", pipeline)), " run report, made ", hub_html(format(Sys.time(), "%Y-%m-%d %H:%M")), ".</p>"),
    "<div class=\"facts\">",
    paste0("<div>Run: <strong>", hub_html(st$state %||% "unknown"), "</strong></div>"),
    paste0("<div>Samples: <strong>", length(qc$samples), "</strong></div>"),
    paste0("<div>Reference: <strong>", hub_html(tryCatch(project_reference_label(p), error = function(e) "")), "</strong></div>"),
    "</div>",
    if (!identical(st$state, "complete")) "<p><strong>This run hasn't finished, so the report shows only the steps done so far.</strong></p>" else "",
    paste0("<p>Full results, interactive plots and every file are in CodeSpringApp: choose RNA-seq, then project <strong>", hub_html(p$name),
           "</strong>, then the Results Explorer tab (or ask Claude to open the results). On the cluster: <code>", hub_html(project_result_dir(p)), "</code></p>"),
    "<section><h2>Steps</h2>", hub_html_table(step_rows, 50L), "</section>",
    qc_html, de_html, db_html, gsea_html,
    "<section><h2>Methods</h2><p class=\"muted\">Generated by CodeSpringApp from the steps that ran; check versions and settings before publishing.</p>",
    paste0("<pre>", hub_html(methods), "</pre>"),
    if (!is.null(tools) && NROW(tools)) hub_html_table(tools, 50L) else "", "</section>",
    left_html,
    "<script>document.querySelectorAll('a[data-img]').forEach(function(a){a.href=document.getElementById(a.dataset.img).src;});</script>",
    "</main></body></html>")

  out <- hub_report_path(p)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(out, ".tmp")
  writeLines(page, tmp, useBytes = TRUE)
  file.rename(tmp, out)
  if (!is.null(qc_csv)) unlink(qc_csv)
  list(ok = TRUE, project_id = p$id, state = st$state, report = out, size_mb = round(file.size(out) / 1024^2, 2),
       plots = lapply(plots, function(f) list(label = f$label, path = f$path, section = f$section, comparison = if (is.na(f$comparison)) NULL else f$comparison)),
       files = lapply(datas, function(f) f$label), left_out = left_out,
       deg_summary = deg_summary, gsea_top = gsea_top, peak_summary = peak_summary, pipeline = pipeline, qc_summary = qc$summary,
       note = paste("Everything in the report comes from CodeSpringApp's own outputs; nothing is recalculated.",
                    "Plots and small files are embedded and downloadable; big files are listed with their cluster paths."))
}
