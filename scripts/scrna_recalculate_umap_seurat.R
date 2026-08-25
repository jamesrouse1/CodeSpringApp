args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6L) stop("Usage: scrna_recalculate_umap_seurat.R object.rds cells.tsv output.tsv n_neighbors min_dist seed")

object_path <- args[[1]]
cells_path <- args[[2]]
output_path <- args[[3]]
n_neighbors <- suppressWarnings(as.integer(args[[4]]))
min_dist <- suppressWarnings(as.numeric(args[[5]]))
seed <- suppressWarnings(as.integer(args[[6]]))
if (is.na(n_neighbors) || n_neighbors < 2L) n_neighbors <- 15L
if (!is.finite(min_dist) || min_dist < 0) min_dist <- 0.3
if (is.na(seed)) seed <- 1234L

for (package in c("Seurat", "uwot")) if (!requireNamespace(package, quietly = TRUE)) stop("Missing R package: ", package)
saved <- readRDS(object_path)
object <- if (inherits(saved, "Seurat")) saved else if (is.list(saved) && inherits(saved$object, "Seurat")) saved$object else NULL
if (is.null(object)) stop("The RDS does not contain a Seurat object.")

requested <- utils::read.delim(cells_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!"cell" %in% names(requested)) stop("The included-cell table does not contain a cell column.")
cells <- intersect(as.character(requested$cell), colnames(object))
if (length(cells) < 10L) stop("At least 10 selected cells must be present in the Seurat object.")

reductions <- names(object@reductions)
reduction <- if ("harmony" %in% reductions) "harmony" else if ("pca" %in% reductions) "pca" else ""
if (!nzchar(reduction)) stop("The processed Seurat object has no PCA or Harmony representation.")
embedding <- Seurat::Embeddings(object, reduction = reduction)
embedding <- embedding[cells, seq_len(min(30L, ncol(embedding))), drop = FALSE]
n_neighbors <- min(n_neighbors, max(2L, nrow(embedding) - 1L))

set.seed(seed)
coordinates <- uwot::umap(
  embedding,
  n_neighbors = n_neighbors,
  min_dist = min_dist,
  n_components = 2L,
  metric = "euclidean",
  init = "spectral",
  n_threads = max(1L, suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))),
  ret_model = FALSE,
  verbose = FALSE,
  seed = seed
)
result <- data.frame(cell = rownames(embedding), UMAP_1 = coordinates[, 1], UMAP_2 = coordinates[, 2], check.names = FALSE)
temporary <- paste0(output_path, ".tmp")
utils::write.table(result, temporary, sep = "\t", row.names = FALSE, quote = FALSE)
if (!file.rename(temporary, output_path)) stop("Could not finalize the selected-cell UMAP table.")
