#!/usr/bin/env python3
import os
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


def main():
    if len(sys.argv) < 7:
        raise SystemExit("Usage: scrna_recalculate_umap_scanpy.py object.h5ad cells.tsv output.tsv n_neighbors min_dist seed")
    object_path, cells_path, output_path = sys.argv[1:4]
    n_neighbors = max(2, int(sys.argv[4]))
    min_dist = max(0.0, float(sys.argv[5]))
    seed = int(sys.argv[6])

    requested = pd.read_csv(cells_path, sep="\t", dtype=str)
    if "cell" not in requested.columns:
        raise SystemExit("The included-cell table does not contain a cell column.")
    source = ad.read_h5ad(object_path, backed="r")
    available = set(map(str, source.obs_names))
    cells = [cell for cell in requested["cell"].astype(str) if cell in available]
    if len(cells) < 10:
        raise SystemExit("At least 10 selected cells must be present in the Scanpy object.")
    subset = source[cells].to_memory()

    integration = str(subset.uns.get("codespring_integration", "none")).lower()
    preferred = {"scvi": "X_scVI", "harmony": "X_harmony"}.get(integration, "X_pca")
    representation = preferred if preferred in subset.obsm else next((key for key in ("X_harmony", "X_scVI", "X_pca") if key in subset.obsm), "")
    if not representation:
        raise SystemExit("The processed Scanpy object has no PCA, Harmony, or scVI representation.")

    n_neighbors = min(n_neighbors, max(2, subset.n_obs - 1))
    n_pcs = min(30, subset.obsm[representation].shape[1])
    sc.pp.neighbors(subset, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=representation, metric="euclidean", random_state=seed)
    sc.tl.umap(subset, min_dist=min_dist, spread=1.0, init_pos="spectral", random_state=seed)
    result = pd.DataFrame({"cell": subset.obs_names.astype(str), "UMAP_1": subset.obsm["X_umap"][:, 0], "UMAP_2": subset.obsm["X_umap"][:, 1]})
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    result.to_csv(temporary, sep="\t", index=False)
    os.replace(temporary, output)


if __name__ == "__main__":
    main()
