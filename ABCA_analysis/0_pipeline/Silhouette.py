import math
import seaborn as sns
import numpy as np
import scanpy as sc
import pandas as pd
from kneed import KneeLocator
from pathlib import Path
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

# Selecting the brain region
select_region = "Isocortex-1"

# Load AnnData object
base_path = Path("/data/scRNA/ABCA/AIBS/AWS/expression_matrices/WMB-10Xv3/20230630/")
expr_path = base_path / f"WMB-10Xv3-{select_region}-raw-wmeta.h5ad"
adata = sc.read_h5ad(expr_path)

# Filter out low-count cell types
n_cells_threshold = 100
class_counts = adata.obs["class"].value_counts()
classes_to_remove = class_counts[class_counts < n_cells_threshold].index.tolist()
adata = adata[~adata.obs["class"].isin(classes_to_remove)].copy()

# Normalize and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # Store raw normalized data

# HVG selection
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)

# PCA
sc.tl.pca(adata, n_comps=11)

# Neighbors + UMAP
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Silhouette score output
scores = []

# Sweep resolutions from 0.001 to 1.0 with interval 0.002
resolutions = np.arange(0.001, 1.002, 0.002)

# Compute silhouette scores for each resolution
for res in resolutions:
    key = f'leiden_{res:.3f}'
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, key_added=key, resolution=res)

    # Subsample
    n = 10000
    idx = np.random.choice(adata.n_obs, size=n, replace=False)
    X_sub = adata.obsm["X_pca"][idx]
    labels_sub = adata.obs[key].values[idx]

    try:
        score = silhouette_score(X_sub, labels_sub)
    except:
        score = np.nan  # in case of errors (e.g., only one cluster)

    print(f"Resolution {res:.3f} → Silhouette Score: {score:.3f}")
    scores.append({'resolution': res, 'silhouette_score': score})

# Save to CSV
scores_df = pd.DataFrame(scores)
scores_df.to_csv("silhouette_scores_by_resolution.csv", index=False)

# df = pd.read_csv("silhouette_scores_by_resolution.csv")
# plt.figure(figsize=(10, 5))
# sns.lineplot(data=df, x="resolution", y="silhouette_score", marker="o")
# plt.title("Silhouette Score vs Leiden Resolution")
# plt.xlabel("Resolution")
# plt.ylabel("Silhouette Score")
# plt.grid(True)
# plt.tight_layout()
# plt.show()