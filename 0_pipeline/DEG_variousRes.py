import math
import seaborn as sns
import numpy as np
import scanpy as sc
import pandas as pd
from kneed import KneeLocator
from pathlib import Path
import matplotlib.pyplot as plt
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

# Loading AnnData object
select_region = "Isocortex-1"
base_path = Path("/data/scRNA/ABCA/AIBS/AWS/expression_matrices/WMB-10Xv3/20230630/")
expr_path = base_path / f"AnnData/WMB-10Xv3-{select_region}-raw-wmeta-DEGpp.h5ad"
adata = sc.read_h5ad(expr_path)

# Create a DataFrame to hold cluster labels
cluster_labels_df = pd.DataFrame(index=adata.obs_names)

# Sweep resolutions from 0.001 to 1.0 with interval 0.002
resolutions = np.arange(0.001, 1.002, 0.002)

for res in resolutions:
    key = f'leiden_{res:.3f}'
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, key_added=key, resolution=res)

    # Save the cluster labels for this resolution
    cluster_labels_df[key] = adata.obs[key].values

    print(f"Finished clustering for resolution {res:.3f}")

# Save cluster labels to CSV
cluster_labels_df.to_csv("leiden_cluster_labels_by_resolution.csv")