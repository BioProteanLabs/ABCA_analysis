import math
import seaborn as sns
import numpy as np
import scanpy as sc
import pandas as pd
from kneed import KneeLocator
from pathlib import Path
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

# Selecting the brain region
select_region = "Isocortex-1"

# Loading AnnData object
base_path = Path("/data/scRNA/ABCA/AIBS/AWS/expression_matrices/WMB-10Xv3/20230630/")
expr_path = base_path / f"WMB-10Xv3-{select_region}-raw-wmeta.h5ad"
adata = sc.read_h5ad(expr_path)

# Count the number of cells in each class
class_counts = adata.obs["class"].value_counts()

# Identify classes with few cells
n_cells_threshold = 100

# Filter out the cells belonging to the identified classes
classes_to_remove = class_counts[class_counts < n_cells_threshold].index.tolist()
adata = adata[~adata.obs["class"].isin(classes_to_remove)].copy()

gene_names_df = adata.var.copy()
adata.var.set_index("gene_symbol", inplace=True)

# Preprocess the data
sc.pp.normalize_total(adata, target_sum=1e4)  # Normalize counts per cell
sc.pp.log1p(adata)  # Log-transform the data
adata.raw = adata # Save raw log-normalized data

sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)  # Select highly variable genes


# Assuming `adata` has run PCA
sc.tl.pca(adata)
variance_ratio = adata.uns['pca']['variance_ratio']
x = range(1, len(variance_ratio)+1)
knee = KneeLocator(x, variance_ratio, curve='convex', direction='decreasing')
sc.tl.pca(adata, n_comps=knee.knee)

sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Using the igraph implementation and a fixed number of iterations can be faster
for res in [0.006,0.2]:
    #Clustering
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, key_added=f'leiden_{res}', resolution=res)
    
    # Perform differential expression analysis
    sc.tl.rank_genes_groups(
        adata,
        key_added=f"rank_genes_leiden_{res}",
        groupby=f'leiden_{res}',
        method="wilcoxon",
        use_raw=True
    )

    # Extract the differential expression results
    deg_results = pd.DataFrame({
        group: pd.DataFrame(adata.uns[f"rank_genes_leiden_{res}"]["names"])[group]
        for group in adata.uns[f"rank_genes_leiden_{res}"]["names"].dtype.names
    })

    # Save the dataframe to a CSV file
    output_path = base_path / "outputs/DEG/" / f"WMB-10Xv3-{select_region}-DEG-leiden_{res}.csv"
    deg_results.to_csv(output_path, index=False)