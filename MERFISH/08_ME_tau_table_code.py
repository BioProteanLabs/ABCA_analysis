from scipy.stats import pearsonr
import seaborn as sns
import numpy as np
import scanpy as sc
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
import concurrent.futures

modules_df = pd.read_csv(
    "/data/scRNA/ABCA/AIBS/AWS/expression_matrices/WMB-10Xv3/20230630/outputs/WMB-10Xv3-Isocortex-1-raw-sc-wgcna-modules.csv",
    index_col = "Unnamed: 0")

merfish_panel_df = pd.read_csv(
    "/data/scRNA/ABCA/AIBS/AWS/expression_matrices/MERFISH-C57BL6J-638850/20230830/gene_panel.csv"
)

# Getting all module names
all_modules = modules_df["color"].unique().tolist()

# Function to compute Tau score from an AnnData object that has PCA computed
def compute_tau(adata_obj, module_obs_key):
    # Create a DataFrame with cell class and the module eigengene value
    tau_df = adata_obj.obs[['class']].copy()
    tau_df[module_obs_key] = adata_obj.obsm['X_pca'][:, 0]
    # Shift the values to be positive
    tau_df[module_obs_key] = tau_df[module_obs_key] - tau_df[module_obs_key].min() + 1e-6
    # Compute average expression per class
    avg_expression = tau_df.groupby('class')[module_obs_key].mean()
    max_expr = avg_expression.max()
    normalized_expr = avg_expression / max_expr
    # Compute Tau: tau = sum(1 - normalized_expr)/(number_of_classes - 1)
    tau = np.sum(1 - normalized_expr) / (len(normalized_expr) - 1)
    return tau

# Define a function for one permutation iteration.
def permutation_tau(seed):
    # Set the random seed for reproducibility in this process.
    np.random.seed(seed)
    # Randomly sample genes without replacement
    random_genes = np.random.choice(all_genes, size=n_genes_in_module, replace=False)
    # Create a copy of adata for the selected genes
    adata_perm = adata[:, random_genes].copy()
    
    # Preprocess the subset (normalization, log transform, scaling, PCA)
    sc.pp.normalize_total(adata_perm, target_sum=1e4)
    sc.pp.log1p(adata_perm)
    sc.pp.scale(adata_perm, max_value=10)
    sc.tl.pca(adata_perm, n_comps=2, svd_solver="arpack")
    
    # Store the first principal component as the temporary module eigengene
    temp_key = "temp_module"
    adata_perm.obs[temp_key] = adata_perm.obsm['X_pca'][:, 0]
    
    # Compute and return the Tau score for this random gene set.
    return compute_tau(adata_perm, temp_key)

# Parameters
n_permutations = 1000  # number of random permutations

'''
MERFISH data
'''

expr_path = '/data/scRNA/ABCA/AIBS/AWS/expression_matrices/MERFISH-C57BL6J-638850/20230830/C57BL6J-638850-raw-wmeta.h5ad'
adata = sc.read_h5ad(expr_path)

# Prepare a list to hold results for each module.
merfish_results = []

for module_name in set(all_modules) - {"grey"}:

    # Observed: Select the genes that belong to your module
    select_genes = modules_df[modules_df["module"] == module_name].index

    # Compute the observed Tau score for the whole module
    common_genes = adata.var_names.intersection(select_genes)
    adata_subset = adata[:, common_genes].copy()

    # Preprocess: normalization, log, scaling, PCA
    sc.pp.normalize_total(adata_subset, target_sum=1e4)
    sc.pp.log1p(adata_subset)
    sc.pp.scale(adata_subset, max_value=10)
    sc.tl.pca(adata_subset, n_comps=2, svd_solver="arpack")

    # Store the eigengene in obs with a key (here, we use the module name)
    adata_subset.obs[module_name] = adata_subset.obsm['X_pca'][:, 0]
    observed_tau = compute_tau(adata_subset, module_name)
    print(f"Observed Tau score for module {module_name}: {observed_tau:.3f}")

    # Get the list of all genes from the AnnData object
    all_genes = adata.var_names

    # Prepare random seeds for each permutation for reproducibility
    seeds = np.random.randint(0, 1000000, size=n_permutations)

    # Run the permutation iterations in parallel using ProcessPoolExecutor.
    with concurrent.futures.ProcessPoolExecutor() as executor:
        tau_null = list(executor.map(permutation_tau, seeds))

    tau_null = np.array(tau_null)

    # Compute the p-value: proportion of permutations with Tau score >= observed_tau.
    r = np.sum(tau_null >= observed_tau)
    p_value = (r + 1) / (n_permutations + 1)
    print(f"Permutation test p-value: {p_value:.5f}")

    merfish_results.append({
        "module": module_name,
        "n_common_genes": len(common_genes),
        "MERFISH_tau": observed_tau,
        "p_value": p_value
    })

merfish_results_df = pd.DataFrame(merfish_results)
merfish_results_df.to_csv("/data/scRNA/ABCA/AIBS/AWS/expression_matrices/WMB-10Xv3/20230630/outputs/C57BL6J-638850-raw-sc-wgcna-ME-merfish-Tau.csv")

del adata