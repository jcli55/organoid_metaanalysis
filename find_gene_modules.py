import scanpy as sc
import hotspot
import numpy as np
import pandas as pd
import joblib

# Find gene modules using Hotspot
# Following basic workflow: https://hotspot.readthedocs.io/en/latest/index.html

# Load object and set output directory
adata = sc.read_h5ad("/dfs3b/ruic20_lab/jeancl2/data/chen_ro_clean.h5ad")
adata = adata[adata.obs.majorclass=="PRPC"].copy()

output = "/dfs3b/ruic20_lab/jeancl2/data/prpc_gene_modules/"

# Subset to highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', subset=True)

# If the matrix is missing total_counts, calculate it again using raw counts assuming data is not normalized
counts = adata.X.toarray()
total_counts = np.sum(counts, axis=1)
adata.obs["total_counts"] = total_counts

# Double check and remove genes with zero variance
gene_variances = np.var(counts, axis=0) 
zero_variance_genes = adata.var_names[gene_variances == 0]
adata = adata[:, ~adata.var_names.isin(zero_variance_genes)].copy()

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log_normalized"] = adata.X.copy()
sc.pp.scale(adata)
sc.tl.pca(adata)

adata.write(f"{output}hotspot_adata.h5ad")

adata.layers["counts_csc"] = adata.layers["counts"].tocsc()
hs = hotspot.Hotspot(
    adata,
    layer_key="counts_csc",
    model='danb',
    latent_obsm_key="X_pca",
    umi_counts_obs_key="total_counts"
)

hs.create_knn_graph(weighted_graph=False, n_neighbors=30)

hs_results = hs.compute_autocorrelations()
hs_results['gene'] = adata.var_names
hs_results.to_csv(f"{output}hs_results.csv", index=False)

hs_genes = hs_results.loc[hs_results.FDR < 0.05].index # Select genes

local_correlations = hs.compute_local_correlations(hs_genes, jobs=4) # jobs for parallelization

modules = hs.create_modules(
    min_gene_threshold=30, core_only=True, fdr_threshold=0.05
)

joblib.dump(hs, f"{output}hotspot_object.h5ad")

# Still need to find a way to save the plot
#hs.plot_local_correlations(vmin=-12, vmax=12)

module_scores = hs.calculate_module_scores()
module_scores.to_csv(f"{output}module_scores.csv")

module_cols = []
for c in module_scores.columns:
    key = f"Module{c}"
    adata.obs[key] = module_scores[c]
    module_cols.append(key)

sc.pl.umap(adata, color=module_cols, frameon=False, vmin=-1, vmax=1, save="_gene_modules.png")
