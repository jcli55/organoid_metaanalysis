import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv

# This is part 2 of the multivelo script after running seurat_wnn.R
# Part 1 is in run_multivelo.py
# Adapted to run on the transition cells found by celltran - 04/25/25

# Part 2 ------------------------------------
adata_rna = sc.read_h5ad("/dfs3b/ruic20_lab/jeancl2/data/celltran/velocity/adata_rna_normalized_transition_cells.h5ad")
adata_atac = sc.read_h5ad("/dfs3b/ruic20_lab/jeancl2/data/celltran/velocity/adata_atac_normalized_transition_cells.h5ad")

# Read in Seurat WNN neighbors.
nn_idx = np.loadtxt("/dfs3b/ruic20_lab/jeancl2/data/celltran/velocity/seurat_wnn/nn_idx.txt", delimiter=',')
nn_dist = np.loadtxt("/dfs3b/ruic20_lab/jeancl2/data/celltran/velocity/seurat_wnn/nn_dist.txt", delimiter=',')

# Make sure cell names match.
#nn_cells = pd.Index(pd.read_csv("/dfs3b/ruic20_lab/jeancl2/data/celltran/velocity/seurat_wnn/nn_cells.txt", header=None)[0])
#np.all(nn_cells == adata_atac.obs_names)

mv.knn_smooth_chrom(adata_atac, nn_idx, nn_dist)

# Running Multi-omic Dynamic Model - This will take a while. Parallelization is high recommended.
# mv.settings.VERBOSITY = 0

adata_result = mv.recover_dynamics_chrom(adata_rna, 
                                         adata_atac, 
                                         max_iter=5, 
                                         init_mode="invert",
                                         parallel=True,
                                         n_jobs=20,
                                         save_plot=False,
                                         rna_only=False,
                                         fit=True,
                                         n_anchors=500, 
                                        )

adata_result.write("/dfs3b/ruic20_lab/jeancl2/data/celltran/velocity/multivelo_result.h5ad")

# Transfer some metadata labels to the object
metadata = pd.read_csv('/dfs3b/ruic20_lab/singlecell/jeanl/organoid/csv/metadata_clean.csv', index_col = 0)
class_list = []
subclass_list = []
age_list = []
for barcode in adata_result.obs_names:
    class_list.append(metadata.loc[barcode, 'majorclass'])
    subclass_list.append(metadata.loc[barcode, 'subclass'])
    age_list.append(metadata.loc[barcode, 'age'])
adata_result.obs['class'] = class_list
adata_result.obs['subclass'] = subclass_list
adata_result.obs['age'] = age_list

# # Moved these from set_macrostates_scvi.py - don't need here
# adata_result.obs["maturityclass"] = adata_result.obs.subclass.replace(
#     {
#         "ML_Cone": "Cone",
#         "S_Cone": "Cone",
#         "OFF-BC": "BC",
#         "ON-BC": "BC",
#         "RBC": "BC",
#         "GABAergic": "AC",
#         "Glycinergic": "AC",
#         "SACs": "AC",
#         "dual ACs": "AC",
#         "HC0": "HC",
#         "HC1": "HC",
#         "OFF_MGC": "RGC",
#         "OFF_PGC": "RGC",
#         "ON_MGC": "RGC",
#         "ON_PGC": "RGC"
#     }
# )

# adata_result.obs['terminalstates'] = adata_result.obs.maturityclass.replace(
#     {
#         "Cone Precursor": np.nan,
#         "Rod Precursor": np.nan,
#         "BC Precursor": np.nan,
#         "AC Precursor": np.nan,
#         "HC Precursor": np.nan,
#         "PRPC": np.nan,
#         "NRPC": np.nan
#     }
# )

# Effectively renaming velo_s to velocity so cellrank can run on the multivelo object
adata_result.layers['velocity'] = adata_result.layers['velo_s']

# Plot velocity stream and latent time
scv.pp.neighbors(adata_result)

mv.velocity_graph(adata_result)
mv.latent_time(adata_result)

# Save the result for use later on
adata_result.write("/dfs3b/ruic20_lab/jeancl2/data/celltran/velocity/multivelo_result.h5ad")

mv.velocity_embedding_stream(adata_result, basis='umap', color='class', save='velocity_stream.png')
scv.pl.scatter(adata_result, color=['latent_time', 'velo_s_norm_pseudotime'], color_map='gnuplot', size=80, save='latent_time.png')