import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv

# This is part 2 of the multivelo script after running seurat_wnn.R
# Part 1 is in run_multivelo.py
# By Jean Li 12/11/2023
# Updated to run on hpc3 and on the full data - 2/13/25

# Part 2 ------------------------------------
adata_rna = sc.read_h5ad("/dfs3b/ruic20_lab/jeancl2/data/multivelo_full/adata_rna_normalized_full_highly_variable.h5ad")
adata_atac = sc.read_h5ad("/dfs3b/ruic20_lab/jeancl2/data/multivelo_full/adata_atac_normalized_full_highly_variable.h5ad")

# Read in Seurat WNN neighbors.
nn_idx = np.loadtxt("/dfs3b/ruic20_lab/jeancl2/data/multivelo_full/seurat_wnn/nn_idx.txt", delimiter=',')
nn_dist = np.loadtxt("/dfs3b/ruic20_lab/jeancl2/data/multivelo_full/seurat_wnn/nn_dist.txt", delimiter=',')

# Make sure cell names match.
#nn_cells = pd.Index(pd.read_csv("/dfs3b/ruic20_lab/jeancl2/data/multivelo_full/seurat_wnn/nn_cells.txt", header=None)[0])
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

adata_result.write("/dfs3b/ruic20_lab/jeancl2/data/multivelo_full/multivelo_result.h5ad")

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

# Effectively renaming velo_s to velocity so cellrank can run on the multivelo object
adata_result.layers['velocity'] = adata_result.layers['velo_s']

# Plot velocity stream and latent time
scv.pp.neighbors(adata_result)

mv.velocity_graph(adata_result)
mv.latent_time(adata_result)

mv.velocity_embedding_stream(adata_result, basis='umap', color='majorclass', save='velocity_stream.png')
scv.pl.scatter(adata_result, color='latent_time', color_map='gnuplot', size=80, save='latent_time.png')

# Save the result for use later on
adata_result.write("/dfs3b/ruic20_lab/jeancl2/data/multivelo_full/multivelo_result.h5ad")
