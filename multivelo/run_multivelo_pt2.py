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

# Part 2 ------------------------------------
adata_rna = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/adata_rna_normalized_downsampled.h5ad")
adata_atac = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/adata_atac_normalized_downsampled.h5ad")

# Read in Seurat WNN neighbors.
nn_idx = np.loadtxt("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/multivelo/seurat_wnn/nn_idx.txt", delimiter=',')
nn_dist = np.loadtxt("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/multivelo/seurat_wnn/nn_dist.txt", delimiter=',')
nn_cells = pd.Index(pd.read_csv("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/multivelo/seurat_wnn/nn_cells.txt", header=None)[0])

# Make sure cell names match.
#np.all(nn_cells == adata_atac.obs_names)

mv.knn_smooth_chrom(adata_atac, nn_idx, nn_dist)

# Running Multi-omic Dynamic Model - This will take a while. Parallelization is high recommended.
# mv.settings.VERBOSITY = 0

adata_result = mv.recover_dynamics_chrom(adata_rna, 
                                         adata_atac, 
                                         max_iter=5, 
                                         init_mode="invert",
                                         parallel=True,
                                         save_plot=False,
                                         rna_only=False,
                                         fit=True,
                                         n_anchors=500, 
                                        )

# Transfer some metadata labels to the object
metadata = pd.read_csv('/storage/singlecell/jeanl/organoid/csv/metadata_clean.csv', index_col = 0)
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

# Plot velocity stream and latent time
scv.pp.neighbors(adata_result)

mv.velocity_graph(adata_result)
mv.latent_time(adata_result)

mv.velocity_embedding_stream(adata_result, basis='umap', color='majorclass', save='velocity_stream.png')
scv.pl.scatter(adata_result, color='latent_time', color_map='gnuplot', size=80, save='latent_time.png')

# Save the result for use later on
adata_result.write("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/multivelo/multivelo_result.h5ad")