import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv

# Runs part 1 of the multivelo pipeline
# Adapted to run on the reannotated data - 3/5/24 (old)
# Adapted to run on hpc3 without downsampling - 2/13/25

# Demo params --------------
# scv.settings.verbosity = 3
# scv.settings.presenter_view = True
# scv.set_figure_params('scvelo')
# pd.set_option('display.max_columns', 100)
# pd.set_option('display.max_rows', 200)
# np.set_printoptions(suppress=True)

# Part 1 -------------------------------------
# Read in and prepare the merged rna and atac data
# RNA
adata_rna = sc.read_h5ad("/dfs3b/ruic20_lab/singlecell/jeanl/organoid/data/Chen/h5ad/merged_spliced_unspliced.h5ad")
adata_rna.var_names_make_unique()
# Added this functionality in the merge_loom.py so no need for this in here anymore
# names = []
# for name in adata_rna.obs_names:
#     names.append(f'{name}-0')
# adata_rna.obs_names = names

# Transfer the cell labels to the spliced + unspliced data
filtered_annotated_rna = sc.read_h5ad("/dfs3b/ruic20_lab/singlecell/jeanl/organoid/data/Chen/h5ad/adata_rna.h5ad")

filtered_cells = pd.Index(np.intersect1d(adata_rna.obs_names, filtered_annotated_rna.obs_names))
adata_rna = adata_rna[filtered_cells]

celltypes = []
for name in adata_rna.obs_names:
    celltypes.append(filtered_annotated_rna.obs.loc[name]['majorclass'])
adata_rna.obs['majorclass'] = celltypes
adata_rna.obs['majorclass'] = adata_rna.obs['majorclass'].astype('category')

# ATAC
adata_atac = sc.read_h5ad("/dfs3b/ruic20_lab/singlecell/jeanl/organoid/data/Chen/h5ad/adata_atac.h5ad")

# Get shared cells and genes between rna and atac and subset on just these
shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))

adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]

# Normalization for rna and atac
scv.pp.normalize_per_cell(adata_rna)
scv.pp.log1p(adata_rna)
scv.pp.moments(adata_rna, n_pcs=30, n_neighbors=50) # Might need to move to after subsetting bc it will corrupt the neighbors graph

mv.tfidf_norm(adata_atac)

'''
# Subset cells and genes (to reduce runtime and memory requirements)
sc.pp.highly_variable_genes(adata_rna, flavor="seurat_v3", n_top_genes=2000, subset=True)
sc.pp.subsample(adata_rna, n_obs=10000)

# Another round of subset to get atac to match
shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))

adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]
'''

# Transfer the umap coordinates to the spliced/unspliced data
# Create a dictionary of sampleid/barcodes mapped to their umap coordinate
dict = {}
for i in range(0, len(filtered_annotated_rna)):
    dict[filtered_annotated_rna.obs_names[i]] = filtered_annotated_rna.obsm['X_umap'][i]

# Recreate the X_umap stacked array using the order the velocyto data was in
array = dict[adata_rna.obs_names[0]]
for i in range(1, len(adata_rna)):
    array = np.vstack([array, dict[adata_rna.obs_names[i]]])

# Set the umap
adata_rna.obsm['X_umap'] = array

adata_rna.write("/dfs3b/ruic20_lab/jeancl2/data/multivelo_full/adata_rna_normalized_full.h5ad")
adata_atac.write("/dfs3b/ruic20_lab/jeancl2/data/multivelo_full/adata_atac_normalized_full.h5ad")

# Generate the filtered cells txt for the R script
adata_rna.obs_names.to_frame().to_csv('/dfs3b/ruic20_lab/jeancl2/data/multivelo_full/seurat_wnn/filtered_cells.txt', header=False, index=False)

# -------------------------------------------
# Run Seurat WNN script to generate WNN files
# -------------------------------------------

# # Moved part 2 to its own script run_multivelo_pt2.py
# # Part 2 ------------------------------------
# adata_rna = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/adata_rna_normalized_downsampled.h5ad")
# adata_atac = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/adata_atac_normalized_downsampled.h5ad")

# # Read in Seurat WNN neighbors.
# nn_idx = np.loadtxt("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/multivelo/seurat_wnn/nn_idx.txt", delimiter=',')
# nn_dist = np.loadtxt("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/multivelo/seurat_wnn/nn_dist.txt", delimiter=',')
# nn_cells = pd.Index(pd.read_csv("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/multivelo/seurat_wnn/nn_cells.txt", header=None)[0])

# # Make sure cell names match.
# #np.all(nn_cells == adata_atac.obs_names)

# mv.knn_smooth_chrom(adata_atac, nn_idx, nn_dist)

# # Running Multi-omic Dynamic Model - This will take a while. Parallelization is high recommended.
# # mv.settings.VERBOSITY = 0

# adata_result = mv.recover_dynamics_chrom(adata_rna, 
#                                          adata_atac, 
#                                          max_iter=5, 
#                                          init_mode="invert",
#                                          parallel=True,
#                                          save_plot=False,
#                                          rna_only=False,
#                                          fit=True,
#                                          n_anchors=500, 
#                                         )

# # Save the result for use later on
# adata_result.write("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/multivelo/multivelo_result.h5ad")
