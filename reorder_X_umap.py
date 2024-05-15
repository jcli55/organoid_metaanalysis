import scanpy as sc
import numpy as np

# This script logs how I transferred the umap embedding from the original data to the velocyto spliced/unspliced data
# This was for the multivelo trajectory analysis for the organoid project
# By Jean Li 12/5/2023

# Update: Added to the run_multivelo.py
# 12/11/2023

# Read in the data, first is the source of the umap, second is the velocyto data to transfer to
adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/adata_rna.h5ad")
adata_rna = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/adata_rna_normalized.h5ad")

# Create a dictionary of sampleid/barcodes mapped to their umap coordinate
dict = {}
for i in range(0, len(adata)):
    dict[adata.obs_names[i]] = adata.obsm['X_umap'][i]

# Recreate the X_umap stacked array using the order the velocyte data was in
array = dict[adata_rna.obs_names[0]]
for i in range(1, len(adata_rna)):
    array = np.vstack([array, dict[adata_rna.obs_names[i]]])

# Set the umap
adata_rna.obsm['X_umap'] = array
