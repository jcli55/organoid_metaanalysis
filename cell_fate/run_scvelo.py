import scanpy as sc
import scvelo as scv
import pandas as pd
import os
import scipy

# Read the processed snRNA-seq data with spliced/unspliced counts
adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/adata_rna_normalized_downsampled.h5ad")

# Add the MG annotation to the object
#metadata = pd.read_csv('/storage/singlecell/jeanl/organoid/csv/metadata_clean.csv', index_col = 0)
#list = []
#for barcode in adata.obs_names:
#    list.append(metadata.loc[barcode, 'majorclass'])
#adata.obs['class'] = list

#adata.write('/storage/singlecell/jeanl/organoid/data/Chen/h5ad/adata_rna_normalized_downsampled.h5ad')

scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode="dynamical")

scv.pp.neighbors(adata)
scv.tl.velocity_graph(adata, basis='umap', color='class', save='velocity_stream.png')

adata.write('/storage/singlecell/jeanl/organoid/data/Chen/h5ad/scvelo/scvelo_result.h5ad')
