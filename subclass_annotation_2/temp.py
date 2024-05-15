import scanpy as sc
import scvi
import torch
import pandas as pd
import numpy as np
import anndata as ad

# Checking if annotated MG from PRPC lies close to the fetal mg

adata = sc.read_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/subclass_annotation/test_PRPC_mg.h5ad')
temp = sc.read_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/subclass_annotation/test_PRPC_mg.h5ad')

sc.pp.highly_variable_genes(
		adata, flavor="seurat_v3", n_top_genes=2000, subset=True
	)
scvi.settings.seed = 0
scvi.model.SCVI.setup_anndata(adata, batch_key='sampleid', labels_key=None)
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()
adata.obsm["X_scVI"] = vae.get_latent_representation()

sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.leiden(adata)
sc.tl.umap(adata)

#Transferring labels
temp.obs["leiden"] = adata.obs.leiden
temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
temp.obsm["X_umap"] = adata.obsm["X_umap"]

sc.pl.umap(
    adata,
    color=['batch', 'leiden', 'majorclass', 'subclass', 'age'],
    frameon=False,
    ncols=2,
    size=3,
    legend_loc='on data',
    save='test_PRPC_mg.png'
)

#Save h5ad
temp.write_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/subclass_annotation/test_PRPC_mg_umap.h5ad')

sc.pp.normalize_total(temp, target_sum=1e4)
sc.pp.log1p(temp)

sc.pl.umap(
    temp,
    color=['RLBP1', 'SLC1A3', 'CRYM','CLU'],
    frameon=False,
    ncols=2,
    save='test_PRPC_mg_markers.png'
)
'''
#-----------------------------After looking at the plots----------------------------------------------

# I also ended up loading the object back in interactively and plotting the same plots with just organoid cells, and also dotplot and age umap
# Assigning the cell types
# Decided organoid cluster 1 and 3 were MG
cell_dict = {'MG': ['1','3'], 'PRPC': ['0','2','4','5','6','7','8','9','10','11','12','13','14','15']} 
temp.obs['subclass'] = np.nan

# Only wanted to annotate clusters 1 and 3 for organoid; didn't want to mess with fetal annotation
ro = temp[temp.obs.batch=='Query']

for i in cell_dict.keys(): 
    ind = pd.Series(ro.obs.leiden).isin(cell_dict[i]) 
    ro.obs.loc[ind,'subclass'] = i 

fet = temp[temp.obs.batch=='Reference']
fet.obs['subclass'] = 'PRPC'

temp = ad.concat([fet, ro])
temp.obs['subclass'] = temp.obs['subclass'].astype('category')

temp.write_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/subclass_annotation/chen_cherry_merged_PRPC_mg_annotated.h5ad')
'''
