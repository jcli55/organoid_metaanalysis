import scanpy as sc
import scvi
import torch

# Subcluster the transition cells using batch correction from scVI

#Set params
numFeatures = 2000
batchKey = "sampleid"
labelsKey = None

#Read in data and subset to transition cells
adata = sc.read_h5ad('/dfs3b/ruic20_lab/jeancl2/data/chen_ro_clean.h5ad')
adata = adata[adata.obs.transition_cell == 'Transition Cell']

#Normalize and get HVGs
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata, n_top_genes=numFeatures, subset=True
)

#Setup and run scVI model
scvi.settings.seed = 0
scvi.model.SCVI.setup_anndata(adata, batch_key=batchKey, labels_key=labelsKey)
vae = scvi.model.SCVI(adata)
vae.train()
adata.obsm["X_scVI"] = vae.get_latent_representation()

#Calculate UMAP using scVI representation
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.leiden(adata)
sc.tl.umap(adata)

#Plot UMAP
sc.pl.umap(
    adata,
    color=['leiden', 'age', 'transition_index'],
    frameon=False,
    save='transition_cells_subclustered_scvi.png'
)
sc.pl.umap(
		adata,
		color=['ac_fwd_probs', 'bc_fwd_probs', 'cone_fwd_probs', 'hc_fwd_probs', 'mg_fwd_probs', 'rgc_precursor_fwd_probs', 'rod_fwd_probs'],
		frameon=False,
        ncols=4,
		save='transition_cells_subclustered_fates_scvi.png'
	)
sc.pl.umap(
		adata,
		color=['latent_time', 'velo_s_norm_pseudotime'],
		frameon=False,
		save='transition_cells_subclustered_pseudotime_scvi.png'
	)

#Save h5ad
adata.write_h5ad('/dfs3b/ruic20_lab/jeancl2/data/celltran/transition_cells_subclustered_scvi.h5ad')