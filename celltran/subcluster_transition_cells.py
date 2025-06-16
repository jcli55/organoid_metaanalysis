import scanpy as sc

# Subcluster the transition cells

# Load in object after running label_transition_cells.py
adata = sc.read_h5ad('/dfs3b/ruic20_lab/jeancl2/data/chen_ro_clean.h5ad')
adata = adata[adata.obs.transition_cell == 'Transition Cell']

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
sc.pp.neighbors(adata, use_rep='X')
sc.tl.leiden(adata)
sc.tl.umap(adata)

sc.pl.umap(
		adata,
		color=['leiden', 'age', 'transition_index'],
		frameon=False,
		save='transition_cells_subclustered.png'
	)
sc.pl.umap(
		adata,
		color=['ac_fwd_probs', 'bc_fwd_probs', 'cone_fwd_probs', 'hc_fwd_probs', 'mg_fwd_probs', 'rgc_precursor_fwd_probs', 'rod_fwd_probs'],
		frameon=False,
        ncols=4,
		save='transition_cells_subclustered_fates.png'
	)
sc.pl.umap(
		adata,
		color=['latent_time', 'velo_s_norm_pseudotime'],
		frameon=False,
		save='transition_cells_subclustered_pseudotime.png'
	)

#Save h5ad
adata.write_h5ad('/dfs3b/ruic20_lab/jeancl2/data/celltran/transition_cells_subclustered.h5ad')