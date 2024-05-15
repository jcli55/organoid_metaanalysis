import scanpy as sc
import pandas as pd

list = ['ac','bc','hc','rgc','cone','rod']

for i in list:
	adata = sc.read_h5ad(f"/storage/singlecell/jeanl/organoid/data/merged_chen/subtype_annotation/cherry_chen_{i}_umap_sampleid.h5ad")
	adata.uns['sampletype_colors'] = ['#ff7f0e', '#1f77b4', '#2ca02c']
	list = sorted(set(pd.Series(adata[adata.obs.batch == 'Reference',].obs.author_cell_type).to_list()))
	'''
	sc.pl.embedding(
		adata,
		basis="X_umap",
		color=['sampletype'],
		frameon=False,
		size=3,
		title=f'Sample Type: {i}',
		save=f"_cherry_chen_{i}_sampletype.png"
	)
	'''
	sc.pl.embedding(
		adata,
		basis="X_umap",
		color=['author_cell_type'],
		frameon=False,
		size=3,
		legend_loc='on data',
		legend_fontoutline=1,
		legend_fontsize='x-small',
		groups=list,
		title=f'Cell Subtypes: {i}',
		save=f"_cherry_chen_{i}_subtype.png"
	)
