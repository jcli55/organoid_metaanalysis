import scanpy as sc
import pandas as pd

#Plot the mature cells by major class

#Load in the data and the markers
adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad")
#adata = adata[adata.obs.source == 'chen']
adata = adata[adata.obs.sampletype == 'organoid',]
#adata = adata[adata.obs.sampletype == 'fetal',]

gene_list = {'RGC': ['RGC Precursor'],
			 'AC': ['AC Precursor'],
			 'HC': ['HC Precursor'],
			 'BC': ['BC Precursor'],
			 'Cone': ['Cone Precursor'],
			 'Rod': ['Rod Precursor']}

#Plot by obs
for i in gene_list.keys():
	sc.pl.umap(adata, color='subclass', frameon=False, title = '', na_in_legend=False, groups=gene_list[i], ncols = 4, size=5, save=f'_{i}_progenitors.png')
