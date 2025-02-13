import scanpy as sc
import pandas as pd

#Plot the mature cells by major class

#Load in the data and the markers
adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad")
#adata = adata[adata.obs.source == 'chen']
#adata = adata[adata.obs.sampletype == 'organoid',]
adata = adata[adata.obs.sampletype == 'fetal',]

gene_list = {'RGC': ['OFF_MGC','ON_MGC','OFF_PGC','ON_PGC'],
			 'AC': ['GABAergic','Glycinergic','dual ACs','SACs'],
			 'HC': ['HC0','HC1'],
			 'BC': ['OFF-BC','ON-BC','RBC'],
			 'Cone': ['S_Cone','ML_Cone'],
			 'Rod': ['Rod']}

#Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#Plot by obs
for i in gene_list.keys():
	sc.pl.umap(adata, color='subclass', frameon=False, groups=gene_list[i], ncols = 4, size=5, save=f'_fet_mature_{i}.png')

