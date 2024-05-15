import scanpy as sc
import pandas as pd

#Reads in an adata object and plots a list of genes after normalizing
#Simple script where you just need to change the object path, the gene list, and potentially any subsetting

#Load in the data and the markers
adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad")
#adata = adata[adata.obs.source == 'chen']
#adata = adata[adata.obs.sampletype == 'organoid',]
adata = adata[adata.obs.sampletype == 'fetal',]

gene_list = {'RGC': ['POU4F2','ISL1','ATOH7','POU4F1'],
			 'AC': ['NEUROD4','PAX6','PTF1A','PRDM13'],
			 'HC': ['ONECUT1','ONECUT2','ONECUT3','PROX1'],
			 'BC': ['OTX2','VSX2','PRDM1','VSX1'],
			 'Cone': ['THRB','CRX','PRDM1','OTX2'],
			 'Rod': ['NRL','CRX','NR2E3','OTX2']}

#Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#Plot
for key in gene_list.keys():
	sc.pl.umap(adata, color=gene_list[key], frameon=False, ncols = 4, size=5, save=f'_{key}_fet_features.png')


'''
#Plot by obs
for age in adata.obs.age.cat.categories:
	data = adata[adata.obs.age==age]
	for gene in gene_list:
		sc.pl.umap(data, color=gene, frameon=False, title=f'{age}_{gene}', save=f'_ro_{gene}_{age}.png')
'''
