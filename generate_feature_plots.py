import scanpy as sc
import pandas as pd

#Reads in an adata object and plots a list of genes after normalizing
#Simple script where you just need to change the object path, the gene list, and potentially any subsetting

#Load in the data and the markers
adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad")
adata = adata[adata.obs.sampletype == 'organoid',]
adata = adata[adata.obs.source == 'chen']

gene_list = ['FABP5','NFKB1','VAX2','TBX5']

#Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


#Plot
for gene in gene_list:
	sc.pl.umap(adata, color=gene, frameon=False, title=f'{gene}', save=f'_ro_{gene}.png')


#Plot by obs
for age in adata.obs.age.cat.categories:
	data = adata[adata.obs.age==age]
	for gene in gene_list:
		sc.pl.umap(data, color=gene, frameon=False, title=f'{age}_{gene}', save=f'_ro_{gene}_{age}.png')

