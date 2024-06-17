import scanpy as sc
import pandas as pd

#Reads in an adata object and plots a list of genes after normalizing
#Simple script where you just need to change the object path, the gene list, and potentially any subsetting
# Used to plot RA related genes by age
# Used to plot Astro/microglia/RPE/Ciliary Margin genes - 5/20/24
# Used to plot majorclass markers on the scRNA data for the scAtlas - 5/28/24

#Load in the data and the markers
#adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad")
adata = sc.read_h5ad("q/storage/singlecell/jeanl/organoid/data/other_data_merged/scAtlas/sc_ro_rna_merged_annotated_ro_ref.h5ad")
adata = adata[adata.obs.batch == 'Query',]
#adata = adata[adata.obs.source == 'chen']

#gene_list = ['FABP5','NFKB1','VAX2','TBX5']
#gene_list = {'Astrocyte': ['GFAP','AQP4','S100B','GLUL'], 'Microglia': ['CD74','DOCK8','CTSB','SLC11A1'], 'RPE': ['RPE65','BEST1','LINC00276','COL8A1'], 'CiliaryMargin': ['ZIC1','WNT2B','KCNJ8','TPM2']}
gene_list = {'Rod': ['RHO','RCVRN','CRX'], 'Cone': ['PROM1','ARR3','OTX2'], 'BC': ['VSX1','VSX2'], 'AC': ['GAD1','GAD2','TFAP2A'], 'HC': ['TFAP2B','ONECUT1','ONECUT2'], 'RGC': ['NEFM','RBPMS','POU4F2'], 'MG': ['RLBP1','SLC1A3','SFRP2']}

#Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


#Plot
'''
for gene in gene_list:
	sc.pl.umap(adata, color=gene, frameon=False, title=f'{gene}', save=f'_ro_{gene}.png')


#Plot by obs
for age in adata.obs.age.cat.categories:
	data = adata[adata.obs.age==age]
	for gene in gene_list:
		sc.pl.umap(data, color=gene, frameon=False, title=f'{age}_{gene}', save=f'_ro_{gene}_{age}.png')
'''
for key in gene_list.keys():
	sc.pl.umap(adata, color=gene_list[key], frameon=False, ncols=3, size=6, save=f'_ro_{key}_genes.png')

#sc.pl.umap(adata, color=['source','majorclass'], frameon=False, ncols=2, size=6, save=f'_ro_source_class.png')
