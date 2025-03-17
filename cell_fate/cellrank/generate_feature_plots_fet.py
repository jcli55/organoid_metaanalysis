import scanpy as sc
import pandas as pd

# Reads in an adata object and plots a list of genes after normalizing
# This one plots known lineage markers on the immature cells in Fetal

# Load in the data and the markers
adata = sc.read_h5ad("/dfs3b/ruic20_lab/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad")
adata = adata[adata.obs.sampletype == 'fetal',]
adata = adata[adata.obs.subclass.isin(['PRPC','NRPC','Cone Precursor','BC Precursor','AC Precursor','HC Precursor','RGC Precursor', 'Rod Precursor'])]

# Known lineage drivers for cell fate
gene_list={'rgc': ['POU4F2','ISL1','ATOH7','POU4F1'],
           'ac': ['NEUROD4','PAX6','PTF1A','PRDM13'],
           'hc': ['ONECUT1','ONECUT2','ONECUT3','PROX1'],
           'bc': ['PRDM1','VSX1','VSX2','OTX2'],
           'cone': ['PRDM1','CRX','THRB','OTX2'],
           'rod': ['NRL','CRX','NR2E3','OTX2']}

precursors = ['Cone Precursor','BC Precursor','AC Precursor','HC Precursor','RGC Precursor', 'Rod Precursor']
for subclass in precursors:
        sc.pl.umap(adata, color='subclass', groups=subclass, frameon=False, size=6, title='', legend_loc='none', save=f'_fet_{subclass}.png')

# Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Plot known lineage drivers for cell fate
for key in gene_list.keys():
        sc.pl.umap(adata, color=gene_list[key], frameon=False, ncols=4, size=6, save=f'_fet_{key}_genes.png')