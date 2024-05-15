import scanpy as sc
import numpy as np

# Script to check some new stem cell markers from a paper from China
# by Jean Li 1/10/24

marker_list = ['MITF','PMEL','TYRP1','TFPI2','SLC6A15','RAB38','GJA1','PCDH7','RELN','MECOM','CPAMD8','COL9A1','ZIC1','TGFB2']

adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/merged_chen/cherry_chen_merged_clean.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.score_genes(adata, gene_list=marker_list)
fet = adata[adata.obs.sampletype == 'fetal']
ro = adata[adata.obs.sampletype == 'organoid']
ro_chen = ro[ro.obs.source == 'chen']
ro_cherry = ro[ro.obs.source == 'cherry']

# assign keys to the subset adata objects
sam_dict = {'fet': fet, 'ro_chen': ro_chen, 'ro_cherry': ro_cherry}

print("Finished loading and preparing data")

# Feature plots for combined data (fet + Chen ro + Cherry)
#for gene in marker_list:
#    sc.pl.umap(adata, color=gene, frameon=False, title=f'{gene}', size=5, save=f'_combined_{gene}.png')

for sample in sam_dict.keys():
    print(f"Starting {sample}")
    sam = sam_dict[sample]
    age_list = sorted(set(sam.obs.age.to_list()))
    for age in age_list:
        temp = sam[sam.obs.age == age]
        sc.pl.umap(temp, color='score', frameon=False, title=f'{age}', size=5, save=f'_{sample}_{age}_sc_marker_score.png')
        # print(f"Starting age {age}")
        # sam.obs[age] = np.nan
        # for barcode in sam.obs_names:
        #     if(sam.obs.loc[barcode, 'score'] != np.nan):
        #         sam.obs.loc[barcode, age] = sam.obs.loc[barcode, 'score']
        # sc.pl.umap(sam, color=age, frameon=False, title=f'{age}', size=5, na_in_legend=False, save=f'_{sample}_{age}_sc_marker_score.png')
