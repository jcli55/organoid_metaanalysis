import scanpy as sc
import numpy as np
import pandas as pd

# This plots a umap of cell types by age, leaving cells that do not fit in that age group gray
# Adapted the original for the reannotation using fetal as reference to use for transition cells
# 04/21/25 - Jean Li

adata = sc.read_h5ad("/dfs3b/ruic20_lab/jeancl2/data/celltran/velocity/multivelo_result.h5ad")

age_list = sorted(set(adata.obs.age.to_list()))
#for age in age_list:
#    adata.obs[age] = np.nan
#    for barcode in adata.obs_names:
#        if(adata.obs.loc[barcode, 'age'] == age):
#            adata.obs.loc[barcode, age] = adata.obs.loc[barcode, 'transition_index']
#    sc.pl.umap(adata, color=age, frameon=False, title=f'{age} Transition Index', size=50, na_in_legend=False, save=f'_transition_index_{age}.png')

for age in age_list:
    adata.obs[age] = np.nan
    for barcode in adata.obs_names:
        if(adata.obs.loc[barcode, 'age'] == age):
            adata.obs.loc[barcode, age] = adata.obs.loc[barcode, 'velo_s_norm_pseudotime']
    sc.pl.umap(adata, color=age, frameon=False, title=f'{age} Pseudotime', size=50, na_in_legend=False, save=f'_transition_pseudotime_{age}_mv.png')
