import scanpy as sc
import pandas as pd

# Get mean gex from fetal photoreceptors for candidate enhancers
# 03/17/2025 - Jean Li

# Generated object by subsetting cone and rod fetal cells from the merged data and concatenating
adata = sc.read_h5ad('/dfs3b/ruic20_lab/jeancl2/data/enhancer_assay_candidate/fetal_photoreceptors.h5ad')

enhancer_list = pd.read_csv('/dfs3b/ruic20_lab/jeancl2/data/enhancer_assay_candidate/enhancer_assay_candidate.csv', header=0)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

gex = []

for i in enhancer_list.Gene:
    gex.append(adata[:,i].X.mean())

enhancer_list['meanGEX'] = gex
enhancer_list.to_csv('/dfs3b/ruic20_lab/jeancl2/data/enhancer_assay_candidate/enhancer_gex_fetal.csv')