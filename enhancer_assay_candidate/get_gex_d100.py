import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('/dfs3b/ruic20_lab/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_ro_only.h5ad')
adata = adata[adata.obs.age == 'D100']

enhancer_list = pd.read_csv('/dfs3b/ruic20_lab/jeancl2/data/enhancer_assay_candidate/enhancer_assay_candidate.csv', header=0)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

gex = []

for i in enhancer_list.Gene:
    gex.append(adata[:,i].X.mean())

enhancer_list['meanGEX'] = gex
enhancer_list.to_csv('/dfs3b/ruic20_lab/jeancl2/data/enhancer_assay_candidate/enhancer_gex_d100.csv')