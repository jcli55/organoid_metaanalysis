import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('/dfs3b/ruic20_lab/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_ro_only.h5ad')
adata = adata[adata.obs.source == 'Chen']

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

enhancer_list = pd.read_csv('/dfs3b/ruic20_lab/jeancl2/data/enhancer_assay_candidate/enhancer_assay_candidate.csv', header=0)
enhancer_list = enhancer_list.Gene.tolist()

D100 = adata[adata.obs.age == 'D100']
D183 = adata[adata.obs.age == 'D183']

sc.pl.dotplot(D100, var_names = enhancer_list, groupby = 'maturityclass', title = 'Candidate Enhancer Expression D100', save = 'candidate_enhancer_d100.png')
sc.pl.dotplot(D183, var_names = enhancer_list, groupby = 'maturityclass', title = 'Candidate Enhancer Expression D183', save = 'candidate_enhancer_d183.png')
sc.pl.dotplot(adata, var_names = enhancer_list, groupby = 'age', title = 'Candidate Enhancer Expression', save = 'candidate_enhancer.png')
