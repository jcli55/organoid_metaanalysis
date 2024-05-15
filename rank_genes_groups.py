import scanpy as sc
import pandas as pd

# Script to run rank_genes_groups on a dataset and save the results to a csv (with percent expression)
# Usage log:
# Used to get degs between majorclass on the merged cherry and chen ro data
# Used to get percent expression of certain TFs for another experiment 11/27/23
# Used to get degs between majorclass on the reannotated organoid data (using fetal reference)

adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad")
adata = adata[adata.obs.sampletype=='organoid']

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.tl.rank_genes_groups(adata, 'majorclass', pts=True)
 
for j in adata.obs.majorclass.cat.categories:
    d = pd.DataFrame() 
    for k in ['names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']:  
        d[k] = adata.uns["rank_genes_groups"][k][j]  
    
    f = pd.DataFrame(adata.uns['rank_genes_groups']['pts'][j]) 
    list = [] 
    for i in d['names']: 
        list.append(f.loc[i][0]) 
    d['pts'] = list 

    f = pd.DataFrame(adata.uns['rank_genes_groups']['pts_rest'][j]) 
    list = [] 
    for i in d['names']: 
        list.append(f.loc[i][0]) 
    d['pts_rest'] = list 

    pd.DataFrame.to_csv(d, f'/storage/singlecell/jeanl/organoid/csv/reannotation/majorclass_degs/{j}_deg.csv')
