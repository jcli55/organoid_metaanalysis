import scanpy as sc
import pandas as pd

# Get DEGs between transition cells and non-transition cells
# Jean Li - 04/16/25

# Load in the data, normalize, find degs, write to a dataframe and save
adata = sc.read_h5ad("/dfs3b/ruic20_lab/jeancl2/data/chen_ro_clean.h5ad")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.tl.rank_genes_groups(adata, 'transition_cell', pts=True)

for j in adata.obs.transition_cell.cat.categories:
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

    pd.DataFrame.to_csv(d, f'/dfs3b/ruic20_lab/jeancl2/data/celltran/{j}_deg.csv')