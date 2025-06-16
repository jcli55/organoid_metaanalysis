import scanpy as sc
import pandas as pd

# Get DEGs between transition cells and non-transition cells
# Jean Li - 04/16/25

# Load in the data, normalize, find degs, write to a dataframe and save
adata = sc.read_h5ad("/dfs3b/ruic20_lab/jeancl2/data/chen_ro_clean.h5ad")
adata = adata[adata.obs.transition_cell == 'Transition Cell']

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Group three youngest ages as 'young' and rest as 'old' - late born cells significantly emerge after D100
adata.obs['age_class'] = adata.obs.age.replace(
    {
    'D35': 'young',
    'D47': 'young',
    'D75': 'young',
    'D100': 'old',
    'D123': 'old',
    'D133': 'old',
    'D183': 'old',
    'D206': 'old',
    'D243': 'old',
    }
)
adata.obs.age_class = adata.obs.age_class.astype('category')

sc.tl.rank_genes_groups(adata, 'age_class', pts=True)

for j in adata.obs.age_class.cat.categories:
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

    pd.DataFrame.to_csv(d, f'/dfs3b/ruic20_lab/jeancl2/data/celltran/{j}_transition_cell_deg.csv')