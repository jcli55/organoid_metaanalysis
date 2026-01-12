import scanpy as sc
import pandas as pd

# Script to run rank_genes_groups on a dataset and save the results to a csv (with percent expression)
# Usage log:
# Used to get degs between majorclass on the merged cherry and chen ro data
# Used to get percent expression of certain TFs for another experiment 11/27/23
# Used to get degs between majorclass on the reannotated organoid data (using fetal reference)
# Updated to get degs between majorclass by source
# Updated to split mature and precursor cells and get deg by source
# Updated to get DEGs between NRPC and each precursor cell population 12/11/25

adata = sc.read_h5ad("/dfs3b/ruic20_lab/jeancl2/data/chen_ro_clean.h5ad")
adata = adata[adata.obs.maturityclass.isin(['NRPC','AC Precursor','BC Precursor','HC Precursor','RGC Precursor','Cone Precursor','Rod Precursor'])].copy()
adata.obs.maturityclass = adata.obs.maturityclass.cat.reorder_categories(['NRPC','AC Precursor','BC Precursor','HC Precursor','RGC Precursor','Cone Precursor','Rod Precursor'])

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

for mclass in adata.obs.maturityclass.cat.categories:
    if mclass == 'NRPC':
        continue
    else:
        temp = adata[adata.obs.maturityclass.isin(['NRPC', mclass])].copy()
        sc.tl.rank_genes_groups(temp, groupby='maturityclass', pts=True)
    
        d = pd.DataFrame() 
        for k in ['names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']:  
            d[k] = temp.uns["rank_genes_groups"][k][mclass]  
        
        f = pd.DataFrame(temp.uns['rank_genes_groups']['pts'][mclass]) 
        list = [] 
        for i in d['names']: 
            list.append(f.loc[i][0]) 
        d['pts'] = list 

        f = pd.DataFrame(temp.uns['rank_genes_groups']['pts_rest'][mclass]) 
        list = [] 
        for i in d['names']: 
            list.append(f.loc[i][0]) 
        d['pts_rest'] = list 

        pd.DataFrame.to_csv(d, f'/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/NRPC_precursor_deg/nrpc_{mclass}_deg.csv')
    
        sc.pl.rank_genes_groups_dotplot(temp, title=f'NRPC {mclass} DEG', n_genes=10, categories_order=['NRPC', mclass], save=f'nrpc_{mclass}_deg.png')