import scanpy as sc
import pandas as pd

# Script to run rank_genes_groups on a dataset and save the results to a csv (with percent expression)
# Usage log:
# Used to get degs between majorclass on the merged cherry and chen ro data
# Used to get percent expression of certain TFs for another experiment 11/27/23
# Recreated (this version) to get percent expression of new stem cell markers across age 1/10/24

adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/adata_rna.h5ad")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.tl.rank_genes_groups(adata, 'age', pts=True)

age = sorted(set(adata.obs.age.to_list()))
marker_list = ['MITF','PMEL','TYRP1','TFPI2','SLC6A15','RAB38','GJA1','PCDH7','RELN','MECOM','CPAMD8','COL9A1','ZIC1','TGFB2']
 
for j in age:
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

    pd.DataFrame.to_csv(d, f'/storage/singlecell/jeanl/organoid/csv/sc_markers/chen_age_degs/{j}_deg.csv')
    
    # Filters out only the genes in the marker list
    custom_dict = {'MITF': 0,'PMEL': 1,'TYRP1': 2,'TFPI2': 3,'SLC6A15': 4,'RAB38': 5,'GJA1': 6,'PCDH7': 7,'RELN': 8,'MECOM': 9,'CPAMD8': 10,'COL9A1': 11,'ZIC1': 12,'TGFB2': 13}
    d = d[d['names'].isin(marker_list)]
    d = d.sort_values(by = ['names'], key = lambda x: x.map(custom_dict))
    pd.DataFrame.to_csv(d, f'/storage/singlecell/jeanl/organoid/csv/sc_markers/chen_age_degs/{j}_deg_filtered.csv')