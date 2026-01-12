import scanpy as sc
import pandas as pd

# Script to run rank_genes_groups on a dataset and save the results to a csv (with percent expression)
# Usage log:
# Used to get degs between majorclass on the merged cherry and chen ro data
# Used to get percent expression of certain TFs for another experiment 11/27/23
# Used to get degs between majorclass on the reannotated organoid data (using fetal reference)
# Updated to get degs between majorclass by source
# Updated to get degs between fetal vs ro by majorclass
# Updated to split mature and precursor cells and get deg between fetal vs ro between majorclass and maturity class
# Updated to run on hpc3 and on just chen data (no cherry data) 1/7/26

adata = sc.read_h5ad("/dfs3b/ruic20_lab/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad")
adata = adata[adata.obs.source=='chen'].copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

for mclass in adata.obs.majorclass.cat.categories:
    temp = adata[adata.obs.majorclass==mclass].copy()

    for maturity in temp.obs.maturityclass.cat.categories:
        temp2 = temp[temp.obs.maturityclass == maturity].copy()

        sc.tl.rank_genes_groups(temp2, 'sampletype', method='wilcoxon', pts=True)

        for j in temp2.obs.sampletype.cat.categories:
            d = pd.DataFrame() 
            for k in ['names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']:  
                d[k] = temp2.uns["rank_genes_groups"][k][j]  
            
            f = pd.DataFrame(temp2.uns['rank_genes_groups']['pts'][j]) 
            list = [] 
            for i in d['names']: 
                list.append(f.loc[i][0]) 
            d['pts'] = list 

            f = pd.DataFrame(temp2.uns['rank_genes_groups']['pts_rest'][j]) 
            list = [] 
            for i in d['names']: 
                list.append(f.loc[i][0]) 
            d['pts_rest'] = list 

            pd.DataFrame.to_csv(d, f'/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_{j}_{mclass}_{maturity.replace(" ","_")}_deg.csv')

        sc.pl.rank_genes_groups_dotplot(temp2, n_genes = 10, title=f'{maturity} DEGs between RO and Fetal Retina', save = f'chen_{maturity.replace(" ","_")}_deg_by_sampletype.png')
