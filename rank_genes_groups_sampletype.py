import scanpy as sc
import pandas as pd

# Script to run rank_genes_groups on a dataset and save the results to a csv (with percent expression)
# Usage log:
# Used to get degs between majorclass on the merged cherry and chen ro data
# Used to get percent expression of certain TFs for another experiment 11/27/23
# Used to get degs between majorclass on the reannotated organoid data (using fetal reference)
# Updated to get degs between majorclass by source
# Updated to get degs between fetal vs ro by majorclass
# Updated to split mature and precursor cells and get deg between fetal vs ro by majorclass

adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

for mclass in adata.obs.majorclass.cat.categories:
    temp = adata[adata.obs.majorclass==mclass]

    for maturity in temp.obs.maturityclass.cat.categories:
        temp2 = temp[temp.obs.maturityclass == maturity]

        sc.tl.rank_genes_groups(temp2, 'sampletype', pts=True)

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

            pd.DataFrame.to_csv(d, f'/storage/singlecell/jeanl/organoid/csv/reannotation/majorclass_degs/degs_by_sampletype/{j}_{mclass}_{maturity}_deg.csv')

        sc.pl.rank_genes_groups_dotplot(temp2, n_genes = 10, title=f'{mclass} {maturity}', save = f'{mclass}_{maturity}_deg_by_sampletype.png')