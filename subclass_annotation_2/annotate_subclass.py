import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad

# This script helps me annotate the subclasses by seeing which fetal subclasses got clustered to each leiden cluster and then I can create a mapping of cluster to subclass to use for organoid

adata = sc.read_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/subclass_annotation/chen_cherry_merged_{class}_umap.h5ad')

# Basic workflow for reclustering at various resolutions
'''
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.leiden(adata, resolution = 1)
sc.pl.umap(
    adata,
    color=['leiden'],
    frameon=False,
    ncols=2,
    legend_loc='on data',
    save='chen_cherry_{class}_leiden_{res}.png'
)
'''

# Get the subclass by leiden dataframe
df  = pd.DataFrame()
fet = adata[adata.obs.batch=='Reference']

for l in fet.obs.leiden.cat.categories:
    df = df.append(fet[fet.obs.leiden==l].obs['subclass'].value_counts(normalize=True, sort=False))
df.to_csv('/storage/singlecell/jeanl/organoid/csv/reannotation/fet_{class}_subclass_leiden.csv')

# I don't have a script for the actual mapping leiden to subclass, but the mappings are in my lab notebook February week of 3/1
# Here's a basic workflow for the mapping
# Addition: I went back and redid PRPC (maybe NRPC too) - check March week of 3/8
'''
cell_dict = {'subclass': ['leiden', 'clusters']} 

ro = adata[adata.obs.batch=='Query']
#fet = adata[adata.obs.batch=='Reference'] #run if didn't create the csv above

ro.obs['subclass'] = np.nan

for i in cell_dict.keys(): 
    ind = pd.Series(ro.obs.leiden).isin(cell_dict[i]) 
    ro.obs.loc[ind,'subclass'] = i 
ro.obs['subclass'] = ro.obs['subclass'].astype('category')
    
temp = ad.concat([ro, fet])

# This is both but I also made umaps for ro and fet only
sc.pl.umap(
    temp,
    color=['subclass'],
    frameon=False,
    ncols=2,
    save='chen_cherry_{class}_subclass.png'
)

temp.write_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/subclass_annotation/chen_cherry_merged_{class}_annotated.h5ad')

sc.pp.normalize_total(ro, target_sum=1e4)
sc.pp.log1p(ro)

sc.pl.dotplot(
    ro,
    var_names=[''],
    groupby='subclass',
    categories_order=[''],
    swap_axes=True,
    save='chen_cherry_{class}_makers.png'
)
'''