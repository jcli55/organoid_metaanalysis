import scanpy as sc
import numpy as np
import pandas as pd

# How I am getting the subclass labels from each individual subset to the original big adata object
# Strat - load in each majorclass subset with subclass labeled, load in original object, merge all of the subsets together (also add NRPC too), go through loop to transfer labels
# I ended up writing a separate short script to merge the subsets together - merge_class.py
# Adapted to update the subclass annotation of the MG from PRPC
subtypes = sc.read_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/subclass_annotation/chen_cherry_merged_subclass_annotated.h5ad')
adata = sc.read_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad')

adata.obs['temp'] = np.nan
list = pd.Series(subtypes.obs_names).to_list()

print("Starting transfer...")
for barcode in adata.obs_names:
    if(barcode in list):
        adata.obs.loc[barcode, 'temp'] = str(subtypes.obs.loc[barcode, 'subclass'])
    else:
        adata.obs.loc[barcode, 'temp'] = str(adata.obs.loc[barcode, 'subclass'])
print("Finished transfer")

print("Writing file...")
adata.write_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad')
