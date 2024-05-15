import scanpy as sc
import numpy as np
import pandas as pd

print("Finished importing, starting reading...")
subtypes = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/merged_chen/subtype_annotation/cherry_chen_subtypes_annotated.h5ad")
subtypes.obs_names_make_unique()
adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/merged_chen/cherry_chen_subtypes_annotated.h5ad")
adata.obs_names_make_unique()
print("Finished reading files")

adata.obs['celltype'] = adata.obs['subtype']
adata.obs['subtype'] = np.nan
list = pd.Series(subtypes.obs_names).to_list()

print("Starting transfer...")
for i in range(0, len(adata)):
    barcode = adata.obs_names[i]
    if(barcode in list):
        adata.obs.iloc[i, 9] = str(subtypes.obs.loc[barcode, 'subtype'])
    else:
        adata.obs.iloc[i, 9] = str(adata.obs.iloc[i, 6])
print("Finished transfer")

print("Writing file...")
adata.write_h5ad("/storage/singlecell/jeanl/organoid/data/merged_chen/cherry_chen_subtypes_annotated.h5ad")