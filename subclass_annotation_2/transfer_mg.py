import scanpy as sc
import numpy as np
import pandas as pd

# Transferring just the MG labels to create a new majorclass column with MG
subtypes = sc.read_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/subclass_annotation/chen_cherry_merged_PRPC_mg_annotated.h5ad')
adata = sc.read_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad')

adata.obs['temp'] = np.nan
list = pd.Series(subtypes.obs_names).to_list()

print("Starting transfer...")
for barcode in adata.obs_names:
    if(barcode in list):
        adata.obs.loc[barcode, 'temp'] = str(subtypes.obs.loc[barcode, 'subclass'])
    else:
        adata.obs.loc[barcode, 'temp'] = str(adata.obs.loc[barcode, 'majorclass'])
print("Finished transfer")

print("Writing file...")
adata.write_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad')
