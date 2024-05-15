import scanpy as sc
import numpy as np
import pandas as pd

# This plots a umap of cell types by age, leaving cells that do not fit in that age group gray
# Did this for Chen organoid and Cherry organoid, and Chen fetal

adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/merged_chen/cherry_chen_merged_clean.h5ad")
fet = adata[adata.obs.sampletype == 'fetal']
ro = adata[adata.obs.sampletype == 'organoid']
ro_chen = ro[ro.obs.source == 'chen']
ro_cherry = ro[ro.obs.source == 'cherry']

sam_dict = {'fet': fet, 'ro_chen': ro_chen, 'ro_cherry': ro_cherry}

print("Finished loading and preparing data")

for sample in sam_dict.keys():
    print(f"Starting {sample}")
    sam = sam_dict[sample]
    age_list = sorted(set(sam.obs.age.to_list()))
    for age in age_list:
        print(f"Starting age {age}")
        sam.obs[age] = np.nan
        for barcode in sam.obs_names:
            if(sam.obs.loc[barcode, 'age'] == age):
                sam.obs.loc[barcode, age] = str(sam.obs.loc[barcode, 'subclass'])
            # else:
            #     sam.obs.loc[barcode, age] = np.nan
        sc.pl.umap(sam, color=age, frameon=False, title=f'{age}', size=5, na_in_legend=False, save=f'_{sample}_{age}.png')
