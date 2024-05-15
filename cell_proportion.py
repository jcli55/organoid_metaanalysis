import scanpy as sc
import pandas as pd

#Gets the number of cells of each cell type for each sampleid in an adata object and returns it as a dataframe in a csv
#Usage log:
#1. Used on majorclass, and subtype (had RGC6 etc)
#2. Updated on new subclass annotations I did for RGC specifically
#3. Updated to work on the reannotated object 2/21/2024
#4. Updated to work on the reannotated object with subclass labeled 2/29/2024
#5. Reran after updating subclass annotation (no mature RGC, and less MG) 3/18/2024
#Adata loading and processing
#adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/merged_chen/cherry_chen_merged_clean.h5ad")
adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad")

#Get the sampleids and cell types, initialize the dictionary
ids = sorted(set(adata.obs.sampleid.to_list()))
cell_types = sorted(set(adata.obs.subclass.to_list()))
dict = {}

#Get the number of cells for each cell type for each sampleid
for i in ids:
    subset = adata[adata.obs.sampleid == i,]
    sub_dict = {}
    for j in cell_types:
        sub_dict[j] = len(subset[subset.obs.subclass == j,])
    dict[i] = sub_dict

#Turn the dictionary into a csv
df = pd.DataFrame.from_dict(dict, orient="index") 
df.to_csv('/storage/singlecell/jeanl/organoid/csv/reannotation/reannotated_subclass_proportions.csv')
