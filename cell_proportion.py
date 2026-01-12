import scanpy as sc
import pandas as pd

#Gets the number of cells of each cell type for each sampleid in an adata object and returns it as a dataframe in a csv
#Copied from organoid meta-analysis
#Adata loading and processing
adata = sc.read_h5ad("/dfs3b/ruic20_lab/jeancl2/data/crb1/crb1_merged_annotated.h5ad")

#Get the sampleids and cell types, initialize the dictionary
ids = sorted(set(adata.obs.sampleid.to_list()))
cell_types = sorted(set(adata.obs.majorclass.to_list()))
dict = {}

#Get the number of cells for each cell type for each sampleid
for i in ids:
    subset = adata[adata.obs.sampleid == i,]
    sub_dict = {}
    for j in cell_types:
        sub_dict[j] = len(subset[subset.obs.majorclass == j,])
    dict[i] = sub_dict

#Turn the dictionary into a csv
df = pd.DataFrame.from_dict(dict, orient="index") 
df.to_csv('/dfs3b/ruic20_lab/jeancl2/data/crb1/crb1_merged_annotated_proportions.csv')
