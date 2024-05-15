import scanpy as sc
import pandas as pd

#Logs what I did to make Jonny's major class heatmap from his MERFISH paper on the organoid data
#Created 10/31/23 by Jean Li
#Adapted to run on reannotated data: 2/21/2024
#Updated to plot top 50 DEGs for each majorclass: 3/20/2024
#Updated to plot mclass markers again but with MG this time: 4/10/2024

#Set Params
gene_list = "/storage/singlecell/jeanl/organoid/csv/retina_majorclass_markers.csv"
figure_name = "_mclass_markers.png"

#gene_list = "/storage/singlecell/jeanl/organoid/csv/reannotation/majorclass_degs/top_50_deg.csv"
#figure_name = "_top_degs.png"

#Load in the data and the markers
adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad")
adata = adata[adata.obs.sampletype == 'organoid', ]
adata.obs['majorclass'] = adata.obs['majorclass'].cat.reorder_categories(new_categories = ['Rod','Cone','BC','AC','HC','RGC','MG','PRPC','NRPC'])

df = pd.read_csv(gene_list, header=None, sep=',')
#df = pd.read_csv(gene_list, header=None, sep='\t')

#Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#Create and fill the dict with markers
dict = {'Rod':[], 'Cone':[], 'BC':[], 'AC':[], 'HC':[], 'RGC':[], 'MG':[], 'PRPC':[], 'NRPC':[]}

for i in range(0,len(df)):
    dict[df.iloc[i,0]].append(df.iloc[i,1].upper())

#Plot
sc.pl.heatmap(adata, var_names=dict, groupby='majorclass', var_group_rotation=90, vmax=5, swap_axes=True, save=figure_name)
