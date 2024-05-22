import scanpy as sc
import scvi
import torch
import anndata as ad

#This script documents how I integrated data from different papers for the meta analysis. (just the ones I couldn't get to)
#It can be adapted to integrate any number of datasets
#It largely follows the atlas integration of lung data tutorial from scvi.
#This version integrates the single-cell data (Clark, Reh, Roska) for the sc organoid atlas

file_name = ["Clark/clark_ro_rna_merged_raw_sampleid.h5ad", "Reh/reh_ro_rna_merged_raw_sampleid.h5ad", "Roska/roska_ro_rna_merged_raw_sampleid.h5ad"]
inPath = "/storage/singlecell/jeanl/organoid/data/"
outPath = "/storage/singlecell/jeanl/organoid/data/other_data_merged/scAtlas/"
outFile = "sc_ro_rna_merged_integrated.h5ad"
numFeatures = 2000
batch_key = "source"

#Merge the samples
data_list = []
for i in file_name:
    data = sc.read_h5ad(f'{inPath}{i}')
    data_list.append(data)
adata = ad.concat(data_list, join='inner')
adata.write_h5ad(f'{outPath}sc_ro_rna_merged_raw.h5ad')

#Load and preprocess the data
adata = sc.read_h5ad(f'{outPath}sc_ro_rna_merged_raw.h5ad')
temp = sc.read_h5ad(f'{outPath}sc_ro_rna_merged_raw.h5ad')

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=numFeatures,
    batch_key=batch_key,
    subset=True,
)

#Set up and train the model
scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key)
model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train()

#Get latent representation from the model
SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

#Calculate neighbors and assign leiden (can change resolution as see fit)
sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata)

#Visualize the data using UMAP
sc.tl.umap(adata)

#Transfer lables (comment out if just want the embedding w/o saving)
temp.obs['leiden'] = adata.obs['leiden']
temp.obsm[SCVI_LATENT_KEY] = adata.obsm[SCVI_LATENT_KEY]
temp.obsm['X_umap'] = adata.obsm['X_umap']
temp.write_h5ad(f'{outPath}{outFile}')

#Plot the embedding
sc.pl.umap(
    adata,
    color=[batch_key, "leiden"],
    frameon=False,
    ncols=3, save=f"_integrated.png"
)
