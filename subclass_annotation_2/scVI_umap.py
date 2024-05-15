import scanpy as sc
import scvi
import torch

#Adapted from Zhen's umap function
#Then further adapted to run on the reannotated object (w/ fetal reference)

list = ['AC','BC','HC','RGC','Cone','Rod']

for i in list:
	#Set params
	inputFile = f"/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/subclass_annotation/chen_cherry_merged_{i}_raw.h5ad"
	outputFile = f"/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/subclass_annotation/chen_cherry_merged_{i}_umap.h5ad"
	numFeatures = 2000
	batchKey = "sampleid"
	labelsKey = None
	embeddingList = ['batch', 'leiden', 'subclass','celltype']
	figureName = f"_cherry_chen_merged_{i}_leiden.png"

	adata = sc.read_h5ad(inputFile)
	temp = sc.read_h5ad(inputFile)

	sc.pp.highly_variable_genes(
		adata, flavor="seurat_v3", n_top_genes=numFeatures, subset=True
	)
	scvi.settings.seed = 0
	scvi.model.SCVI.setup_anndata(adata, batch_key=batchKey, labels_key=labelsKey)
	vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
	vae.train()
	adata.obsm["X_scVI"] = vae.get_latent_representation()

	sc.pp.neighbors(adata, use_rep="X_scVI")
	sc.tl.leiden(adata)
	sc.tl.umap(adata)

	#Transferring labels
	temp.obs["leiden"] = adata.obs.leiden
	temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
	temp.obsm["X_umap"] = adata.obsm["X_umap"]

	#Embedding
	sc.pl.umap(
	adata,
	color=embeddingList,
	frameon=False,
	ncols=2,
	legend_loc='on data',
	save=figureName
	)

	#Save h5ad
	temp.write_h5ad(outputFile)