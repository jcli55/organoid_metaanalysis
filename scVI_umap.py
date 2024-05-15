import scanpy as sc
import scvi
import anndata as ad
import torch

#Adapted from Zhen's umap function

list = ['bc', 'hc', 'ac', 'rgc', 'cone', 'rod']

for i in list:
	#Set params
	inputFile = f"/storage/singlecell/jeanl/organoid/data/merged_chen/subtype_annotation/cherry_chen_{i}_raw.h5ad"
	outputFile = f"/storage/singlecell/jeanl/organoid/data/merged_chen/subtype_annotation/cherry_chen_{i}_umap_sampleid.h5ad"
	numFeatures = 2000
	batchKey = "sampleid"
	labelsKey = None
	embeddingList = ['batch', 'author_cell_type','sampletype','leiden']
	figureName = f"_cherry_chen_{i}.png"

	adata = sc.read_h5ad(inputFile)
	#temp = sc.read_h5ad(inputFile)

	#'''
	#Some processing to get sampleid to the reference
	adata_sub = adata[adata.obs.batch == "Reference", ]
	adata_sub.obs.sampleid = adata_sub.obs.sample_uuid
	adata_sub.obs.sampletype = "adult"
	adata = ad.concat([adata[adata.obs.batch == 'Query',], adata_sub])
	temp = adata
	adata.obs.sampleid
	#'''

	sc.pp.highly_variable_genes(
		adata, flavor="seurat_v3", n_top_genes=numFeatures, subset=True
	)
	adata
	temp
	scvi.settings.seed = 0
	scvi.model.SCVI.setup_anndata(adata, batch_key=batchKey, labels_key=labelsKey)
	vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
	vae.train()
	adata.obsm["X_scVI"] = vae.get_latent_representation()

	sc.pp.neighbors(adata, use_rep="X_scVI")
	sc.tl.leiden(adata)
	sc.tl.umap(adata)

	#adata.obsm["X_scVI_MDE"] = scvi.model.utils.mde(adata.obsm["X_scVI"])

	#Transferring labels
	temp.obs["leiden"] = adata.obs.leiden
	temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
	temp.obsm["X_umap"] = adata.obsm["X_umap"]
	#temp.obsm["X_scVI_MDE"] = adata.obsm["X_scVI_MDE"]

	'''
	#scANVI pipeline
	scanvi_model = scvi.model.SCANVI.from_scvi_model(
		vae,
		adata=adata,
		labels_key="majorclass",
		unlabeled_category="Unknown",
	)

	scanvi_model.train(max_epochs=100)

	adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()

	sc.pp.neighbors(adata, use_rep="X_scANVI")
	sc.tl.leiden(adata)
	adata.obsm["X_scANVI_MDE"] = scvi.model.utils.mde(adata.obsm["X_scANVI"])

	temp.obsm["X_scANVI"] = adata.obsm["X_scANVI"]
	temp.obsm["X_scANVI_MDE"] = adata.obsm["X_scANVI_MDE"]
	'''

	#Embedding
	sc.pl.embedding(
		adata,
		basis='X_umap',
		color=embeddingList,
		frameon=False,
		ncols=2,
		size=3,
		legend_loc='on data',
		save=figureName
	)

	#Save h5ad
	temp.write_h5ad(outputFile)
