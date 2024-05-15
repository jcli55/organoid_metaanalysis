import scanpy as sc
import scvi
import torch
import anndata as ad

#This script documents how I integrated data from different papers for the meta analysis. (just the ones I couldn't get to)
#It can be adapted to integrate any number of datasets
#It largely follows the atlas integration of lung data tutorial from scvi.

file_name = ["chen_fetal_ro_merged_raw"]
inPath = "/storage/singlecell/jeanl/organoid/data/Chen/"
outPath = "/storage/singlecell/jeanl/organoid/data/sanity_check/"
numFeatures = 2000
batch_key = "sampletype"

#Merge the samples
#adata = ad.concat([ , ], join='inner')

for i in file_name:
    adata = sc.read_h5ad(f"{inPath}{i}.h5ad")
    temp = sc.read_h5ad(f"{inPath}{i}.h5ad")

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

    #Visualize the data using this GPU-accelerated UMAP
    SCVI_MDE_KEY = "X_scVI_MDE"
    adata.obsm[SCVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCVI_LATENT_KEY])


    #Transfer lables (comment out if just want the embedding w/o saving)
    temp.obs['leiden'] = adata.obs['leiden']
    temp.obsm[SCVI_LATENT_KEY] = adata.obsm[SCVI_LATENT_KEY]
    temp.obsm[SCVI_MDE_KEY] = adata.obsm[SCVI_MDE_KEY]
    temp.write_h5ad(f"{outPath}{i}_integrated.h5ad")


    #Plot the embedding
    sc.pl.embedding(
        adata,
        basis=SCVI_MDE_KEY,
        color=[batch_key, "leiden"],
        frameon=False,
        ncols=3, save=f"_{i}.png"
    )

    #Sometimes the default cluster colors are ugly. If so use these:
    #('#000000', '#aaaa77', '#7777aa', '#ee8800', '#cccccc', '#004400',
    #       '#009900', '#00cc00', '#666600', '#449944', '#00ff00', '#ffff00',
    #       '#999900', '#cccc00', '#88cc88', '#ff3366', '#ff0000', '#cc3366',
    #       '#993366', '#cc0000', '#443366', '#ff4444', '#cc4444', '#994444',
    #       '#ff6666', '#cc6666', '#996666', '#ff8888', '#cc8888', '#cccccc',
    #       '#cccccc', '#998888', '#ccaaaa', '#ffcccc', '#ffeeee', '#000066',
    #       '#000099', '#0000ff', '#0000aa', '#aa00ff')

    #Need to make a list of color codes with length equal to the number of clusters and add to adata.uns as '(obs_name)_colors'
