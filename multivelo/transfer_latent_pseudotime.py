import scanpy as sc

# Transfers the latent and pseudotime calculated using Multivelo to the obs of the original chen ro clean adata
# Written by Jean Li 04/11/25

# Load the multivelo result object
adata = sc.read_h5ad('/dfs3b/ruic20_lab/jeancl2/data/multivelo_full/multivelo_result.h5ad')
# Load the clean ro object
clean = sc.read_h5ad('/dfs3b/ruic20_lab/jeancl2/data/chen_ro_clean.h5ad')

latent = []
pseudo = []

# Fill out lists with the probabilities
for barcode in clean.obs_names.to_list():
    latent.append(adata.obs.loc[barcode, 'latent_time'])
    pseudo.append(adata.obs.loc[barcode, 'velo_s_norm_pseudotime'])

# Save the lists into adata as metadata columns
clean.obs['latent_time'] = latent
clean.obs['velo_s_norm_pseudotime'] = pseudo


# Save data
clean.obs.to_csv('/dfs3b/ruic20_lab/jeancl2/data/chen_ro_clean_metadata.csv')
clean.write('/dfs3b/ruic20_lab/jeancl2/data/chen_ro_clean.h5ad')