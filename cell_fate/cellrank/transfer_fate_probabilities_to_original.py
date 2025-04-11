import scanpy as sc

# Transfers the fate probabilites calculated using Cellrank to the obs of the original chen ro clean adata
# Written by Jean Li 04/11/25

# Load the cellrank object with fate probabilities
adata = sc.read_h5ad('/dfs3b/ruic20_lab/jeancl2/data/cellrank_full/cellrank_scvi_states.h5ad')
# Load the clean ro object
clean = sc.read_h5ad('/dfs3b/ruic20_lab/jeancl2/data/chen_ro_clean.h5ad')

ac = []
bc = []
cone = []
hc = []
mg = []
rgc = []
rod = []

# Fill out lists with the probabilities
for barcode in clean.obs_names.to_list():
    ac.append(adata.obs.loc[barcode, 'ac_fwd_probs'])
    bc.append(adata.obs.loc[barcode, 'bc_fwd_probs'])
    cone.append(adata.obs.loc[barcode, 'cone_fwd_probs'])
    hc.append(adata.obs.loc[barcode, 'hc_fwd_probs'])
    mg.append(adata.obs.loc[barcode, 'mg_fwd_probs'])
    rgc.append(adata.obs.loc[barcode, 'rgc_precursor_fwd_probs'])
    rod.append(adata.obs.loc[barcode, 'rod_fwd_probs'])

# Save the lists into adata as metadata columns
clean.obs['ac_fwd_probs'] = ac
clean.obs['bc_fwd_probs'] = bc
clean.obs['cone_fwd_probs'] = cone
clean.obs['hc_fwd_probs'] = hc
clean.obs['mg_fwd_probs'] = mg
clean.obs['rgc_precursor_fwd_probs'] = rgc
clean.obs['rod_fwd_probs'] = rod

# Save data
clean.obs.to_csv('/dfs3b/ruic20_lab/jeancl2/data/chen_ro_clean_metadata.csv')
clean.write('/dfs3b/ruic20_lab/jeancl2/data/chen_ro_clean.h5ad')