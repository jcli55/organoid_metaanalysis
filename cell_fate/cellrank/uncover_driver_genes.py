import cellrank as cr
import scanpy as sc
import joblib

# Gets driver genes for cell type lineages - run after computing fate probabilities
# Written by Jean Li 02/13/24
# Adpated to run on hpc3 02/27/25

#adata = sc.read_h5ad('/storage/singlecell/jeanl/organoid/data/Chen/h5ad/multivelo/multivelo_result.h5ad') # Will need if plan to plot some embeddings - Note: can just use g.adata instead

g = joblib.load("/dfs3b/ruic20_lab/jeancl2/data/cellrank_full/gpcca_estimator_multivelo_scvi_states.h5ad")

driver_df = g.compute_lineage_drivers()
driver_df.to_csv("/dfs3b/ruic20_lab/jeancl2/data/cellrank_full/csv/lineage_drivers_scvi_states.csv")

# More plots and functionality in the multivelo tutorial on the website
