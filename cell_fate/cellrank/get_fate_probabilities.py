import cellrank as cr
import scanpy as sc
import joblib

# Computes fate probabilities for cell type lineages
# Written by Jean Li 02/13/24
# Adpated to run on hpc3 02/27/25

g = joblib.load("/dfs3b/ruic20_lab/jeancl2/data/cellrank_full/gpcca_estimator_multivelo_scvi_states.h5ad")

g.compute_fate_probabilities()
g.plot_fate_probabilities(same_plot=False, save='fate_probabilities_sep.png')
g.plot_fate_probabilities(same_plot=True, save='fate_probabilities.png')

joblib.dump(g, "/dfs3b/ruic20_lab/jeancl2/data/cellrank_full/gpcca_estimator_multivelo_scvi_states.h5ad")
