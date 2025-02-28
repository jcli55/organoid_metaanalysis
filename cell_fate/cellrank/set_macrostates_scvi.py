import scanpy as sc
import joblib
import cellrank as cr

# Sets macrostates for cell fate prediction based on previously annotated data
# Written by Jean Li 02/13/24
# Adpated to run on hpc3 02/26/25

adata = sc.read('/dfs3b/ruic20_lab/jeancl2/data/multivelo_full/multivelo_result.h5ad')

vk = joblib.load("/dfs3b/ruic20_lab/jeancl2/data/cellrank_full/velocity_kernel_multivelo.h5ad")

ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
g = cr.estimators.GPCCA(combined_kernel)

g.compute_macrostates(n_states=14, cluster_key="maturityclass")
g.plot_macrostates(which="all", legend_loc="right", s=100, save='macrostates.png')

g.plot_macrostate_composition(key="maturityclass", figsize=(7, 4), save='macrostate_composition_maturityclass.png')
g.plot_macrostate_composition(key="age", figsize=(7, 4), save='macrostate_composition_age.png')

g.plot_coarse_T(annotate=False, save='coarse_t.png')

# After looking at above plots, determine initial states based on known biology
g.set_initial_states(states=['PRPC_3','PRPC_4'])
g.plot_macrostates(which="initial", legend_loc="right", s=100, save='macrostates_initial.png')

g.set_terminal_states(states = adata.obs.terminalstates)
g.plot_macrostates(which="terminal", legend_loc="right", s=100, save='macrostates_terminal.png')

g.plot_macrostates(which="terminal", discrete=False, save='macrostates_terminal_cont.png')

# write macrostates to AnnData
adata.obs["macrostates"] = g.macrostates
adata.uns["macrostates_colors"] = g.macrostates_memberships.colors

adata.write('/dfs3b/ruic20_lab/jeancl2/data/multivelo_full/multivelo_result.h5ad')

joblib.dump(g, "/dfs3b/ruic20_lab/jeancl2/data/cellrank_full/gpcca_estimator_multivelo_scvi_states.h5ad")