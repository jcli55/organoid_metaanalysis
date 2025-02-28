import joblib
import cellrank as cr
import scanpy as sc

# CellRank meets RNA velocity
# Written by Jean Li 02/13/24
# Adpated to run on hpc3 02/26/25

adata = scv.read('/dfs3b/ruic20_lab/jeancl2/data/multivelo_full/multivelo_result.h5ad')

vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()

joblib.dump(vk, "/dfs3b/ruic20_lab/jeancl2/data/cellrank_full/velocity_kernel_multivelo.h5ad")
