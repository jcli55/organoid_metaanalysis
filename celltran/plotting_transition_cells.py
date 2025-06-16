import scanpy as sc

# Some additional plotting of the transition cells after subclustering
# Jean Li - 04/16/25
# Add more functions as more plots are created

adata = sc.read_h5ad('/dfs3b/ruic20_lab/jeancl2/data/celltran/transition_cells_subclustered.h5ad')
#adata = sc.read_h5ad('/dfs3b/ruic20_lab/jeancl2/data/celltran/transition_cells_subclustered_scvi.h5ad')

# Reorder the ages so they will be in chronological order when plotting
adata.obs.age = adata.obs.age.cat.reorder_categories(['D35','D47','D75','D100','D123','D133','D183','D206','D243'])

# Plot transition index by organoid age
sc.pl.violin(
    adata,
    keys = ['transition_index'],
    groupby = 'age',
    jitter = 0.4,
    save = 'transition_index_by_age.png'
)

