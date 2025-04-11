import scanpy as sc
import numpy as np
import pandas as pd

# Transfer the transition index and label transition cells in the h5ad object obtained from CellTran
# Jean Li - 04/09/25

# Load in the data object and the celltran results
adata = sc.read_h5ad('/dfs3b/ruic20_lab/jeancl2/data/chen_ro_clean.h5ad')
transition_index = pd.read_csv('/dfs3b/ruic20_lab/jeancl2/data/celltran/transition_index_na_rm.csv', header=0, index_col = 0, names = ['index'])
transition_cells = pd.read_csv('/dfs3b/ruic20_lab/jeancl2/data/celltran/transition_cells.csv', header=0, index_col = 0, names = ['name'])

# Create new metadata columns
adata.obs['transition_index'] = np.nan
adata.obs['transition_cell'] = np.nan

# Get the list of cell names and iterate through them filling out with the results from celltran
list = pd.Series(adata.obs_names).to_list()
for i in range(0, len(adata)):
    barcode = adata.obs_names[i]
    # Transfer the transition indices
    if(barcode in transition_index.index):
        adata.obs.loc[barcode, 'transition_index'] = transition_index.loc[barcode, 'index']
    # Denote whether a cell was labeled as a transition cell or not
    if(barcode in transition_cells.name.to_list()):
        adata.obs.loc[barcode, 'transition_cell'] = 'Transition Cell'
adata.obs.transition_cell = adata.obs.transition_cell.astype('category')

# Save the updated object
adata.write_h5ad('/dfs3b/ruic20_lab/jeancl2/data/chen_ro_clean.h5ad')

# (Optional) Plot the transition indices and cells
sc.pl.umap(adata, color = ['transition_index','transition_cell'], frameon = False, na_in_legend = False, size = 6, save = 'transition_cells.png')

# (Optional) Plot the transition indices between transition and non-transition cells as violins
sc.pl.violin(
    adata[adata.obs.transition_cell == 'Transition Cell'],
    ["transition_index"],
    jitter=0.4,
    save='transition_index_of_transition_cells.png'
)
sc.pl.violin(
    adata[adata.obs.transition_cell != 'Transition Cell'],
    ["transition_index"],
    jitter=0.4,
    save='transition_index_non_transition_cells.png'
)

# (Optional) Save the transition indices of transition cells and non-transition cells
tcells = adata[adata.obs.transition_cell == 'Transition Cell']
tcells.obs.transition_index.to_csv('/dfs3b/ruic20_lab/jeancl2/data/celltran/transition_index_of_transition_cells.csv')

ntcells = adata[adata.obs.transition_cell != 'Transition Cell']
ntcells = ntcells[ntcells.obs.transition_index > 0]
ntcells.obs.transition_index.to_csv('/dfs3b/ruic20_lab/jeancl2/data/celltran/transition_index_non_transition_cells.csv')