import scanpy as sc
import matplotlib.pyplot as plt

adata = sc.read_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad')
adata = adata[adata.obs.sampletype=='organoid']
adata = adata[adata.obs.source=='chen']
adata.obs['age'] = adata.obs.age.cat.reorder_categories(['D35', 'D47', 'D75', 'D100', 'D123', 'D133', 'D183', 'D206', 'D243'])

sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

markers = ['ARR3','RHO']
df = sc.get.obs_df(adata, markers + ["age"])
grouped_df = df.groupby(["age"])
result = grouped_df.mean().reset_index()

# Set the 'Date' column as the index
result.set_index('age', inplace=True)

# Plot multiple lines using plt
plt.clf()
fig = plt.gcf()
fig.set_size_inches(8, 8)
#plt.plot(result.index, result['ALDH1A1'], marker='o', linestyle='-', label='ALDH1A1')
#plt.plot(result.index, result['ALDH1A3'], marker='^', linestyle='-', label='ALDH1A3')
#plt.plot(result.index, result['CYP26A1'], marker='s', linestyle='-', label='CYP26A1')
#plt.plot(result.index, result['FABP5'], marker='v', linestyle='-', label='FABP5')
#plt.plot(result.index, result['NFKB1'], marker='8', linestyle='-', label='NFKB1')
#plt.plot(result.index, result['VAX2'], marker='p', linestyle='-', label='VAX2')
#plt.plot(result.index, result['TBX5'], marker='*', linestyle='-', label='TBX5')

#plt.plot(result.index, result['THRB'], marker='o', linestyle='-', label='THRB')
#plt.plot(result.index, result['CRX'], marker='^', linestyle='-', label='CRX')
#plt.plot(result.index, result['PRDM1'], marker='s', linestyle='-', label='PRDM1')
#plt.plot(result.index, result['OTX2'], marker='v', linestyle='-', label='OTX2')

plt.plot(result.index, result['ARR3'], marker='o', linestyle='-', label='ARR3')
plt.plot(result.index, result['RHO'], marker='^', linestyle='-', label='RHO')

# Customize the plot (optional)
plt.title('')
plt.xlabel('')
plt.ylabel('log CPM')
plt.grid(True)
plt.legend()  # Add legend to distinguish lines
#plt.ylim(0, 5)
plt.tight_layout()
plt.savefig(
    "/storage/singlecell/jeanl/organoid/figures/plot5.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
