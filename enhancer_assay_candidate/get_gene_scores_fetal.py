import scanpy as sc
import pandas as pd

# Get mean gene score from fetal photoreceptors for candidate enhancers
# 04/01/2025 - Jean Li

# Gene scores from Zhen's object
adata = sc.read_h5ad('/dfs3b/ruic20_lab/singlecell/tingtiny/adult_dev_compare/adult_dev_compare/cellbygene/gene_score.h5ad')
adata = adata[adata.obs['sampleid'].isin(['Multi_Fetal_11w2d_FR', 'Multi_Fetal_11w2d_FR_2',
       'Multi_Fetal_11w2d_NR', 'Multi_Fetal_13W_FR', 'Multi_Fetal_13W_NR',
       'Multi_Fetal_14w5d_FR', 'Multi_Fetal_14w5d_NR', 'Multi_Fetal_19W4d_FR',
       'Multi_Fetal_19W4d_NR', 'Multi_Fetal_20W2d_FR', 'Multi_Fetal_20W2d_NR',
       'Multi_Fetal_23w1d_FR', 'Multi_Fetal_23w1d_NR', 'Multiome_10w_FR',
       'Multiome_10w_NR', 'Multiome_12w3d_FR', 'Multiome_12w3d_NR',
       'Multiome_14w2d_FR', 'Multiome_14w2d_NR', 'Multiome_16w4d_FR',
       'Multiome_16w4d_NR', 'Multiome_20w1d_FR', 'Multiome_20w1d_NR',
       'Multiome_23w4d_FR', 'Multiome_23w4d_NR'])]
adata = adata[adata.obs['majorclass'].isin(['Cone', 'Rod'])]  

enhancer_list = pd.read_csv('/dfs3b/ruic20_lab/jeancl2/data/enhancer_assay_candidate/enhancer_assay_candidate.csv', header=0)

# Normalize again after subsetting
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

genescore = []

for i in enhancer_list.Gene:
    genescore.append(adata[:,i].X.mean())

enhancer_list['meanGS'] = genescore
enhancer_list.to_csv('/dfs3b/ruic20_lab/jeancl2/data/enhancer_assay_candidate/enhancer_genescore_fetal.csv')