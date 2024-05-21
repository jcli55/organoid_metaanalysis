import scanpy as sc
import anndata as ad

# Take the cellqc outputs for each Roska sample and add metadata and merge them

# Set the samples
samples = ['HOPAW6_1','HOPAW6_2','HOPBW6_1',
           'HOPBW6_2','HOPAW12_1','HOPAW12_2',
           'HOPBW12_1','HOPBW12_2','HOPAW18_1',
           'HOPAW18_2','HOPEx17W18_1','HOPEx17W18_2',
           'HOEx17W24_1','HOEx17W24_2','HOSEx13W24_1',
           'HOSEx13W24_2','HOEx13W30_1','HOEx13W30_2',
           'HOEx17W30_1','HOEx17W30_2','HOEx17W30_3',
           'HOEx17W30_4','M_w30_l1_1','M_w30_l1_2',
           'M_w30_l1_3','M_w30_l2_1','M_w30_l2_2',
           'M_w30_l2_3','HOEx17W38_a1','HOEx17W38_a2',
           'HOEx17W38_b1','HOEx17W38_b2','HOEx17W38_c1',
           'HOEx17W38_c2','M_w38_1_1','M_w38_1_2',
           'M_w38_1_3','M_w38_2_1','M_w38_2_2',
           'M_w38_2_3','S-191126-00103_1','S-191126-00103_2',
           'S-191126-00103_3','S-191126-00103_4']

# Add metadata
adata_list = []
for sam in samples:
    adata = sc.read_h5ad(f'/storage/chentemp/u250758/mrrdir/organoid_metaanalysis/Roska/cellqc_skipdk/result/{sam}.h5ad')
    age = 'NA'
    if sam in ['HOPAW6_1','HOPAW6_2','HOPBW6_1','HOPBW6_2']:
        age = 'wk6'
    elif sam in ['HOPAW12_1','HOPAW12_2','HOPBW12_1','HOPBW12_2']:
        age = 'wk12'
    elif sam in ['HOPAW18_1','HOPAW18_2','HOPEx17W18_1','HOPEx17W18_2']:
        age = 'wk18'
    elif sam in ['HOEx17W24_1','HOEx17W24_2','HOSEx13W24_1','HOSEx13W24_2']:
        age = 'wk24'
    elif sam in ['HOEx13W30_1','HOEx13W30_2','HOEx17W30_1','HOEx17W30_2','HOEx17W30_3','HOEx17W30_4','M_w30_l1_1','M_w30_l1_2','M_w30_l1_3','M_w30_l2_1','M_w30_l2_2','M_w30_l2_3']:
        age = 'wk30'
    elif sam in ['HOEx17W38_a1','HOEx17W38_a2','HOEx17W38_b1','HOEx17W38_b2','HOEx17W38_c1','HOEx17W38_c2','M_w38_1_1','M_w38_1_2','M_w38_1_3','M_w38_2_1','M_w38_2_2','M_w38_2_3']:
        age = 'wk38'
    elif sam in ['S-191126-00103_1','S-191126-00103_2','S-191126-00103_3','S-191126-00103_4']:
        age = 'wk46'

    adata.obs['age'] = age
    adata.obs['source'] = 'roska'
    adata.obs['sampleid'] = sam

    adata_list.append(adata)

# Merge the samples
merged = ad.concat(adata_list)
merged.write_h5ad('/storage/singlecell/jeanl/organoid/data/Roska/roska_ro_rna_merged_raw.h5ad')

# Update the obs_names
list = []
for i in range(0, len(merged)):
    list.append(f'{merged.obs.sampleid[i]}_{merged.obs_names[i]}')
merged.obs_names = list

merged.write_h5ad('/storage/singlecell/jeanl/organoid/data/Roska/roska_ro_rna_merged_raw_sampleid.h5ad')