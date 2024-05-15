import scanpy as sc
import pandas as pd
import numpy as np

# Script logging how I made the csv's comparing the annotation of the old object to the new annotation (using fetal as reference)

adata = sc.read_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/cherry_chen_merged_clean.h5ad')
ro = sc.read_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_ro_merged_annotated.h5ad')

df = pd.DataFrame()
names = []
for name in ro[ro.obs.batch=='Query'].obs_names.to_list():
    #name = name.replace('-1-1-0', '-1-0')
    name = name.replace('-1-0-1', '-1-0')
df['cell_names'] = ro[ro.obs.batch=='Query'].obs_names
df['age'] = ro[ro.obs.batch=='Query'].obs.age.to_list()
df['new_class'] = ro[ro.obs.batch=='Query'].obs.majorclass.to_list()
old_class = []
for barcode in df.cell_names:
    try:
        old_class.append(adata.obs.loc[barcode, 'subclass'])
    except KeyError:
        old_class.append('NA')
df['old_class'] = old_class
df.to_csv('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/ro_annotation_compare.csv')

# This compares the fetal annotations, but I commented out because I needed to rerun ro annotation again (above)
#del df
#del old_class
#del names
#
#df = pd.DataFrame()
#names = []
#for name in adata[adata.obs.sampletype=='fetal'].obs_names.to_list():
#    name = name.replace('-1-0', '-1-1')
#    names.append(name)
#df = pd.DataFrame()
#df['cell_names'] = names
#df['age'] = adata[adata.obs.sampletype=='fetal'].obs.age.to_list()
#df['old_class'] = adata[adata.obs.sampletype=='fetal'].obs.subclass.to_list()
#new_class = []
#for barcode in df.cell_names:
#    try:
#        new_class.append(ro.obs.loc[barcode, 'majorclass'])
#    except KeyError:
#        new_class.append('NA')
#df['new_class'] = new_class
#df.to_csv('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/fet_annotation_compare.csv')
