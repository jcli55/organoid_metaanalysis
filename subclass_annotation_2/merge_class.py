import scanpy as sc
import anndata as ad

# Merges the separate annotated classes together

list = ['AC','BC','HC','RGC','Cone','Rod','PRPC_mg','NRPC']

adatas = []
for mclass in list:
    adata = sc.read_h5ad(f'/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/subclass_annotation/chen_cherry_merged_{mclass}_annotated.h5ad')
    adatas.append(adata)

merged = ad.concat(adatas)

merged.write_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/subclass_annotation/chen_cherry_merged_subclass_annotated.h5ad')