import scanpy as sc

# Subsets out each majorclass for downstream subclass annotation
# Ran on the reannotated object (using fetal reference)
adata = sc.read_h5ad('/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_ro_merged_annotated.h5ad')

list = ['AC','BC','HC','RGC','Cone','Rod','PRPC','NRPC']

for mclass in list:
    adata[adata.obs.majorclass == mclass].write(f'/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/subclass_annotation/chen_cherry_merged_{mclass}_raw.h5ad')