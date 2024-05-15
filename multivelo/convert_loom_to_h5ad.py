import scvelo as scv
import anndata

ldata_file = "/storage/chentemp/u250758/organoid_metaanalysis/Chen/merged.loom"
output_file = "/storage/singlecell/jeanl/organoid/data/Chen/h5ad/merged_spliced_unspliced.h5ad"

ldata = anndata.read_loom(ldata_file)

ldata.obs.index = [
    x.split(":", 1)[0] + "_" + x.split(":", 1)[1][:-1] + "-1-0" for x in ldata.obs.index
]
index = list(ldata.obs.index)

ldata.obs.index = index
ldata.raw = ldata

ldata.write(output_file)
