import scanpy as sc
import multivelo as mv

path = "/storage/singlecell/zz4/multi_organoid/data/"
samples = ("Multi_H9_D35",
            "Multi_H9_D47",
            "Multi_organoid_NRL_D075",
            "Multi_H9_D100",
            "Multi_NRL_GFP_D123",
           "Multi_organoid_D133",
           "Multi_Organoid_D183",
           "Multi_organoid_NRL_D206",
           "Multi_organoid_NRL_D243")
list = []

for sam in samples:
    temp = sc.read_10x_mtx(f'{path}{sam}/outs/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=True, gex_only=False)
    temp = temp[:,temp.var['feature_types'] == "Peaks"]
    temp = mv.aggregate_peaks_10x(temp, f'{path}{sam}/outs/atac_peak_annotation.tsv', f'{path}{sam}/outs/analysis/feature_linkage/feature_linkage.bedpe')

    names = []
    if (sam == "Multi_Organoid_D183"):
        for name in temp.obs_names:
                names.append(f'Multi_organoid_D183_{name}-0')
        temp.obs_names = names
    else:
        for name in temp.obs_names:
                    names.append(f'{sam}_{name}-0')
        temp.obs_names = names
    list.append(temp)

adata_atac = sc.concat(list)

adata_atac.write_h5ad("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/adata_atac.h5ad")
