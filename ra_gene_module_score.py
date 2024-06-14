import scanpy as sc

# Adapted from Zhen's cell cycle score code
# Scores Retinoic Acid related genes

adata = sc.read(
	"/storage/singlecell/jeanl/organoid/data/merged_chen/reannotate_w_fetal/chen_cherry_merged_clean.h5ad"
)
adata = adata[adata.obs.sampletype == 'organoid',]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

regulating = ['ALDH1A1','ALDH1A3','CYP26A1']
tfs = ['FABP5','NFKB1','VAX2','TBX5']

sc.tl.score_genes(adata, gene_list=regulating, score_name = 'ra_regulating_genes_score')
sc.tl.score_genes(adata, gene_list=tfs, score_name = 'ra_tfs_score')

sc.pl.umap(
    adata,
    color=["ra_regulating_genes_score", "ra_tfs_score"],
    frameon=False,
    size=6,
    save='_ra_genes_score.png'
)