import scanpy as sc

#Generate a dotplot given a list of genes and a groups to categorize the data

#Load in the data and set the params
adata = sc.read_h5ad("/dfs3b/ruic20_lab/jeancl2/data/chen_ro_clean.h5ad")
gene_list = ['ISL1','NRL','CRX','NFIB','PRDM8','NFIA','HES6','RORA','PPARA','ONECUT3','ID2']
groupby = 'age'
categories_order = ['D35','D47','D75','D100','D123','D133','D183','D206','D243']
title = ''
filename = 'candidate_expression_by_age.png'

#Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pl.dotplot(adata, var_names=gene_list, groupby=groupby, categories_order=categories_order, title = title, save = filename)
