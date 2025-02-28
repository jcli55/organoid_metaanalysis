import scanpy as sc
import cellrank as cr
import joblib

# Transfers the fate probabilites calculated using Cellrank to the obs of the adata object
# Written by Jean Li 02/07/25

# Load the object generated after get fate probabilities 
g = joblib.load('/dfs3b/ruic20_lab/jeancl2/data/cellrank_full/gpcca_estimator_multivelo_scvi_states.h5ad')

list = g.fate_probabilities.tolist()
ac = []
bc = []
cone = []
hc = []
mg = []
rgc = []
rod = []

for i in range(0,len(list)):
    ac.append(list[i][0])
    bc.append(list[i][1])
    cone.append(list[i][2])
    hc.append(list[i][3])
    mg.append(list[i][4])
    rgc.append(list[i][5])
    rod.append(list[i][6])

g.adata.obs['ac_fwd_probs'] = ac
g.adata.obs['bc_fwd_probs'] = bc
g.adata.obs['cone_fwd_probs'] = cone
g.adata.obs['hc_fwd_probs'] = hc
g.adata.obs['mg_fwd_probs'] = mg
g.adata.obs['rgc_precursor_fwd_probs'] = rgc
g.adata.obs['rod_fwd_probs'] = rod

g.adata.obs.to_csv('/dfs3b/ruic20_lab/jeancl2/data/cellrank_full/csv/cellrank_scvi_states.csv')
g.adata.write('/dfs3b/ruic20_lab/jeancl2/data/cellrank_full/cellrank_scvi_states.h5ad')
#joblib.dump(g, '')
