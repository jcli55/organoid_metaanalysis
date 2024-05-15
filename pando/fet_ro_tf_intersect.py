import pandas as pd
import numpy as np

# Gets the shared TFs between ro and fetal by major class
# by Jean Li 1/5/24 
# Updated for reannotation 3/12/24

files = ['AC_grn.csv', 'BC_grn.csv', 'HC_grn.csv', 'RGC_grn.csv', 'Cone_grn.csv', 'Rod_grn.csv', 'MG_grn.csv', 'PRPC_grn.csv', 'NRPC_grn.csv']
fet_files = ['AC_fet_grn.csv', 'BC_fet_grn.csv', 'HC_fet_grn.csv', 'RGC_fet_grn.csv', 'Cone_fet_grn.csv', 'Rod_fet_grn.csv', 'MG_fet_grn.csv', 'PRPC_fet_grn.csv', 'NRPC_fet_grn.csv']

for i in range(len(files)):
    grn = pd.read_csv(f'/storage/singlecell/jeanl/organoid/csv/reannotation/grns/{files[i]}')
    fet_grn = pd.read_csv(f'/storage/singlecell/jeanl/organoid/data/Chen/rds/fetal_grns/{fet_files[i]}')
    ro = np.array(grn.tf)
    fet = np.array(fet_grn.tf)
    same = np.intersect1d(ro, fet)
    pd.DataFrame(same.tolist()).to_csv(f'/storage/singlecell/jeanl/organoid/csv/reannotation/shared_tfs_grns/shared_{files[i]}')
