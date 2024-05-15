import pandas as pd

# This script gets the number of occurances of each tf found by Pando in either the organoids or fetal major classes
# by Jean Li 1/4/24 
# Updated to run on the reannotation 3/12/24

#files = ['AC_grn.csv', 'BC_grn.csv', 'HC_grn.csv', 'RGC_grn.csv', 'Cone_grn.csv', 'Rod_grn.csv', 'MG_grn.csv', 'PRPC_grn.csv', 'NRPC_grn.csv']
files = ['AC_fet_grn.csv', 'BC_fet_grn.csv', 'HC_fet_grn.csv', 'RGC_fet_grn.csv', 'Cone_fet_grn.csv', 'Rod_fet_grn.csv', 'MG_fet_grn.csv', 'PRPC_fet_grn.csv', 'NRPC_fet_grn.csv']
tf_list = []
count_dict = {}

for file in files:
    #grn = pd.read_csv(f'/storage/singlecell/jeanl/organoid/csv/reannotation/grns/{file}')
    grn = pd.read_csv(f'/storage/singlecell/jeanl/organoid/data/Chen/rds/fetal_grns/{file}')
    tfs = set(grn.tf.to_list())
    for tf in tfs:
        tf_list.append(tf)

for tf in set(tf_list):
    count_dict[tf] = tf_list.count(tf)

df = pd.DataFrame.from_dict(count_dict, orient='index')
#df.to_csv('/storage/singlecell/jeanl/organoid/csv/reannotation/grns/tf_occurances.csv')
df.to_csv('/storage/singlecell/jeanl/organoid/data/Chen/rds/fetal_grns/tf_occurances.csv')
