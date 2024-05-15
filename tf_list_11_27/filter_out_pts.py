import pandas as pd

# Get the TF list
df = pd.read_csv('/storage/singlecell/jeanl/organoid/csv/tf_list.csv')
list = []
for genes in df['genes']:
    list.append(genes)

# Filter the deg dataframes to only have the tfs of interest
majorclass = ['AC', 'BC', 'HC', 'RGC', 'Rod', 'Cone', 'MG', 'PRPC', 'NRPC']
for mclass in majorclass:
    df = pd.read_csv(f'/storage/singlecell/jeanl/organoid/csv/chen_majorclass_progs_deg/{mclass}_deg.csv')
    df = df[df['names'].isin(list)]
    df = df.sort_values('names')
    df.to_csv(f'/storage/singlecell/jeanl/organoid/csv/chen_majorclass_progs_deg/{mclass}_deg_filtered.csv')
