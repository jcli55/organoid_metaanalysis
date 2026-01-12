import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Plot volcano plot based on output of Scanpy's rank_genes_groups()

filepaths = ['/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_AC_AC_deg.csv',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_AC_AC_Precursor_deg.csv',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_BC_BC_deg.csv',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_BC_BC_Precursor_deg.csv',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_Cone_Cone_deg.csv',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_Cone_Cone_Precursor_deg.csv',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_HC_HC_deg.csv',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_HC_HC_Precursor_deg.csv',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_MG_MG_deg.csv',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_NRPC_NRPC_deg.csv',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_PRPC_PRPC_deg.csv',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_RGC_RGC_Precursor_deg.csv',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_Rod_Rod_deg.csv',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_Rod_Rod_Precursor_deg.csv']

savepaths = ['/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_AC_AC_deg.png',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_AC_AC_Precursor_deg.png',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_BC_BC_deg.png',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_BC_BC_Precursor_deg.png',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_Cone_Cone_deg.png',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_Cone_Cone_Precursor_deg.png',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_HC_HC_deg.png',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_HC_HC_Precursor_deg.png',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_MG_MG_deg.png',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_NRPC_NRPC_deg.png',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_PRPC_PRPC_deg.png',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_RGC_RGC_Precursor_deg.png',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_Rod_Rod_deg.png',
            '/dfs3b/ruic20_lab/jeancl2/data/csv/reannotation/majorclass_degs/degs_by_sampletype/chen_organoid_Rod_Rod_Precursor_deg.png']

for i in range(0,len(filepaths)):
    filepath = filepaths[i]
    savepath = savepaths[i]

    # Sort DEGs by LFC and adjusted pvalue
    # Sometimes, rank_genes_groups() returns genes that are so significant the pval rounds to 0, so -log10 is undefined
    # Set the -log10 as -log10(1e-308) (the limit) for these values
    data = pd.read_csv(filepath)
    df = pd.DataFrame(data)
    log10_p = []
    for i in range(len(df)):
        if df['pvals_adj'][i] == 0:
            log10_p.append(-np.log10(1e-308))
        else:
            log10_p.append(-np.log10(df['pvals_adj'][i]))
    df['-log10_pvalue'] = log10_p
    colorlist=[]
    for i in range(len(df)):
        if df['logfoldchanges'][i] < -1 and df['-log10_pvalue'][i] > -np.log10(0.05):
            colorlist.append('red')
        elif df['logfoldchanges'][i] > 1 and df['-log10_pvalue'][i] > -np.log10(0.05):
            colorlist.append('green')
        else:
            colorlist.append('gray')
    df['color'] = colorlist
    
    plt.figure(figsize=(8, 6))
    plt.scatter(df['logfoldchanges'], df['-log10_pvalue'], s=10, alpha=0.7, color=df['color'], label=df['names'])
    plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--')
    plt.axvline(x=1, color='gray', linestyle='--')
    plt.axvline(x=-1, color='gray', linestyle='--')
    for i in range(len(df)):
        if df['-log10_pvalue'][i] > 150 and abs(df['logfoldchanges'][i]) > 5:
            plt.annotate(df['names'][i], (df['logfoldchanges'][i], df['-log10_pvalue'][i]),  # Point to annotate
                        xytext=(df['logfoldchanges'][i] + 0.05, df['-log10_pvalue'][i] + 0.05), # Position of the text
                        #arrowprops=dict(facecolor='gray', shrink=0.05), # Arrow properties
                        fontsize=5, color='black')
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-log10(p-value)')
    #plt.title('Volcano Plot')
    #plt.grid(True)
    plt.savefig(savepath)
