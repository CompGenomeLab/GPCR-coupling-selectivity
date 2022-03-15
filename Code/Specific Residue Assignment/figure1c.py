import scipy
import pandas as pd
binary_table=pd.read_table(r"D:\Users\suuser\Desktop\gitdocs\GPCR-ClassA\G protein coupling\binary_coupling_table.tsv", sep='\t', header=0)
# print(binary_table["GNAS"])
values_table=binary_table[["GNAS","GNAI1","GNAI2","GNAZ","GNAO1","GNAQ","GNA11","GNA14","GNA15","GNA12","GNA13"]]

import seaborn as sns; sns.set(color_codes=True,font_scale = 8)

Z_columns=scipy.cluster.hierarchy.linkage(values_table, method='complete', metric='euclidean', optimal_ordering=True)
Z_rows=scipy.cluster.hierarchy.linkage(values_table.T, method='complete', metric='euclidean', optimal_ordering=True)

g = sns.clustermap(values_table.T,cmap="#00441B",figsize=(80, 30),row_linkage=Z_rows,col_linkage=Z_columns,
                   vmax=1,vmin=0,xticklabels=binary_table["Residues"],linewidths=5,linecolor="black",
                   dendrogram_ratio=(.08, .2),cbar_pos=(0, 1, .03, .030))
for a in g.ax_row_dendrogram.collections:
    a.set_linewidth(10)

for a in g.ax_col_dendrogram.collections:
    a.set_linewidth(10)
#%%
g.savefig(r"D:\Users\suuser\Desktop\GNALL_all_residues_binary.svg")