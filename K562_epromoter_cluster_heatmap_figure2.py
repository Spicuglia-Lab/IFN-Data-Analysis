import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager


df = pd.read_csv('K562_InducedGene_InducedEpromoter_Cluster_final_heatmap_data_final_plot.real.tsv.csv',sep="\t",header=0,
                usecols=['cluster number', 'noninduced_epromoter_induced_gene',
                         'induced_epromoter_induced_gene', 'induced_epromoter_noninduced_gene',
                         'constitutive_epromoter'])

#df = pd.read_csv('K562_InducedGene_InducedEpromoter_Cluster_final_heatmap_data_final_plot3.tsv',sep="\t",header=0)
#print(df.shape);exit(0)
#only till cluster 5 rows.
df = df.head(20)
df.set_index('cluster number', inplace=True)
#df = df.sort_values(by=['induced_epromoter_induced_gene','induced_epromoter_noninduced_gene'], ascending=False)
df1 = df[['noninduced_epromoter_induced_gene']]
df2 = df[['induced_epromoter_induced_gene']]
df3 = df[['induced_epromoter_noninduced_gene']]
df4 = df[['constitutive_epromoter']]

fig, (ax1,ax2,ax3,ax4) = plt.subplots(ncols=4)
fig.set_size_inches(3.5, 4.5)
#fig.set_size_inches(3.5, 11.5)
fig.subplots_adjust(wspace=0.0)

sns.heatmap(df1, cmap="Blues",annot=True,fmt='g', annot_kws={'size':10},ax=ax1, cbar=False, linewidths=.5,robust=True, vmin=0, vmax=5)
ax1.set_xticklabels(ax1.get_xticklabels(),rotation=30,horizontalalignment='right')
ax1.set_ylabel('')

sns.heatmap(df2, cmap="Greens",annot=True,fmt='g', annot_kws={'size':10},ax=ax2, yticklabels=[], cbar=False, linewidths=.5,robust=True, vmin=0, vmax=5)
ax2.set_xticklabels(ax2.get_xticklabels(),rotation=30,horizontalalignment='right')
ax2.set_ylabel('')

sns.heatmap(df3, cmap="Oranges",annot=True,fmt='g', annot_kws={'size':10},ax=ax3, yticklabels=[], cbar=False, linewidths=.5,robust=True, vmin=0, vmax=5)
ax3.set_xticklabels(ax3.get_xticklabels(),rotation=30,horizontalalignment='right')
ax3.set_ylabel('')

sns.heatmap(df4, cmap="Greys", annot=False, ax=ax4, yticklabels=[], cbar=False, linewidths=.5,robust=True, vmin=0, vmax=150)

for i, c in enumerate(df4.columns):
    for j, v in enumerate(df4[c]):
        if v > 0:
            ax4.text(i + 0.5, j + 0.5, '★', color='black', size=10, ha='center', va='center')
ax4.set_xticklabels(ax4.get_xticklabels(),rotation=30,horizontalalignment='right')
ax4.set_ylabel('')

#plt.show();exit(0)
#fig.savefig("gene_category_heatmap_all.png", bbox_inches='tight', pad_inches=1,dpi=150)
fig.savefig("gene_category_heatmap.png", bbox_inches='tight', pad_inches=0.5,dpi=150)
