import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

df = pd.read_csv('heatmap/K562_InducedGene_InducedEpromoter_Cluster_final_heatmap_data.tsv',sep="\t",header=0)
df['clust'], df['genes'] = df['cluster number'].str.split(':', 1).str
df = df[['clust','noninduced_epromoter_induced_gene',
       'induced_epromoter_induced_gene', 'induced_epromoter_noninduced_gene',
       'constitutive_epromoter']]

#df.set_index('cluster number', inplace=True)
df.set_index('clust', inplace=True)
df = df.sort_values(by=['induced_epromoter_induced_gene', 'noninduced_epromoter_induced_gene'], ascending=False)

#print(df.head(4));exit(0)
df1 = df[['induced_epromoter_induced_gene']]
df2 = df[['induced_epromoter_noninduced_gene']]
df3 = df[['constitutive_epromoter']]
df4 = df[['noninduced_epromoter_induced_gene']]

fig, (ax1,ax2,ax3,ax4) = plt.subplots(ncols=4)
fig.set_size_inches(10.0, 18.27)
fig.subplots_adjust(wspace=0.0)

sns.heatmap(df1, cmap="Greens",annot=True,fmt='g', annot_kws={'size':10},ax=ax1, cbar=False, linewidths=.5,robust=True, vmin=0, vmax=5)
ax1.set_xticklabels(ax1.get_xticklabels(),rotation=30,horizontalalignment='right')
ax1.set_ylabel('')

sns.heatmap(df2, cmap="Oranges",annot=True,fmt='g', annot_kws={'size':10},ax=ax2, yticklabels=[], cbar=False, linewidths=.5,robust=True, vmin=0, vmax=5)
ax2.set_xticklabels(ax2.get_xticklabels(),rotation=30,horizontalalignment='right')
ax2.set_ylabel('')

#sns.heatmap(df3, cmap="Greys",annot=True,fmt='g', annot_kws={'size':10},ax=ax3, yticklabels=[], cbar=False, linewidths=.5,robust=True, vmin=0, vmax=5)
sns.heatmap(df3, cmap="Greys",annot=False,fmt='g', annot_kws={'size':10},ax=ax3, yticklabels=[], cbar=False, linewidths=.5,robust=True, vmin=0, vmax=5)
ax3.set_xticklabels(ax3.get_xticklabels(),rotation=30,horizontalalignment='right')
ax3.set_ylabel('')

sns.heatmap(df4, cmap="Purples",annot=True,fmt='g', annot_kws={'size':10},ax=ax4, yticklabels=[], cbar=False, linewidths=.5,robust=True, vmin=0, vmax=5)
ax4.set_xticklabels(ax4.get_xticklabels(),rotation=30,horizontalalignment='right')
ax4.set_ylabel('')

plt.show()

#fig.savefig("output.png")
