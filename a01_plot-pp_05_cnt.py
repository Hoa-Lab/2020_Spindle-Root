#plot cnt bar plot

import scanpy as sc
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#--------------variables-------------------
fd_in='./out/a00_preprocess_00_pp'
fd_out='./out/a01_plot-pp_05_cnt'

cmap=['#73b4e6', '#e6a973', '#73e6e0','#e6d773', '#60db5a']
l_labels=['raw', 'filter_by_genes (>200)', 'filter_by_counts (<8000)',  'filter_by_Mt% (<10%)', 'de-doublet']
l_legends=['Raw', 'Gene Number >200', 'Total Counts <8000', 'Mitochondrial Genes (%) < 10', 'After De-doublet']

dic_ticklabels={'Ctrl': 'Ctrl', 'MethFix': 'MethFix', 'RNAlater': 'RNAlater'}

#---------------setup------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
df=pd.read_csv(f'{fd_in}/cell_cnts.csv')



#########################################################################
#------------------main-----------------------
#1. clean & normalize
l_pp=df['pp'].tolist()
df=df.set_index('pp')
l_samples=df.columns.tolist()
l_samples=[dic_ticklabels[i] for i in l_samples]

df=df.div(df.max())*100
df=df.stack().reset_index()
df.columns=['pp', 'sample', 'count']

#2. replace legend name
df['pp']=df['pp'].replace(l_labels, l_legends)

#3. setup plot
fig, ax=plt.subplots(figsize=(10,5))
sns.despine()
ax=sns.barplot(x='sample', y='count', hue='pp', data=df, palette=cmap)

#3. adjust plot
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=15)
ax.set_xlabel('')
ax.set_ylabel('Cell Percentage', fontsize=30, labelpad=8)
ax.set_xticklabels(l_samples, fontsize=20)
ax.tick_params(axis='x', which='major', pad=15)
plt.legend(loc=(1.04,0), fontsize=12, frameon=False)

#4. save plot
plt.tight_layout()
plt.savefig(f'{fd_out}/cell_cnts.png', dpi=300, bbox_inches="tight")
plt.close()

