import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

#--------------------variable------------------------
fmt='png'

f_mtx='./out/a10_cellphonedb_00_clean/heatmap/count_network.txt'
fd_out='./out/a10_cellphonedb_01_hm-pair'

l_cell=['Marginal', 'Intermediate', 'Basal', 'Spindle', 'Root', 'Fibrocyte', 'Reissner', 'Macrophage']

#----------------------setup--------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)


##################################################################
#1. load
df=pd.read_csv(f_mtx, sep='\t')
df=pd.pivot_table(df, values='count', index=['SOURCE'], columns=['TARGET'])

#2. clean
df=df.reindex(l_cell).loc[:, l_cell]

#3. plot
fig, ax=plt.subplots(figsize=(13,11))
ax=sns.heatmap(df, cmap='Purples')

ax.xaxis.tick_top()
plt.xlabel('')
plt.ylabel('')
plt.xticks(fontsize=32, rotation=90, weight='medium')
plt.yticks(fontsize=32, rotation=0, weight='medium')

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=15)

plt.tight_layout()
plt.savefig(f'{fd_out}/heatmap.{fmt}', dpi=300)
plt.close()

