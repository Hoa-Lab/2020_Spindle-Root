import scanpy as sc
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

#--------------------variable----------------------------
fmt='png'

fd_in='./out/a02_modularity_01_anno'
fd_out='./out/a03_plot-modu_00_summary'

l_sample=['Ctrl', 'MethFix', 'RNAlater']
l_cell=['Marginal', 'Intermediate', 'Basal','Spindle-Root', 'Fibrocyte', 'Reissner', 'RBC', 'B Cell', 'Neutrophil', 'Macrophage']

#----------------------setup----------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)


#--------------------function---------------------------
def load_adata(fname):
	adata=sc.read(fname)
	adata.obs['cell']=adata.obs['anno'].astype('str')
	adata.obs['cell']=adata.obs['cell'].replace(['Spindle-Root-1', 'Spindle-Root-2'], ['Spindle-Root', 'Spindle-Root'])
	return adata

########################################################################
#-------------------heatmap------------------------------
#1. make df
adata1=load_adata(f'{fd_in}/Ctrl.h5ad')
adata2=load_adata(f'{fd_in}/MethFix.h5ad')
adata3=load_adata(f'{fd_in}/RNAlater.h5ad')

l_data=[]
for cell in l_cell:
	ctrl=int(cell in adata1.obs['cell'].tolist())
	meth=int(cell in adata2.obs['cell'].tolist())
	later=int(cell in adata3.obs['cell'].tolist())
	l_data.append((cell, ctrl, meth, later))
df=pd.DataFrame(l_data, columns=['cell']+l_sample)
df=df.set_index('cell')

df_labels=df.replace([0,1], ['No', 'Yes'])

#2. plot
#plt.rcParams['xtick.bottom']=plt.rcParams['xtick.labelbottom']=False
#plt.rcParams['xtick.top']=plt.rcParams['xtick.labeltop']=True
fig, ax=plt.subplots(figsize=(8,10))
ax=sns.heatmap(df, linewidths=2, cmap='icefire', annot=df_labels.values, fmt='', annot_kws={"size": 22}, cbar=False)

ax.set_title('Cells Identfied in Datasets', fontsize=24, pad=20, weight='semibold')
plt.xticks(fontsize=20, rotation=45)
plt.yticks(fontsize=24, rotation=0)
plt.ylabel('')
ax.xaxis.set_ticks_position('none') 
ax.yaxis.set_ticks_position('none') 
plt.tight_layout()
plt.savefig(f'{fd_out}/heatmap.{fmt}', dpi=300)
plt.close()











