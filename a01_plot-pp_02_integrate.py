# plot sample label in integrated adata

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path

#---------------------variable-----------------------------
fmt='tif'

fd_in='./out/a00_preprocess_01_harmony'
fd_out='./out/a01_plot-pp_02_integrate'

l_sample=['Ctrl', 'MethFix', 'RNAlater']
dic_cmap={'Ctrl': '#4287f5', 'MethFix': '#f5a142', 'RNAlater': '#4bf542'}

#--------------------setup----------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#------------------function---------------------------------
def plot_batch(adata, f_out, title, dic_cmap=dic_cmap):
	cmap=[dic_cmap[i] for i in l_sample]
	#1. plot
	fig, ax=plt.subplots(figsize=(10,10))
	ax=sc.pl.umap(adata, color=['sample'], show=False, palette=cmap, frameon=False, s=8, alpha=0.8)
	#2. adjust
	ax.set_title(title, fontsize=16, pad=25, weight='medium')
	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, prop={'size': 13})

	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return



######################################################################
#-----------------concat data-------------------------------
fname='concat_all'
title='Combined Data Sets (Not Integrated)'

#1. load
adata=sc.read(f'{fd_in}/{fname}.h5ad')
adata.obs['sample']=pd.Categorical(adata.obs['sample'], categories=l_sample, ordered=True)

#2. plot
f_out=f'{fd_out}/{fname}.{fmt}'
plot_batch(adata, f_out, title)


#-----------------harmony data-------------------------------
fname='harmony_all'
title='Combined Data Sets (Integrated)'

#1. load
adata=sc.read(f'{fd_in}/{fname}.h5ad')
adata.obs['sample']=pd.Categorical(adata.obs['sample'], categories=l_sample, ordered=True)

#2. plot
f_out=f'{fd_out}/{fname}.{fmt}'
plot_batch(adata, f_out, title)











