import pandas as pd
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt

#----------------------variable-----------------------------
fmt='png'

fd_h5ad='./out/a09_meth-later_00_merge'
f_cmap='./raw/ref/color_code.csv'
fd_out='./out/a09_meth-later_01_umap'

l_cell=['Marginal', 'Intermediate', 'Basal', 'Spindle-Root-1', 'Spindle-Root-2', 'Reissner', 'Fibrocyte', 'Macrophage', 'Unknown']

#-----------------------setup--------------------------------
Path(fd_out).mkdir(exist_ok=True)

#color map
df_cmap=pd.read_csv(f_cmap, index_col=0).reindex(l_cell)
dic_cmap=df_cmap.to_dict()['color']

#-----------------function-----------------------------------
def plot_anno(adata, f_out, l_cell=l_cell, title=None, dic_cmap=dic_cmap, fon=False, fs=(15,17), col='anno'):
	cmap=[dic_cmap[i] for i in l_cell]
	#2. plot
	plt.figure(figsize=fs)
	ax=sc.pl.umap(adata, color=[col], s=8, alpha=1, show=False, legend_loc='right margin', legend_fontsize=14, legend_fontweight='medium', frameon=fon, palette=cmap)
	ax.set_title(title, fontsize=16, pad=20, weight='medium')
	ax.xaxis.labelpad = 10
	ax.yaxis.labelpad = 10
	
	#2. save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return


#######################################################################
#1. load
ad_concat=sc.read(f'{fd_h5ad}/concat_merged.h5ad')
ad_harm=sc.read(f'{fd_h5ad}/harmony_merged.h5ad')

#2. plot
f_out=f'{fd_out}/concat.{fmt}'
plot_anno(ad_concat, f_out, title='MethFix-RNAlater (Without Integration)')

f_out=f'{fd_out}/harmony.{fmt}'
plot_anno(ad_harm, f_out, title='MethFix-RNAlater (With Integration)')


#3. plot sample
f_out=f'{fd_out}/concat_sample.{fmt}'
plot_anno(ad_concat, f_out, title='MethFix-RNAlater (Without Integration)', col='sample')

f_out=f'{fd_out}/harmony_sample.{fmt}'
plot_anno(ad_harm, f_out, title='MethFix-RNAlater (With Integration)', col='sample')












