import scanpy as sc
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

#--------------------variable----------------------------
fmt='tif'

fd_in='./out/a02_modularity_01_anno'
fd_out='./out/a03_plot-modu_01_pieplot'

l_sample=['Ctrl', 'MethFix', 'RNAlater']
l_cell=['Marginal', 'Intermediate', 'Basal','Spindle-Root'] 
cmap=['#F8766D', '#00BA38', '#619CFF', 'orange']

#----------------------setup----------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
df_cnt=pd.read_csv(f'{fd_out}/cnt.csv')

#--------------------function--------------------------
def load_adata(fname):
	adata=sc.read(fname)
	adata.obs['cell']=adata.obs['anno'].astype('str')
	adata.obs['cell']=adata.obs['cell'].replace(['Spindle-Root-1', 'Spindle-Root-2'], ['Spindle-Root', 'Spindle-Root'])
	return adata

def count_cells(adata, l_cell, col='cell'):
	#count cell number in adata
	l_data=[]
	for cell in l_cell:
		n=adata.obs.loc[adata.obs[col].str.contains(cell),:].shape[0]
		l_data.append((cell, n))
	df=pd.DataFrame(l_data, columns=['cell', 'count'])
	return df

def pieplot_sum(df, title, f_out, cmap=cmap, col_cnts='count', col_cells='cell', figsize=(10,10), labels=None):
	#labels=df[col_cells], if need cell label
	total=df[col_cnts].sum()
	fig, ax=plt.subplots(figsize=figsize)
	ax.pie(df[col_cnts], colors=cmap, labels=labels, textprops={'fontsize':45}, autopct='%1.1f%%')
	ax.set_title(title, fontsize=60, weight='medium')
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return
	
#######################################################################
##--------------------pie plot---------------------------
#for sample in l_sample:
#	#1. set var
#	adata=load_adata(f'{fd_in}/{sample}.h5ad')
#	
#	#2. count
#	df=count_cells(adata, l_cell)
#	
#	#3. pie plot
#	f_out=f'{fd_out}/pie_{sample}.{fmt}'
#	pieplot_sum(df, sample, f_out)
#	
#	#5. plot legend
#	fig, ax=plt.subplots(figsize=(10, 10))
#	sns.despine()
#	ax=sns.barplot(x='cell', y='count', data=df, hue='cell', palette=cmap)
#	plt.ylim([0, 1000])
#	ax.legend(frameon=False, prop={'size': 25})
#	plt.tight_layout()
#	plt.savefig(f'{fd_out}/legend.{fmt}', dpi=300)
#	plt.close()
	
#-------------pieplot on counting data--------------------
f_out=f'{fd_out}/pie_esimated.{fmt}'
title='Esitmated %'
pieplot_sum(df_cnt, title, f_out)


