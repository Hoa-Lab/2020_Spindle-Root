import pandas as pd
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt

#----------------------variable-----------------------------
fmt='tif'

fd_clus='./out/a05_cc_02_cluster'
f_cmap='./raw/ref/color_code.csv'
fd_out='./out/a05_cc_03_anno'

l_cell=['Marginal', 'Intermediate', 'Basal', 'Spindle-Root', 'Spindle-Root-1', 'Spindle-Root-2', 'Reissner', 'Fibrocyte', 'Macrophage', 'B Cell', 'RBC', 'Neutrophil', 'Unknown']

#-----------------------setup--------------------------------
Path(fd_out).mkdir(exist_ok=True)

#color map
df_cmap=pd.read_csv(f_cmap, index_col=0).reindex(l_cell)
dic_cmap=df_cmap.to_dict()['color']
dic_cmap['Spindle-Root']='#ebb931'

#-----------------function-----------------------------------
def anno_cell(adata, l_clus, l_anno, l_cell=l_cell):
	#1. edit l_cell
	l_c=[i for i in l_cell if i in l_anno]
	#2. anno
	adata.obs['anno']=adata.obs['leiden'].astype('int')
	adata.obs['anno']=adata.obs['anno'].replace(l_clus, l_anno)
	adata.obs['anno']=pd.Categorical(adata.obs['anno'], categories=l_c, ordered=True)
	return adata

def plot_anno(adata, l_anno, f_out, l_cell=l_cell, title=None, dic_cmap=dic_cmap, fon=False, fs=(10,10)):
	#1. edit l_cell
	l_c=[i for i in l_cell if i in l_anno]
	cmap=[dic_cmap[i] for i in l_c]
	#2. plot
	plt.figure(figsize=fs)
	ax=sc.pl.umap(adata, color=['anno'], s=8, alpha=1, show=False, legend_loc='right margin', legend_fontsize=10, legend_fontweight='medium', frameon=fon, palette=cmap)
	ax.set_title(title, fontsize=24, pad=20, weight='medium')
	ax.xaxis.labelpad = 10
	ax.yaxis.labelpad = 10
	#2. save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return


########################################################################
#------------------------Ctrl-------------------------------
sample='Ctrl'
n=20

l_anno=['Intermediate', 'Marginal', 'Intermediate', 'Marginal', 'Marginal', 'Marginal', 'RBC', 'Unknown', 'Basal', 'Marginal', 'B Cell', 'Spindle-Root', 'Unknown', 'RBC', 'Reissner', 'Neutrophil', 'Unknown', 'Marginal', 'Unknown', 'B Cell', 'Unknown']

#1. setup
l_clus=list(range(n+1))
adata=sc.read(f'{fd_clus}/{sample}.h5ad')

#2. annotation
anno_cell(adata, l_clus, l_anno)
adata.write(f'{fd_out}/{sample}.h5ad')

#3. plot
f_out=f'{fd_out}/{sample}.{fmt}'
plot_anno(adata, l_anno, f_out, title=sample)



#------------------------MethFix-------------------------------
sample='MethFix'
n=13

l_anno=['Intermediate', 'Marginal', 'Basal', 'Marginal', 'Fibrocyte', 'Reissner', 'Basal', 'Spindle-Root-2', 'Spindle-Root-1', 'Macrophage', 'Unknown', 'Unknown', 'Unknown', 'Unknown']

#1. setup
l_clus=list(range(n+1))
adata=sc.read(f'{fd_clus}/{sample}.h5ad')

#2. annotation
anno_cell(adata, l_clus, l_anno)
adata.write(f'{fd_out}/{sample}.h5ad')

#3. plot
f_out=f'{fd_out}/{sample}.{fmt}'
plot_anno(adata, l_anno, f_out, title=sample)


#------------------------RNAlater-------------------------------
sample='RNAlater'
n=12

l_anno=['Intermediate', 'Basal', 'Marginal', 'Marginal', 'Reissner', 'Basal', 'Fibrocyte', 'Spindle-Root-1', 'Spindle-Root-2', 'Unknown', 'Macrophage', 'Unknown', 'Unknown']

#1. setup
l_clus=list(range(n+1))
adata=sc.read(f'{fd_clus}/{sample}.h5ad')

#2. annotation
anno_cell(adata, l_clus, l_anno)
adata.write(f'{fd_out}/{sample}.h5ad')

#3. plot
f_out=f'{fd_out}/{sample}.{fmt}'
plot_anno(adata, l_anno, f_out, title=sample)


