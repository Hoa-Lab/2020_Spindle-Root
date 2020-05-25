import pandas as pd
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt

#----------------------variable-----------------------------
fmt='png'

fd_merged='./out/a00_preprocess_01_harmony/'
fd_anno='./out/a02_modularity_01_anno'
f_cmap='./raw/ref/color_code.csv'
fd_out='./out/a02_modularity_02_anno-harmony'

l_sample=['Ctrl', 'MethFix', 'RNAlater']

l_cell=['Marginal', 'Intermediate', 'Basal', 'Spindle-Root-1', 'Spindle-Root-2', 'Reissner', 'Fibrocyte', 'Macrophage', 'B Cell', 'RBC', 'Neutrophil', 'Unknown']

#-----------------------setup--------------------------------
Path(fd_out).mkdir(exist_ok=True)
ad_concat=sc.read(f'{fd_merged}/concat_all.h5ad')
ad_harm=sc.read(f'{fd_merged}/harmony_all.h5ad')

#color map
df_cmap=pd.read_csv(f_cmap, index_col=0).reindex(l_cell)
dic_cmap=df_cmap.to_dict()['color']

#-----------------function-----------------------------------
def plot_anno(adata, f_out, title, dic_cmap=dic_cmap, l_cell=l_cell):
	cmap=[dic_cmap[i] for i in l_cell]
	#1. plot
	fig, ax=plt.subplots(figsize=(10,10))
	ax=sc.pl.umap(adata, color=['anno'], show=False, palette=cmap, frameon=False, s=8, alpha=0.8)
	#2. adjust
	ax.set_title(title, fontsize=16, pad=25, weight='medium')
	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, prop={'size': 13})

	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return


########################################################################
#1. make df_anno (clean index first)
df1=sc.read(f'{fd_anno}/Ctrl.h5ad').obs.loc[:, ['anno']]
df1.index=df1.index.map(lambda x: x+'-1')

df2=sc.read(f'{fd_anno}/MethFix.h5ad').obs.loc[:, ['anno']]
df2.index=df2.index.map(lambda x: x+'-0-0')

df3=sc.read(f'{fd_anno}/RNAlater.h5ad').obs.loc[:, ['anno']]
df3.index=df3.index.map(lambda x: x+'-1-0')

df_lbl=pd.concat([df1, df2, df3])

#2. add anno
ad_concat.obs=ad_concat.obs.merge(df_lbl, left_index=True, right_index=True, how='left')
ad_concat.obs['anno']=pd.Categorical(ad_concat.obs['anno'], categories=l_cell, ordered=True)

ad_harm.obs=ad_harm.obs.merge(df_lbl, left_index=True, right_index=True, how='left')
ad_harm.obs['anno']=pd.Categorical(ad_harm.obs['anno'], categories=l_cell, ordered=True)

#3. save
ad_concat.write(f'{fd_out}/concat_all.h5ad')
ad_harm.write(f'{fd_out}/harmony_all.h5ad')

#4. plot
f_out=f'{fd_out}/concat_all.{fmt}'
title='Combined Data Sets (Not Integrated)'
plot_anno(ad_concat, f_out, title)

f_out=f'{fd_out}/harmony_all.{fmt}'
title='Combined Data Sets (Integrated)'
plot_anno(ad_harm, f_out, title)






