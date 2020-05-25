import scanpy as sc
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import re

#--------------------variable------------------------
fmt='png'

fd_h5ad='./out/a00_preprocess_00_pp'
fd_out='./out/a01_plot-pp_04_batch-diff'

l_pair=['Ctrl-MethFix', 'Ctrl-RNAlater', 'MethFix-RNAlater']
l_sample=['Ctrl', 'MethFix', 'RNAlater']
dic_cmap={'Ctrl': '#4287f5', 'MethFix': '#f5a142', 'RNAlater': '#4bf542'}

nn=10
npc=30
n_genes=4000

#-----------------setup----------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#-------------function------------------------------
def plot_batch(adata, f_out, title, dic_cmap=dic_cmap):
	l_cell=[i for i in l_sample if i in adata.obs['sample'].unique()]
	cmap=[dic_cmap[i] for i in l_cell]
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
	
##################################################################
for pair in l_pair:
#	#1. setup
#	sample1=pair.split('-')[0]
#	sample2=pair.split('-')[1]
#	
#	ad1=sc.read(f'{fd_h5ad}/clean_{sample1}.h5ad')
#	ad2=sc.read(f'{fd_h5ad}/clean_{sample2}.h5ad')
#	
#	#2. get common gene
#	l_gene=[i for i in ad1.var.index if i in ad2.var.index]
#	l_gene=[i for i in ad1.var.index if i in ad2.var.index]
#	l_gene=[i for i in l_gene if not re.match('^mt-*', i)]
#	l_gene=[i for i in l_gene if not re.match('^[A-Z][A-Z]+', i)]
#	l_gene=[i for i in l_gene if not ('Rik' in i)]
#	l_gene=[i for i in l_gene if not ('-' in i)]
#	l_gene=[i for i in l_gene if len(i)>1]
#	l_gene=[i for i in l_gene if not re.match('^Gm\d+', i)]
#	l_gene.sort()
#	
#	#3. concat adata
#	ad1=ad1[:, l_gene].copy()
#	ad2=ad2[:, l_gene].copy()
#	adata=ad1.concatenate(ad2)
#	
#	#4. normalization, HVG
#	sc.pp.normalize_total(adata, exclude_highly_expressed=True)
#	sc.pp.log1p(adata)
#	adata.raw=adata
#	
#	sc.pp.highly_variable_genes(adata, n_top_genes=n_genes)
#	adata=adata[:, adata.var['highly_variable']].copy()

#	#5. PCA
#	sc.pp.regress_out(adata, ['n_counts', 'perc_others'])
#	sc.tl.pca(adata, svd_solver='arpack')

#	#6. calculate neighbor
#	sc.pp.neighbors(adata, n_neighbors=nn, n_pcs=npc)
#	
#	#7. embed- umap
#	sc.tl.umap(adata, n_components=2, random_state=42)
#	adata.write(f'{fd_out}/{pair}.h5ad')
	
	#8. plot
	adata=sc.read(f'{fd_out}/{pair}.h5ad')

	f_out=f'{fd_out}/{pair}.{fmt}'
	title=pair
	plot_batch(adata, f_out, title)
	
	















