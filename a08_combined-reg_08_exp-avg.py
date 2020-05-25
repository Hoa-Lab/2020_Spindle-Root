#plot top reg in combined RSS

import pandas as pd
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt

#---------------------------variable--------------------------
fmt='tif'

fd_h5ad='./out/a08_combined-reg_05_pp-exp'
f_reg='./out/a08_combined-reg_07_barplot-rss-avg/reg.txt'
fd_out='./out/a08_combined-reg_08_exp-avg'

#------------------------setup--------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
l_fname=list(Path(fd_h5ad).glob('*.h5ad'))

#------------------------function------------------------------
def txt_to_gl(fname):
	#read gene names in txt file
	txt=Path(fname).read_text()
	l_gl=txt.strip().split('\n')
	l_gene=[x.split(' ')[0] for x in l_gl]
	return l_gene
	
def scanpy_exp(adata, gene, f_out, vmin=-0.1, vmax=None, raw=True):
	ax=sc.pl.umap(adata, color=[gene], show=False, s=8, alpha=1, cmap='BuPu', use_raw=raw, vmin=vmin, vmax=vmax)
	ax.set_title(gene, fontsize=20, pad=20, weight='medium')
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return
	
##############################################################################
##1. get reg list
#l_reg=txt_to_gl(f_reg)

##2. plot
#for fname in l_fname:
#	#1. setup
#	sample=Path(fname).stem
#	Path(f'{fd_out}/{sample}').mkdir(exist_ok=True)
#	adata=sc.read(fname)
#	
#	#2. plot
#	for reg in l_reg:
#		try:
#			f_out=f'{fd_out}/{sample}/{reg}.{fmt}'
#			scanpy_exp(adata, reg, f_out)
#		except Exception as e:
#			print(str(e))
#			continue
			
#-------------------------replot---------------------------------
#3. replot
sample='harmony_all'
reg='Zfp160(+)'
vmax=0.3

adata=sc.read(f'{fd_h5ad}/{sample}.h5ad')
f_out=f'{fd_out}/{sample}/{reg}.{fmt}'
scanpy_exp(adata, reg, f_out, vmax=vmax)

