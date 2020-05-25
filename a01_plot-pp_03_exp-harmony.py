import pandas as pd
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt

#---------------------------variable--------------------------
fmt='png'

f_h5ad='./out/a00_preprocess_01_harmony/harmony_all.h5ad'
f_gl='./raw/ref/gl_stria.txt'
fd_out='./out/a01_plot-pp_03_exp-harmony'

#------------------------setup--------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
adata=sc.read(f_h5ad)

#------------------------function------------------------------
def txt_to_gl(fname):
	#read gene names in txt file
	txt=Path(fname).read_text()
	l_gl=txt.strip().split('\n')
	l_gene=[x.split(' ')[0] for x in l_gl]
	return l_gene

def scanpy_exp(adata, gene, f_out, vmin=-0.5, vmax=None, raw=True):
	ax=sc.pl.umap(adata, color=[gene], show=False, s=8, alpha=1, cmap='BuPu', use_raw=raw, vmin=vmin, vmax=vmax, frameon=False)
	ax.set_title(gene, fontsize=20, pad=20, weight='medium')
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return


##############################################################################
#---------------------plot----------------------------------------
#1. get gene list
l_gene=txt_to_gl(f_gl)
l_gene=['Otog']

#2. plot
for gene in l_gene:
	try:
		f_out=f'{fd_out}/{gene}.{fmt}'
		scanpy_exp(adata, gene, f_out)
	except Exception as e:
		print(str(e))
		continue


#-------------------------replot---------------------------------
##3. replot
#gene='Coch'
#vmax=2

#f_out=f'{fd_out}/{gene}.{fmt}'
#scanpy_exp(adata, gene, f_out, vmax=vmax)




