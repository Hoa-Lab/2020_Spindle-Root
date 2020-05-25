import pandas as pd
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt

#---------------------------variable--------------------------
fmt='png'

fd_h5ad='./out/a00_preprocess_00_pp'
f_gl='./raw/ref/gl_stria.txt'
fd_out='./out/a01_plot-pp_01_exp-single'

#------------------------setup--------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
l_fname=list(Path(fd_h5ad).glob('embed_*.h5ad'))

#------------------------function------------------------------
def txt_to_gl(fname):
	#read gene names in txt file
	txt=Path(fname).read_text()
	l_gl=txt.strip().split('\n')
	l_gene=[x.split(' ')[0] for x in l_gl]
	return l_gene

def scanpy_exp(adata, gene, f_out, vmin=-0.5, vmax=None, raw=True):
	ax=sc.pl.umap(adata, color=[gene], show=False, s=8, alpha=1, cmap='BuPu', use_raw=raw, vmin=vmin, vmax=vmax)
	ax.set_title(gene, fontsize=20, pad=20, weight='medium')
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return


##############################################################################
#---------------------plot----------------------------------------
#1. get gene list
l_gene=txt_to_gl(f_gl)

l_gene=['Tshz3']

#2. plot
for fname in l_fname:
	#1. setup
	sample=Path(fname).stem.split('_')[1]
	Path(f'{fd_out}/{sample}').mkdir(exist_ok=True)
	adata=sc.read(fname)
	
	#2. plot
	for gene in l_gene:
		try:
			f_out=f'{fd_out}/{sample}/{gene}.{fmt}'
			scanpy_exp(adata, gene, f_out)
		except Exception as e:
			print(str(e))
			continue
			
#-------------------------replot---------------------------------
##3. replot
#sample='Ctrl'
#gene='Coch'
#vmax=2

#adata=sc.read(f'{fd_h5ad}/embed_{sample}.h5ad')
#f_out=f'{fd_out}/{sample}/{gene}.{fmt}'
#scanpy_exp(adata, gene, f_out, vmax=vmax)




