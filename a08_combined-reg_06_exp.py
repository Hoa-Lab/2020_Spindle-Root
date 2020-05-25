#plot top reg in combined RSS

import pandas as pd
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt

#---------------------------variable--------------------------
fmt='tif'
n=20

fd_h5ad='./out/a08_combined-reg_05_pp-exp'
f_rss='./out/a08_combined-reg_02_rss/rss_all.csv'
fd_out='./out/a08_combined-reg_06_exp'

#------------------------setup--------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
l_fname=list(Path(fd_h5ad).glob('*.h5ad'))
df_rss=pd.read_csv(f_rss, index_col=0)

#------------------------function------------------------------
def scanpy_exp(adata, gene, f_out, vmin=-0.1, vmax=None, raw=True):
	ax=sc.pl.umap(adata, color=[gene], show=False, s=8, alpha=1, cmap='BuPu', use_raw=raw, vmin=vmin, vmax=vmax)
	ax.set_title(gene, fontsize=20, pad=20, weight='medium')
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return
	
##############################################################################
#1. get reg list
df=df_rss.sort_values('Spindle-Root-1', ascending=False).iloc[0:n, :]
l_root=df.index.tolist()
df=df_rss.sort_values('Spindle-Root-2', ascending=False).iloc[0:n, :]
l_spin=df.index.tolist()
l_reg=list(set(l_root+l_spin))

#2. plot
for fname in l_fname:
	#1. setup
	sample=Path(fname).stem
	Path(f'{fd_out}/{sample}').mkdir(exist_ok=True)
	adata=sc.read(fname)
	
	#2. plot
	for reg in l_reg:
		try:
			f_out=f'{fd_out}/{sample}/{reg}.{fmt}'
			scanpy_exp(adata, reg, f_out)
		except Exception as e:
			print(str(e))
			continue
			
#-------------------------replot---------------------------------
##3. replot
#sample='Ctrl'
#reg='Coch'
#vmax=2

#adata=sc.read(f'{fd_h5ad}/{sample}.h5ad')
#f_out=f'{fd_out}/{sample}/{reg}.{fmt}'
#scanpy_exp(adata, reg, f_out, vmax=vmax)

