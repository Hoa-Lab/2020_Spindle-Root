import pandas as pd
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt

#---------------------------variable--------------------------
fmt='tif'

f_h5ad='./out/a09_meth-later_00_merge/harmony_merged.h5ad'
f_reg='./out/a09_meth-later_06_barplot-rss/all_reg.txt'
f_auc='./out/a09_meth-later_02_scenic/auc.csv'
fd_out='./out/a09_meth-later_07_reg-exp'

#------------------------setup--------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
adata=sc.read(f_h5ad)
df_auc=pd.read_csv(f_auc, index_col=0)

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


#########################################################################
#1. load reg list
l_reg=txt_to_gl(f_reg)

#2. add auc to adata
adata.obs=adata.obs.merge(df_auc, left_index=True, right_index=True)

##3. plot
#for reg in l_reg:
#	try:
#		f_out=f'{fd_out}/{reg}.{fmt}'
#		scanpy_exp(adata, reg, f_out)
#	except Exception as e:
#		print(str(e))
#		continue



#-------------------------replot---------------------------------
#3. replot
reg='Six1(+)'
vmax=0.25

f_out=f'{fd_out}/{reg}.{fmt}'
scanpy_exp(adata, reg, f_out, vmax=vmax)









