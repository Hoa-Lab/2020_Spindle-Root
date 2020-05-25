import pandas as pd
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt

#---------------------------variable--------------------------
fmt='tif'

fd_h5ad='./out/a00_preprocess_00_pp'
fd_gene='./out/a06_deg_04_barplot-avg-uniq'
fd_out='./out/a06_deg_06_exp'

l_sample=['Ctrl', 'MethFix', 'RNAlater']

#------------------------setup--------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#------------------------function------------------------------
def scanpy_exp(adata, gene, f_out, vmin=-0.5, vmax=None, raw=True):
	ax=sc.pl.umap(adata, color=[gene], show=False, s=8, alpha=1, cmap='BuPu', use_raw=raw, vmin=vmin, vmax=vmax)
	ax.set_title(gene, fontsize=20, pad=20, weight='medium')
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return


##################################################################
#1. get gene list
df1=pd.read_csv(f'{fd_gene}/spin1.csv', index_col=0)
df2=pd.read_csv(f'{fd_gene}/spin2.csv', index_col=0)
df=pd.concat([df1, df2])
l_gene=df.index.tolist()

#2. plot
for sample in l_sample:
	#1. setup
	Path(f'{fd_out}/{sample}').mkdir(exist_ok=True)
	adata=sc.read(f'{fd_h5ad}/embed_{sample}.h5ad')
	
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


