import pandas as pd
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt

#---------------------------variable--------------------------
fmt='png'

f_h5ad='./out/a09_meth-later_00_merge/harmony_merged.h5ad'
f_info='./out/a09_meth-later_04_weight/weight.csv'
fd_out='./out/a09_meth-later_08_gene-exp'

l_reg=['Sall2(+)','Zfp160(+)', 'Mta3(+)', 'Ovol2(+)']
l_reg=['Rorb(+)']

#------------------------setup--------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
adata=sc.read(f_h5ad)
df_info=pd.read_csv(f_info)

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
##1. get gene list
l_gene=[]
for reg in l_reg:
	dfi=df_info.loc[df_info['regulon']==reg, :].copy()
	gene=dfi['gene'].tolist()
	l_gene.extend(gene)
l_gene=list(set(l_gene))

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
#sample='Ctrl'
#gene='Coch'
#vmax=2

#adata=sc.read(f'{fd_h5ad}/embed_{sample}.h5ad')
#f_out=f'{fd_out}/{sample}/{gene}.{fmt}'
#scanpy_exp(adata, gene, f_out, vmax=vmax)
