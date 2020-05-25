import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import scrublet as scr
import anndata as ad

#---------------------------------variable----------------------------------
fd_in='./raw/count'
fd_out='./out/a00_preprocess_00_pp'

min_cells=3
min_genes=200
max_counts=8000
max_mt=0.1
nn=10
npc=30
n_genes=4000

l_sample=['Ctrl', 'MethFix', 'RNAlater']

#--------------------------------setup----------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#-------------------------------function-------------------------------
#1. functions
def filter_by_mt(adata, max_mt, l_prefix=('mt-', 'ERCC-', 'tg-', 'tg_')):
	other_genes = adata.var_names.str.startswith(l_prefix)
	adata.obs['perc_others'] = np.sum(adata[:, other_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
	adata=adata[adata.obs['perc_others'] < max_mt, :].copy()
	return adata

def dedoublets(adata, edr=0.1, npc=30, pctl=85, pl=False, f_out_fig=None, dpi=300):
	scrub=scr.Scrublet(adata.X, expected_doublet_rate=edr)
	doublet_scores, predicted_doublets=scrub.scrub_doublets(min_gene_variability_pctl=pctl, n_prin_comps=npc)
	#1. remove doublets
	adata.obs['doublets']=predicted_doublets
	adata=adata[adata.obs['doublets']==False,:].copy()
	#2. drop doublets column in obs
	adata.obs=adata.obs.drop('doublets', axis=1)
	#3. plot
	if pl:
		scrub.plot_histogram()
		plt.savefig(f_out_fig, dpi=dpi)
		plt.close()
	return adata


######################################################################################
l_number=[[], [], []]

#1. main loop
for i, sample in enumerate(l_sample):
	#1. load
	adata=sc.read_10x_mtx(f'{fd_in}/{sample}', var_names='gene_symbols', cache=False)
	adata.obs['sample']=sample
	l_number[i].append(adata.X.shape[0])
	
	#2. save raw adata
	adata.write(f'{fd_out}/raw_{sample}.h5ad')
	
	#3. filter genes by min cells
	sc.pp.filter_genes(adata, min_cells=min_cells)

	#4. filter cells by minimum genes
	sc.pp.filter_cells(adata, min_genes=min_genes)
	l_number[i].append(adata.X.shape[0])

	#5. filter cells by maximum counts
	sc.pp.filter_cells(adata, max_counts=max_counts)
	l_number[i].append(adata.X.shape[0])
	
	#6. filter by mt%
	adata=filter_by_mt(adata, max_mt)
	l_number[i].append(adata.X.shape[0])

	#8. de-doublet
	adata=dedoublets(adata)
	l_number[i].append(adata.X.shape[0])

	#9. save data
	adata.write(f'{fd_out}/clean_{sample}.h5ad')
	
	#10. normalization
	sc.pp.normalize_total(adata, exclude_highly_expressed=True)
	adata.write(f'{fd_out}/norm_{sample}.h5ad')
	
	#11. HVG
	sc.pp.log1p(adata)
	adata.raw=adata
	
	sc.pp.highly_variable_genes(adata, n_top_genes=n_genes)
	adata=adata[:, adata.var['highly_variable']].copy()
	
	#12. PCA
	sc.pp.regress_out(adata, ['n_counts', 'perc_others'])
	sc.tl.pca(adata, svd_solver='arpack')
	
	#13. calculate neighbor
	sc.pp.neighbors(adata, n_neighbors=nn, n_pcs=npc)
	
	#14. embed- umap
	sc.tl.umap(adata, n_components=2, random_state=42)
	sc.external.tl.phate(adata, n_components=3, random_state=42)
	
	#15. save
	adata.write(f'{fd_out}/embed_{sample}.h5ad')

#----------------------------------------------------------------------------
#2. create count df
l_pp=['raw', 'filter_by_genes (>200)', 'filter_by_counts (<8000)', 'filter_by_Mt% (<10%)', 'de-doublet']
dic_data={i:j for i, j in zip(l_sample, l_number)}

df=pd.DataFrame(dic_data)
df=df.loc[:, ['Ctrl', 'MethFix', 'RNAlater']]
df['pp']=l_pp
df.to_csv(f'{fd_out}/cell_cnts.csv', index=False)


