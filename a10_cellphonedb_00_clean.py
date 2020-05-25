#make two files:
#1. count: row is gene, column is cell id
#2. eta: columns: ['cell id', 'cell type']

import pandas as pd
import scanpy as sc
import numpy as np
from pathlib import Path

#-------------------------variable-----------------------
f_h5ad='./out/a09_meth-later_00_merge/clean_merged.h5ad'
fd_out='./out/a10_cellphonedb_00_clean'

l_cell=['Marginal', 'Intermediate', 'Basal', 'Root', 'Spindle','Fibrocyte', 'Reissner', 'Macrophage', 'Unknown']

#---------------------setup-----------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
adata=sc.read(f_h5ad)


#################################################################
#1. get df
df=pd.DataFrame(adata.X.toarray(), index=adata.obs.index, columns=adata.var.index).fillna(0)
df=df.astype('float')
df=df.T
df.index.name='gene'
df.index=df.index.str.upper()
df.to_csv(f'{fd_out}/cnt.txt', sep='\t')

#2. get meta
df=adata.obs.loc[:, ['anno']].copy()
df=df.fillna('Unknown')
df.columns=['cell_type']
df.index.name='Cell'

df['cell_type']=df['cell_type'].replace(['Spindle-Root-1', 'Spindle-Root-2'], ['Root', 'Spindle'])
#df['cell_type']=pd.Categorical(df['cell_type'], categories=l_cell, ordered=True)  #optional
df.to_csv(f'{fd_out}/meta.txt', sep='\t')

#################################################################
#run cellphonedb (gene name must be upper case for mouse data)
#https://github.com/Teichlab/cellphonedb

'''
cellphonedb method statistical_analysis meta.txt cnt.txt --counts-data=gene_name
cellphonedb plot heatmap_plot meta.txt --pvalues-path ./out/pvalues.txt --output-path ./heatmap
'''


