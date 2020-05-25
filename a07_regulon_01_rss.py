import pandas as pd
import scanpy as sc
from pathlib import Path
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
from pyscenic.binarization import binarize

#-----------------------variable-----------------------------
fd_lbl='./out/a02_modularity_01_anno'
fd_auc='./out/a07_regulon_00_scenic'
fd_out='./out/a07_regulon_01_rss'

l_sample=['Ctrl', 'MethFix', 'RNAlater']

#---------------------setup----------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)


####################################################################
for sample in l_sample:
	#1. load
	df=pd.read_csv(f'{fd_auc}/auc_{sample}.csv', index_col=0)
	df_rss=df.copy()
	adata=sc.read(f'{fd_lbl}/{sample}.h5ad')
	
	#2. get anno
	df=df.merge(adata.obs.loc[:, ['anno']], left_index=True, right_index=True)
	l_anno=df['anno']
	
	#3. RSS
	df_rss=regulon_specificity_scores(df_rss, l_anno)
	df_rss=df_rss.T
	
	df_rss.to_csv(f'{fd_out}/rss_{sample}.csv')
	
