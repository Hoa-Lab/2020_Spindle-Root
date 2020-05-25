#auc index maybe duplicated, extra code is required to split auc into different files

import pandas as pd
import scanpy as sc
from pathlib import Path
from pyscenic.rss import regulon_specificity_scores

#-----------------------variable-----------------------------
fd_lbl='./out/a02_modularity_01_anno'
fd_auc='./out/a08_combined-reg_01_scenic'
fd_out='./out/a08_combined-reg_02_rss'

l_sample=['Ctrl', 'MethFix', 'RNAlater']

#---------------------setup----------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
df_auc=pd.read_csv(f'{fd_auc}/auc.csv', index_col=0)


######################################################################
##1. get annotation list
#l_anno=[]
#for sample in l_sample:
#	adata=sc.read(f'{fd_lbl}/{sample}.h5ad')
#	df=adata.obs.loc[:, ['anno']].copy().fillna('Unknown')
#	l_anno.extend(df['anno'].tolist())
#l_anno=pd.Series(l_anno)

##2. RSS- all
#df_rss=regulon_specificity_scores(df_auc, l_anno)
#df_rss=df_rss.T
#df_rss.to_csv(f'{fd_out}/rss_all.csv')

##3. split auc  (this is extra step, since index maybe duplicate in different samples)
#for sample in l_sample:
#	#1. get shape
#	adata=sc.read(f'{fd_lbl}/{sample}.h5ad')
#	n=adata.obs.shape[0]
#	
#	#2. split auc
#	df=df_auc.iloc[0:n,:].copy()
#	df_auc=df_auc.iloc[n:, :].copy()
#	
#	#3. save
#	df.to_csv(f'{fd_out}/auc_{sample}.csv')
#	print(df.shape)

#4. RSS- each
for sample in l_sample:
	#1. load
	dfi_auc=pd.read_csv(f'{fd_out}/auc_{sample}.csv', index_col=0)
	adata=sc.read(f'{fd_lbl}/{sample}.h5ad')
	df=adata.obs.loc[:, ['anno']].copy().fillna('Unknown')
	l_anno=df['anno']
	
	#2. rss
	dfi_rss=regulon_specificity_scores(dfi_auc, l_anno).T
	dfi_rss.to_csv(f'{fd_out}/rss_{sample}.csv')

















