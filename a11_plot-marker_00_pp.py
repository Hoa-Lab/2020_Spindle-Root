# violin plot on spindle-root marker genes

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import ttest_ind

#--------------------------variable------------------------
fmt='png'

fd_cnt='./out/a00_preprocess_00_pp'
fd_anno='./out/a02_modularity_01_anno'
fd_out='./out/a11_plot-marker_00_pp'

l_gene=['Lgr5', 'Epyc', 'Anxa1', 'Dpp10']
l_sample=['Ctrl', 'MethFix', 'RNAlater']
cmap=['#c2eb7c', '#ebbfa4']   #root, spindle

#--------------------setup---------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)


#############################################################################
#-------------------prepare df----------------------------
#1. df list
l_df=[]
for sample in l_sample:
	#1. get exp df
	ad_cnt=sc.read(f'{fd_cnt}/embed_{sample}.h5ad')
	ad_anno=sc.read(f'{fd_anno}/{sample}.h5ad')
	
	df=pd.DataFrame(ad_cnt.X, index=ad_cnt.obs.index, columns=ad_cnt.var.index)
	df_anno=ad_anno.obs.loc[:, ['anno']]
	
	#2. clean df
	df=df.loc[:, l_gene]
	df=df.merge(df_anno, left_index=True, right_index=True)
	df=df.loc[(df['anno']=='Spindle-Root-1') | (df['anno']=='Spindle-Root-2'), :].copy()
	df['sample']=sample
	l_df.append(df)
	
#2. concat
df=pd.concat(l_df)

#3. stack genes
l_df=[]
for gene in l_gene:
	dfi=df.loc[:, [gene, 'anno', 'sample']].copy()
	dfi.columns=['exp', 'anno', 'sample']
	dfi['gene']=gene
	l_df.append(dfi)
	
df=pd.concat(l_df)
df['anno']=pd.Categorical(df['anno'], categories=['Spindle-Root-1', 'Spindle-Root-2'], ordered=True)
df.to_csv(f'{fd_out}/embed-exp_all.csv')






