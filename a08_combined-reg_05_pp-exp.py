#make adata for regulon expressopn plot, need add suffix to each sample first

import pandas as pd
import scanpy as sc
from pathlib import Path

#---------------------variable------------------------------
fd_anno='./out/a02_modularity_01_anno'
f_harm='./out/a02_modularity_02_anno-harmony/harmony_all.h5ad'
fd_auc='./out/a08_combined-reg_02_rss'
fd_out='./out/a08_combined-reg_05_pp-exp'

l_sample=['Ctrl', 'MethFix', 'RNAlater']
dic_suf={'Ctrl':'-1', 'MethFix':'-0-0', 'RNAlater':'-1-0'}

#---------------------setup---------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
ad_harm=sc.read(f_harm)


#########################################################################
#------------------add reg into combined adata------------------------
#1. make auc df
l_df=[]
for sample in l_sample:
	df=pd.read_csv(f'{fd_auc}/auc_{sample}.csv')
	df['Cell']=df['Cell']+dic_suf[sample]
	df=df.set_index('Cell')
	l_df.append(df)
df=pd.concat(l_df)

#2. merge to ad_harm
ad_harm.obs=ad_harm.obs.merge(df, left_index=True, right_index=True)
ad_harm.write(f'{fd_out}/harmony_all.h5ad')


#------------------add reg into each adata------------------------
for sample in l_sample:
	#1. load
	adata=sc.read(f'{fd_anno}/{sample}.h5ad')
	df=pd.read_csv(f'{fd_auc}/auc_{sample}.csv', index_col=0)
	
	#2. merge
	adata.obs=adata.obs.merge(df, left_index=True, right_index=True)
	adata.write(f'{fd_out}/{sample}.h5ad')
	
	print(adata.obs.head())
