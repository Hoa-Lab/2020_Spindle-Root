import pandas as pd
import scanpy as sc
from pathlib import Path

#-------------------------variable---------------------------
fd_cnt='./out/a00_preprocess_00_pp'
fd_lbl='./out/a02_modularity_01_anno'
fd_out='./out/a06_deg_00_pp'

l_sample=['Ctrl', 'MethFix', 'RNAlater']

#------------------------setup-------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)


################################################################################
#------------------spindle-root vs others--------------------
for sample in l_sample:
	#1. setup
	ad_cnt=sc.read(f'{fd_cnt}/clean_{sample}.h5ad')
	ad_lbl=sc.read(f'{fd_lbl}/{sample}.h5ad')
	ad_lbl.obs['anno']=ad_lbl.obs['anno'].astype('str').fillna('Unknown')
	
	#2. split cell tags
	df_spin=ad_lbl.obs.loc[ad_lbl.obs['anno'].str.contains('Spindle-Root'), :]
	l_spin=df_spin.index.tolist()
	
	df_other=ad_lbl.obs.loc[~ad_lbl.obs['anno'].str.contains('Spindle-Root'), :]
	l_other=df_other.index.tolist()
	
	#3. extract counts
	df=pd.DataFrame(ad_cnt.X.toarray(), index=ad_cnt.obs.index, columns=ad_cnt.var.index)
	df1=df.reindex(l_spin, axis=0)
	df2=df.reindex(l_other, axis=0)
	
	#4. concat
	df=pd.concat([df1, df2]).T
	df.to_csv(f'{fd_out}/{sample}_spin-other_{df1.shape[0]}_{df2.shape[0]}.csv')

	
#------------------spindle vs root----------------------------
for sample in l_sample:
	#1. setup
	ad_cnt=sc.read(f'{fd_cnt}/clean_{sample}.h5ad')
	ad_lbl=sc.read(f'{fd_lbl}/{sample}.h5ad')
	ad_lbl.obs['anno']=ad_lbl.obs['anno'].astype('str').fillna('Unknown')
	
	#2. split cell tags
	df_sr1=ad_lbl.obs.loc[ad_lbl.obs['anno']=='Spindle-Root-1', :]
	l_sr1=df_sr1.index.tolist()
	
	df_sr2=ad_lbl.obs.loc[ad_lbl.obs['anno']=='Spindle-Root-2', :]
	l_sr2=df_sr2.index.tolist()
	
	#3. extract counts
	df=pd.DataFrame(ad_cnt.X.toarray(), index=ad_cnt.obs.index, columns=ad_cnt.var.index)
	df1=df.reindex(l_sr1, axis=0)
	df2=df.reindex(l_sr2, axis=0)
	
	#4. concat
	df=pd.concat([df1, df2]).T
	df.to_csv(f'{fd_out}/{sample}_spin-root_{df1.shape[0]}_{df2.shape[0]}.csv')



