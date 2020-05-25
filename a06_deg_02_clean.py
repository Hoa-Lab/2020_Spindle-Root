import scanpy as sc
from pathlib import Path
import pandas as pd
import numpy as np

#--------------------variable--------------------------
fd_in='./out/a06_deg_01_de'
fd_out='./out/a06_deg_02_clean'

#-------------------setup-------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
l_fname=Path(fd_in).glob('*.csv')

#---------------------function------------------------
def clean_deg(df):
	#1. load	
	df=df.loc[:, ['norm_foldChange','pvalue.adj.FDR']]
	df.columns=['fc', 'p']
	#2. conver inf fc
	x=df.loc[df['fc']!=np.inf]['fc'].max()
	df.loc[df['fc']==np.inf, ['fc']]=x
	#3. convert 0 to min
	y=df.loc[df['fc']!=0]['fc'].min()
	df.loc[df['fc']==0, ['fc']]=y
	#4. add log fc
	df['logfc']=df['fc'].apply(np.log2)
	return df


################################################################
for fname in l_fname:
	#1. setup
	name=Path(fname).stem
	df=pd.read_csv(fname, index_col=0)
	
	#2. clean
	df=clean_deg(df)
	df.to_csv(f'{fd_out}/{name}.csv')
	


