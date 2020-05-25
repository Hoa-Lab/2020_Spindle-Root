#make table of top reg info for each cell types

import pandas as pd
from pathlib import Path
import pickle

#------------------variable-----------------------
n=20

f_rss='./out/a09_meth-later_03_rss/rss_merged.csv'
f_reg='./out/a09_meth-later_02_scenic/regulon.pkl'
f_prune='./out/a09_meth-later_02_scenic/prune.csv'
fd_out='./out/a09_meth-later_09_reg-info'

l_cell=['Spindle', 'Marginal', 'Intermediate', 'Basal', 'Root', 'Reissner']
l_col=['Regulon', 'TF', 'Target', 'Weight', 'nMotifs', 'bestMotif', 'NES']

#---------------setup-------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

df_rss=pd.read_csv(f_rss, index_col=0)

l_col=[('Root', i)[i!='Spindle-Root-1'] for i in df_rss.columns.tolist()]
l_col=[('Spindle', i)[i!='Spindle-Root-2'] for i in l_col]
df_rss.columns=l_col

df_prune=pd.read_csv(f_prune, index_col=0, header=2)

with open(f_reg, "rb") as f:
	l_reg=pickle.load(f)

#---------------function----------------------------
def get_motif(reg, df=df_prune):
	#1. select df
	dfi=df.loc[df.index==reg, :].copy()
	#2. count motif
	n=dfi.shape[0]
	#3. get max value motif
	m=dfi['Unnamed: 6'].max()  #nes score
	best=dfi.loc[dfi['Unnamed: 6']==m,:]['MotifID'].iloc[0]
	#4. get nes
	nes=dfi.loc[dfi['Unnamed: 6']==m,:]['Unnamed: 6'].iloc[0]
	return n, best, nes



###############################################################
##--------------parse regulon------------------------
#l_data=[]
#for reg in l_reg:
#	#1. get some info
#	name=reg.name
#	tf=reg.transcription_factor
#	n_motif, best_motif, nes=get_motif(tf)
#	#2. loop on target
#	dic_weight=reg.gene2weight
#	for target in dic_weight.keys():
#		weight=dic_weight[target]
#		l_data.append((name, tf, target, weight, n_motif, best_motif, nes))

#df=pd.DataFrame(l_data, columns=l_col)
#df.to_csv(f'{fd_out}/regulon.csv', index=False)


####################################################################
df_reg=pd.read_csv(f'{fd_out}/regulon.csv')

#-----------top reg in each cell types--------------------
for cell in l_cell:
	#1. get top reg
	dfi=df_rss.loc[:, [cell]].copy().sort_values(cell, ascending=False)
	l_reg=dfi.index.tolist()[0:n]
	
	#2. get df
	df=df_reg.loc[df_reg['Regulon'].isin(l_reg), :].copy()
	df['Regulon']=pd.Categorical(df['Regulon'], categories=l_reg, ordered=True)
	df=df.sort_values('Regulon')
	df.to_csv(f'{fd_out}/{cell}.csv', index=False)
	





