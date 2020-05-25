import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

#----------------------variable------------------------
fmt='tif'
n=20  #rows to plot
o=20  #overlap check

fd_rss='./out/a08_combined-reg_02_rss'
fd_out='./out/a08_combined-reg_07_barplot-rss-avg'

l_sample=['Ctrl', 'MethFix', 'RNAlater']
l_cell=['Root', 'Spindle']
l_color=['#c2eb7c', '#ebbfa4']

#--------------------setup----------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#-------------------function--------------------------
def load_df(fname):
	'''load df, and split to spindle and root
	   df1 is root, df2 is spindle
	'''
	#1. load
	df=pd.read_csv(fname, index_col=0)
	#2. split 
	df1=df.loc[:, ['Spindle-Root-1']].copy()
	df1.columns=[f'{Path(fname).stem.lstrip("rss_")}']
	df2=df.loc[:, ['Spindle-Root-2']].copy()
	df2.columns=[f'{Path(fname).stem.lstrip("rss_")}']
	return df1, df2

def merge_df(l_df):
	#1. merge
	df=l_df[0]
	for dfi in l_df[1:]:
		df=df.merge(dfi, left_index=True, right_index=True)
	#2. calculate avg and std
	df['avg']=df.loc[:, l_sample].mean(axis=1)
	df['std']=df.loc[:, l_sample].std(axis=1)
	#3. sort
	df=df.sort_values('avg', ascending=False)
	return df	
	
def barplot_rss(df, title, f_out, n=n, c='#ebbfa4'):
	#1. clean df
	n=min(n, df.shape[0])
	df=df.iloc[0:n, :]
	df.index.name='Regulon'
	df=df.reset_index()
	low=np.zeros(n)
	#2. plot
	sns.set()
	fig, ax=plt.subplots(figsize=(6, 6))
	ax=sns.barplot(x='avg', y='Regulon', data=df, color=c)
	ax.errorbar(df['avg'], df.index, xerr=[low, df['std']], linestyle='none', color='Grey')
	#3. adjust
	plt.title(title, fontsize=24, pad=25, weight='medium')
	plt.ylim([n-0.5, -0.5])
	plt.xlabel('', fontsize=22, labelpad=10)
	plt.ylabel('', fontsize=18, labelpad=10)
	ax.xaxis.tick_top()
	plt.xticks(fontsize=13)
	plt.yticks(fontsize=14, rotation=0, weight='medium')
	#4. save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return


#######################################################################
#1. load
l_root=[]  
l_spin=[]
for sample in l_sample:
	df1, df2=load_df(f'{fd_rss}/rss_{sample}.csv')
	l_root.append(df1)
	l_spin.append(df2)

#2. merge
df_root=merge_df(l_root)
df_root.to_csv(f'{fd_out}/root.csv')
df_spin=merge_df(l_spin)
df_spin.to_csv(f'{fd_out}/spin.csv')

#3. get uniq regulons
l_overlap=[i for i in df_root.index.tolist()[0:o] if i in df_spin.index.tolist()[0:o]]
df_root=df_root.loc[~df_root.index.isin(l_overlap), :]
df_spin=df_spin.loc[~df_spin.index.isin(l_overlap), :]

#4. save reg for plot
l_reg=df_root.index.tolist()[0:n]+df_spin.index.tolist()[0:n]
Path(f'{fd_out}/reg.txt').write_text('\n'.join(l_reg))


##5. plot
#title=l_cell[0]
#f_out=f'{fd_out}/{title}.{fmt}'
#barplot_rss(df_root, title, f_out, c=l_color[0])

#title=l_cell[1]
#f_out=f'{fd_out}/{title}.{fmt}'
#barplot_rss(df_spin, title, f_out, c=l_color[1])





