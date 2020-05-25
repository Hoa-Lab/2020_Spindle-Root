import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

#----------------------variable------------------------
fmt='png'
n=10  #rows to plot
o=20  #overlap check

f_rss='./out/a09_meth-later_03_rss/rss_merged.csv'
fd_out='./out/a09_meth-later_06_barplot-rss'

l_cell=['Root', 'Spindle']
l_color=['#c2eb7c', '#ebbfa4', '#a452e3']

#--------------------setup----------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#-------------------function--------------------------
def barplot_rss(df, title, f_out, c='#ebbfa4', n=n):
	#1. clean df
	df.columns=['rss']
	df.index.name='Regulon'
	df=df.reset_index()
	#2. plot
	sns.set()
	fig, ax=plt.subplots(figsize=(5, 6))
	ax=sns.barplot(x='rss', y='Regulon', data=df, color=c)
	#ax.errorbar(df['avg'], df.index, xerr=[low, df['std']], linestyle='none', color='Grey')
	#3. adjust
	plt.title('', fontsize=24, pad=25, weight='medium')
	plt.ylim([n-0.5, -0.5])
	plt.xlabel('', fontsize=22, labelpad=10)
	plt.ylabel('', fontsize=18, labelpad=10)
	ax.xaxis.tick_top()
	plt.xticks(fontsize=13)
	plt.yticks(fontsize=16, rotation=0, weight='semibold')
	#4. save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return


########################################################################
#1. load
df_rss=pd.read_csv(f_rss, index_col=0)
l_col=df_rss.columns.tolist()
l_col=[('Root', i)[i!='Spindle-Root-1'] for i in l_col]
l_col=[('Spindle', i)[i!='Spindle-Root-2'] for i in l_col]
df_rss.columns=l_col

#2. find overlap reg in spin and root
df_root=df_rss.sort_values('Root', ascending=False).loc[:, ['Root']].copy()
df_spin=df_rss.sort_values('Spindle', ascending=False).loc[:, ['Spindle']].copy()

df_root.to_csv(f'{fd_out}/root.csv')
df_spin.to_csv(f'{fd_out}/spin.csv')

l_overlap=[i for i in df_root.index.tolist()[0:o] if i in df_spin.index.tolist()[0:o]]
Path(f'{fd_out}/shared_reg.txt').write_text('\n'.join(l_overlap))

#3. save reg for exp plot
l_reg=df_root.index.tolist()[0:o]+df_spin.index.tolist()[0:o]
l_reg=list(set(l_reg))
Path(f'{fd_out}/all_reg.txt').write_text('\n'.join(l_reg))

#4. filter overlap regulons (optional)
#l_overlap=[]  
df_overlap=df_root.loc[df_root.index.isin(l_overlap), :].copy()
df_root=df_root.loc[~df_root.index.isin(l_overlap), :].iloc[0:n]
df_spin=df_spin.loc[~df_spin.index.isin(l_overlap), :].iloc[0:n]

Path(f'{fd_out}/root_reg.txt').write_text('\n'.join(df_root.index.tolist()))
Path(f'{fd_out}/spin_reg.txt').write_text('\n'.join(df_spin.index.tolist()))

#4. plot
f_out=f'{fd_out}/root.{fmt}'
c=l_color[0]
title='Root'
barplot_rss(df_root, title, f_out, c=c)

f_out=f'{fd_out}/spin.{fmt}'
c=l_color[1]
title='Spindle'
barplot_rss(df_spin, title, f_out, c=c)

f_out=f'{fd_out}/shared.{fmt}'
c=l_color[2]
title='Shared'
barplot_rss(df_overlap, title, f_out, c=c, n=9)














