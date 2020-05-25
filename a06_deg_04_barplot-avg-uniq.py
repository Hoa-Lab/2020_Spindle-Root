import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

#------------------------variable------------------------------
fmt='tif'
p=0.05
n=20

fd_in='./out/a06_deg_03_barplot-avg'
fd_out='./out/a06_deg_04_barplot-avg-uniq'

#----------------------setup----------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#---------------------function-------------------------------
def barplot_de(df, title, f_out, n=n, c='#a452e3', fs=(5,10)):
	#1. clean df
	n=min(n, df.shape[0])
	df=df.iloc[0:n, :]
	df.index.name='gene'
	df=df.reset_index()
	low=np.zeros(n)
	#2. plot
	sns.set()
	fig, ax=plt.subplots(figsize=fs)
	ax=sns.barplot(x='Log FC', y='gene', data=df, color=c)
	ax.errorbar(df['Log FC'], df.index, xerr=[low, df['std']], linestyle='none', color='Grey')
	#3. adjust
	plt.title(title, fontsize=24, pad=25, weight='medium')
	plt.ylim([n-0.1, -0.5])
	plt.xlabel('', fontsize=22, labelpad=10)
	plt.ylabel('', fontsize=22, labelpad=10)
	ax.xaxis.tick_top()
	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16, rotation=0, weight='medium')
	#4. save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return
	
###########################################################################
#1. load
df_uniq=pd.read_csv(f'{fd_in}/spin-other.csv', index_col=0)
df=pd.read_csv(f'{fd_in}/spin-root.csv', index_col=0)

#2. remove non-specific genes
l_gene=df_uniq.index.tolist()
df=df.loc[df.index.isin(l_gene), :]

#3. split df
df1=df.loc[df['Log FC']>0, :].copy()
df1=df1.sort_values('Log FC', ascending=False)
df1.to_csv(f'{fd_out}/spin1.csv')

df2=df.loc[df['Log FC']<0, :].copy()
df2['Log FC']=df2['Log FC']*(-1)
df2=df2.sort_values('Log FC', ascending=False)
df2.to_csv(f'{fd_out}/spin2.csv')


#4. barplot
title='Spindle-Root-1'
f_out=f'{fd_out}/spin1.{fmt}'
color='#c2eb7c'
barplot_de(df1, title, f_out, c=color)

title='Spindle-Root-2'
f_out=f'{fd_out}/spin2.{fmt}'
color='#ebbfa4'
barplot_de(df2, title, f_out, c=color, fs=(6, 10))







