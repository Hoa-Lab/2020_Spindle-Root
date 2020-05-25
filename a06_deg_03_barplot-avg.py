import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

#------------------------variable------------------------------
fmt='png'
p=0.05
n=30

fd_in='./out/a06_deg_02_clean'
fd_out='./out/a06_deg_03_barplot-avg'

l_name=['spin-other', 'spin-root']
l_sample=['Ctrl', 'MethFix', 'RNAlater']

dic_title={'spin-other': ['Spindle-Root', 'Others'], 'spin-root': ['Spindle-Root-1', 'Spindle-Root-2']}

dic_cmap={'spin-root': ['#c2eb7c', '#ebbfa4'],
          'spin-other': ['#a452e3', '#4287f5']}

#----------------------setup----------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#----------------------function--------------------------------
def load_df(fname, p=p):
	#1. load
	df=pd.read_csv(fname, index_col=0)
	#2. remove high p value
	df=df.loc[df['p']<=p, ['logfc']]
	#3. clean columns
	col=Path(fname).stem.split('_')[0]
	df.columns=[col]
	return df


def barplot_de(df, title, f_out, n=n, c='#a452e3'):
	#1. clean df
	n=min(n, df.shape[0])
	df=df.iloc[0:n, :]
	df.index.name='gene'
	df=df.reset_index()
	low=np.zeros(n)
	#2. plot
	sns.set()
	fig, ax=plt.subplots(figsize=(6, 8))
	ax=sns.barplot(x='Log FC', y='gene', data=df, color=c)
	ax.errorbar(df['Log FC'], df.index, xerr=[low, df['std']], linestyle='none', color='Grey')
	#3. adjust
	plt.title(title, fontsize=24, pad=25, weight='medium')
	plt.ylim([n-0.5, 0.5])
	plt.xlabel('', fontsize=22, labelpad=10)
	plt.ylabel('', fontsize=22, labelpad=10)
	ax.xaxis.tick_top()
	plt.xticks(fontsize=16)
	plt.yticks(fontsize=15, rotation=0, weight='medium')
	plt.ylim([n-0.5, -0.5])
	#4. save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return
	

######################################################################
for name in l_name:   #loop with spin vs other, and spin vs root
	#1. load
	df1=load_df(f'{fd_in}/Ctrl_{name}.csv')
	df2=load_df(f'{fd_in}/MethFix_{name}.csv')
	df3=load_df(f'{fd_in}/RNAlater_{name}.csv')
	
	#2. merge
	df=df1.merge(df2, left_index=True, right_index=True)
	df=df.merge(df3, left_index=True, right_index=True)
	
	#3. sort
	df['Log FC']=df.mean(axis=1)
	df['std']=df.loc[:, l_sample].std(axis=1)
	df=df.sort_values('Log FC', ascending=False)
	df.to_csv(f'{fd_out}/{name}.csv')
	
	#4. split df
	df1=df.loc[df['Log FC']>0, :].copy()
	
	df2=df.loc[df['Log FC']<0, :].copy()
	df2['Log FC']=df2['Log FC']*(-1)
	df2=df2.sort_values('Log FC', ascending=False)
	
	#5. barplot
	title=dic_title[name][0]
	f_out=f'{fd_out}/{name}_1.{fmt}'
	color=dic_cmap[name][0]
	barplot_de(df1, title, f_out, c=color)

	title=dic_title[name][1]
	f_out=f'{fd_out}/{name}_2.{fmt}'
	color=dic_cmap[name][1]
	barplot_de(df2, title, f_out, c=color)
	






