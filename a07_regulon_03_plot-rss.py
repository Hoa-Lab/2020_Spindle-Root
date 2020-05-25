import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from adjustText import adjust_text
import numpy as np

#-----------------variable-----------------
fd_rss='./out/a07_regulon_01_rss'
fd_out='./out/a07_regulon_03_plot-rss'

fmt='tif'
n=2
n_bar=20

l_sample=['Ctrl', 'MethFix', 'RNAlater']
l_cell=['Root', 'Spindle']
cmap=['#477ba6', '#fa0000']

#---------------setup----------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)


#----------------function------------------
def plot_rss(df, cell, f_out, cmap=cmap, x='rank', fs=(10,8), n=n, hue='hue', title=None):
	sns.set()
	#1.plot
	fig, ax=plt.subplots(figsize=fs)
	ax=sns.scatterplot(x=x, y=cell, hue=hue, data=df, s=30, linewidth=0, alpha=0.9, palette=cmap)
	#2. adjust
	if not title:
		title=cell
	ax.set_title(title, fontsize=30, pad=25, weight='medium')
	ax.xaxis.set_ticklabels([])
	plt.xlabel('')
	plt.ylabel('Regulon Specificity Score', fontsize=25, labelpad=20, weight='semibold')
	#3.text
	l_x=df[x].tolist()[0:n]
	l_y=df[cell].tolist()[0:n]
	l_txt=df.index.tolist()[0:n]
	texts=[]
	for i in range(n):
		x=l_x[i]+0.1
		y=l_y[i]
		txt=l_txt[i]
		texts.append(plt.text(x, y, txt, fontsize=15, weight='semibold'))
	adjust_text(texts)		
	ax.get_legend().remove()
	#4. save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return


###############################################################
for sample in l_sample:
	#1. load
	df=pd.read_csv(f'{fd_rss}/rss_{sample}.csv', index_col=0)
	
	#2. rename column
	l_col=df.columns
	l_col=[(i, 'Root')[i=='Spindle-Root-1'] for i in l_col]
	l_col=[(i, 'Spindle')[i=='Spindle-Root-2'] for i in l_col]
	df.columns=l_col
	
	#3. plot rank
	for cell in l_cell:
		#1. setup
		dfi=df.loc[:, [cell]].copy()
		dfi=dfi.sort_values(cell, ascending=False)
		dfi['rank']=np.arange(dfi.shape[0])
		dfi['hue']=(dfi['rank']<n)   #only top n reg are True
		
		#2. plot
		f_out=f'{fd_out}/{cell}_{sample}.{fmt}'
		plot_rss(dfi, cell, f_out, title=cell)




