import scanpy as sc
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler

#-------------------variable----------------------------
fmt='png'

fd_norm='./out/a00_preprocess_00_pp'
fd_lbl='./out/a02_modularity_01_anno'
fd_gene='./out/a06_deg_04_barplot-avg-uniq'
fd_out='./out/a06_deg_05_hm-uniq'

l_sample=['Ctrl', 'MethFix', 'RNAlater']
cmap=['#c2eb7c', '#ebbfa4']

#---------------------setup----------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
scaler=MinMaxScaler(feature_range=(0, 0.95))

#------------------function----------------------------
def hm_topbar(l_cnt, f_out, cmap=cmap):
	#prepare df
	l=[]
	for i in range(len(l_cnt)):
		l.extend([i]*l_cnt[i])
	df=pd.DataFrame({'n': l}).T
	#plot
	fig, ax=plt.subplots()
	ax=sns.heatmap(df, cmap=cmap, cbar=False)
	#adjust plot
	ax.xaxis.label.set_visible(False)
	ax.yaxis.label.set_visible(False)
	plt.axis('off')
	#save
	plt.tight_layout()
	plt.savefig(f_out)
	plt.close()
	return

##################################################################
#1. get gene list
df1=pd.read_csv(f'{fd_gene}/spin1.csv', index_col=0)
df2=pd.read_csv(f'{fd_gene}/spin2.csv', index_col=0)
df=pd.concat([df1, df2])
l_gene=df.index.tolist()

#2. loop on sample
for sample in l_sample:
	#1. get count data
	ad_cnt=sc.read(f'{fd_norm}/norm_{sample}.h5ad')
	sc.pp.log1p(ad_cnt)
	df=pd.DataFrame(ad_cnt.X.toarray(), index=ad_cnt.obs.index, columns=ad_cnt.var.index)
	df=df.loc[:, l_gene]
	
	#2. get spindle root cells
	ad_lbl=sc.read(f'{fd_lbl}/{sample}.h5ad')
	l1=ad_lbl.obs.loc[ad_lbl.obs['anno']=='Spindle-Root-1', :].index.tolist()
	l2=ad_lbl.obs.loc[ad_lbl.obs['anno']=='Spindle-Root-2', :].index.tolist()
	l_cell=l1+l2
	df=df.reindex(l_cell, axis=0)
	
	#2. scale data
	X=scaler.fit_transform(df.values)
	df=pd.DataFrame(X, index=df.index, columns=df.columns)
	df=df.T
	
	#3. heatmap
	fig, ax=plt.subplots(figsize=(9, 15))
	ax=sns.heatmap(df, cmap='Purples', cbar=False)
	
	#4. adjust plot
	ax.xaxis.label.set_visible(False)
	ax.yaxis.label.set_visible(False)
	plt.xticks([])	
	plt.yticks(fontsize=13, rotation=0, weight='medium')
	
	#5. save
	plt.tight_layout()
	plt.savefig(f'{fd_out}/{sample}.{fmt}', dpi=300)
	plt.close()
	
	#6. plot topbar
	l_cnt=[len(l1), len(l2)]
	f_out=f'{fd_out}/{sample}_topbar.{fmt}'
	hm_topbar(l_cnt, f_out)
	




