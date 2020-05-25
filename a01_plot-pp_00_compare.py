# compare gene numbers in different samples

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import ttest_ind
import scipy.stats as stats
import scikit_posthocs as sp

#-------------------variable--------------------------------
fmt='png'

fd_in='./out/a00_preprocess_00_pp'
fd_out='./out/a01_plot-pp_00_compare'

l_sample=['Ctrl', 'MethFix', 'RNAlater']
dic_cmap={'Ctrl': '#4287f5', 'MethFix': '#f5a142', 'RNAlater': '#4bf542'}

#--------------------setup---------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)


#--------------------function----------------------------
def get_gene_cnt(prefix, l_sample=l_sample, fd_in=fd_in):
	'''count genes in each adata, and concat count dfs
	'''
	#1. add gene count df to list
	l_df=[]
	for sample in l_sample:
		adata=sc.read(f'{fd_in}/{prefix}_{sample}.h5ad')
		sc.pp.filter_cells(adata, min_genes=0)  #this will count genes in each cell
		l_df.append(adata.obs)	
	#2. concat df
	df=pd.concat(l_df)
	return df


def plot_gene(df, f_out, title, dic_cmap=dic_cmap, ylim=None):
	#1. plot
	sns.set()
	fig, ax=plt.subplots(figsize=(8, 5))
	ax=sns.violinplot(x='sample', y='n_genes', data=df, hue='sample', linewidth=0.5, width=1.5, palette=dic_cmap)
	#2. adjust
	ax.set_title(title, fontsize=20, pad=15, weight='medium')
	plt.xlabel('')
	plt.ylabel('Gene Number', fontsize=22, labelpad=15, weight='medium')
	plt.xticks([-0.5, 1, 2.5], fontsize=22, rotation=0, va='center')
	ax.tick_params(axis='x', which='major', pad=15)
	plt.xlim([-1, 3])
	plt.ylim(ylim)
	ax.get_legend().remove()
	#3. save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return


############################################################################
#----------------------raw data------------------------------
prefix='raw'

#1. count df
df=get_gene_cnt(prefix)

##2. plot
#f_out=f'{fd_out}/{prefix}_gene.{fmt}'
#title=f'Gene Numbers ({prefix.capitalize()})'
#plot_gene(df, f_out, title, ylim=[-1000, 16000])

##3. calculate p value
#ctrl=df.loc[df['sample']=='Ctrl']['n_genes']
#meth=df.loc[df['sample']=='MethFix']['n_genes']
#later=df.loc[df['sample']=='RNAlater']['n_genes']
#t1, p1=ttest_ind(ctrl, meth)
#print(p1)   #0
#t2, p2=ttest_ind(ctrl, later)
#print(p2)   #0
#t3, p3=ttest_ind(meth, later)
#print(p3)   #0.03836


#---------------------anova------------------------------------
##1. get data
#l_ctrl=df.loc[df['sample']=='Ctrl', ['n_genes']]['n_genes'].tolist()
#l_meth=df.loc[df['sample']=='MethFix', ['n_genes']]['n_genes'].tolist()
#l_later=df.loc[df['sample']=='RNAlater', ['n_genes']]['n_genes'].tolist()

#l_all=[l_ctrl, l_meth, l_later]

##2. avova
#fvalue, pvalue=stats.f_oneway(l_ctrl, l_meth, l_later)
#print(fvalue, pvalue)  #4148.3173795985 0.0

##3. post hoc ttest
#p=sp.posthoc_conover(l_all, p_adjust='holm')
#print(p)

'''
     1         2         3
1 -1.0  0.000000  0.000000
2  0.0 -1.000000  0.880754
3  0.0  0.880754 -1.000000
'''

#########################################################################
##----------------------cleaned data------------------------------
#prefix='clean'

##1. count df
#df=get_gene_cnt(prefix)

##2. plot
#f_out=f'{fd_out}/{prefix}_gene.{fmt}'
#title=f'Gene Numbers ({prefix.capitalize()})'
#plot_gene(df, f_out, title, ylim=[0, 4500])

##3. calculate p value
#ctrl=df.loc[df['sample']=='Ctrl']['n_genes']
#meth=df.loc[df['sample']=='MethFix']['n_genes']
#later=df.loc[df['sample']=='RNAlater']['n_genes']
#t1, p1=ttest_ind(ctrl, meth)
#print(p1)   #0
#t2, p2=ttest_ind(ctrl, later)
#print(p2)   #0
#t3, p3=ttest_ind(meth, later)
#print(p3)   #6.629720921642305e-78












