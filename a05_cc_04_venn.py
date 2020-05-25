import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import scanpy as sc

#---------------------------variable-------------------------------
fmt='tif'

fd_before='./out/a02_modularity_01_anno'
fd_after='./out/a05_cc_03_anno'
fd_out='./out/a05_cc_04_venn'

l_sample=['Ctrl', 'MethFix', 'RNAlater']
l_cell=['Marginal', 'Intermediate', 'Basal']

#---------------------------setup-------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#------------------function-------------------------
def plot_venn(subset, f_out, title=None):
	#1. setup
	l_lbl=['Original', 'Calibrated']
	#2. plot
	fig, ax=plt.subplots(figsize=(4,4))
	ax=venn2(subsets=subset, set_labels=l_lbl, set_colors=('#c173f5', 'skyblue'), alpha=0.7)
	#3. adjust
	plt.title(title, fontsize=20, pad=6, weight='medium')
	ax.get_patch_by_id('10').set_edgecolor('none')
	ax.get_patch_by_id('01').set_edgecolor('none')
	ax.get_patch_by_id('11').set_edgecolor('none')
	for text in ax.set_labels:     #this is set label
		text.set_fontsize(16)
	for text in ax.subset_labels:  #this is number in circle
		text.set_fontsize(16)
	#4. save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return


###############################################################################
#1. main loop
for sample in l_sample:
	#1. load
	df1=sc.read(f'{fd_before}/{sample}.h5ad').obs.loc[:, ['anno']]
	df2=sc.read(f'{fd_after}/{sample}.h5ad').obs.loc[:, ['anno']].fillna('Unknown')  # why has nan value?
	
	#2. make cell count df
	l_data=[]
	for cell in l_cell:
		#1. cell list
		l_c1=df1.loc[df1['anno']==cell, :].index.tolist()
		l_c2=df2.loc[df2['anno']==cell, :].index.tolist()
		l_overlap=[i for i in l_c1 if i in l_c2]
		n_c1=len(l_c1)
		n_c2=len(l_c2)
		n_overlap=len(l_overlap)
		l_data.append((cell, n_c1, n_c2, n_overlap))
	df=pd.DataFrame(l_data, columns=['cell', 'before', 'after', 'overlap'])
	df=df.set_index('cell')
	df.to_csv(f'{fd_out}/{sample}.csv')
	
	#3. venn
	df['before']=df['before']-df['overlap']
	df['after']=df['after']-df['overlap']
	
	for cell, row in df.iterrows():
		subset=row.values.tolist()
		f_out=f'{fd_out}/{sample}_{cell}.{fmt}'
		title=f'{cell} ({sample})'
		plot_venn(subset, f_out, title=title)
	
	
	
	
	
	
	
	
