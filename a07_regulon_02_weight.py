import pandas as pd
import numpy as np
from pathlib import Path
import pickle

#---------------------variable------------------------------
fd_in='./out/a07_regulon_00_scenic'
fd_out='./out/a07_regulon_02_weight'

l_sample=['Ctrl', 'MethFix', 'RNAlater']

#---------------------setup----------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)


####################################################################
for sample in l_sample:
	#1. load
	with open(f'{fd_in}/regulon_{sample}.pkl', "rb") as f:
		l_reg=pickle.load(f)
	
	#2. make df
	l_data=[]
	for reg in l_reg:
		name=reg.name
		t_genes=reg.genes
		d_weight=reg.gene2weight
		for gene in t_genes:
			l_data.append((name, gene, d_weight[gene]))
	
	df=pd.DataFrame(l_data, columns=['regulon', 'gene', 'weight'])
	
	#3. save
	df.to_csv(f'{fd_out}/weight_{sample}.csv', index=False)
	
