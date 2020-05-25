import pandas as pd
import numpy as np
from pathlib import Path
import pickle
import pprint

#---------------------variable------------------------------
fd_in='./out/a09_meth-later_02_scenic'
fd_out='./out/a09_meth-later_04_weight'

#---------------------setup----------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
pp = pprint.PrettyPrinter(indent=4)

################################################################
#1. load
with open(f'{fd_in}/regulon.pkl', "rb") as f:
	l_reg=pickle.load(f)


##2. make df
#l_data=[]
#for reg in l_reg:
#	name=reg.name
#	t_genes=reg.genes
#	d_weight=reg.gene2weight
#	for gene in t_genes:
#		l_data.append((name, gene, d_weight[gene]))

#df=pd.DataFrame(l_data, columns=['regulon', 'gene', 'weight'])

##3. save
#df.to_csv(f'{fd_out}/weight.csv', index=False)
