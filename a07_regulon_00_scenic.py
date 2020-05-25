import scanpy as sc
import pandas as pd
from pathlib import Path
import os
import time
import pickle
import numpy as np
from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

#---------------------variable--------------------------
sample='RNAlater'

fd_in='./out/a00_preprocess_00_pp'
fd_db='./raw/ref/scenic/db'
f_tf='./raw/ref/scenic/mm_mgi_tfs.txt'
f_motif='./raw/ref/scenic/motifs-v9-nr.mgi-m0.001-o0.0.tbl'
fd_out='./out/a07_regulon_00_scenic'

#--------------------setup-----------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#---------------functions-----------------
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
    
#################################################################
#1. make df
adata=sc.read(f'{fd_in}/clean_{sample}.h5ad')
df_cnt=pd.DataFrame(adata.X.toarray(), index=adata.obs.index, columns=adata.var.index)

#2. tf genes
tf_name=load_tf_names(f_tf)

#3. ranking databases (only 2 mm10 dbs)
l_fname=list(Path(fd_db).glob('*.feather'))
l_db=[RankingDatabase(fname=fname, name=name(fname)) for fname in l_fname]

#3. run
if __name__ =='__main__':
	#1. Inference of co-expression modules
	print('Inference...')
	df_adj=grnboost2(df_cnt, tf_names=tf_name, verbose=True)
	df_adj.to_csv(f'{fd_out}/adj_{sample}.csv', index=False)
	
	#2. prune
	df_adj=pd.read_csv(f'{fd_out}/adj_{sample}.csv')  #if missing, always stuck at 98%
	print('Prune...')
	l_mod=list(modules_from_adjacencies(df_adj, df_cnt))

	with ProgressBar():
		df_prune = prune2df(l_db, l_mod, f_motif)
	df_prune.to_csv(f'{fd_out}/prune_{sample}.csv')
	
	#3. create regulon
	print('Regulon...')
	regulon=df2regulons(df_prune)

	#4. Save the enriched motifs and the discovered regulons
	with open(f'{fd_out}/regulon_{sample}.pkl', "wb") as f:
		pickle.dump(regulon, f)
	
	#5. auc
	print('AUC...')
	with open(f'{fd_out}/regulon_{sample}.pkl', "rb") as f:   #if missing, always stuck
		regulon=pickle.load(f)
		
	df_auc=aucell(df_cnt, regulon, num_workers=10)
	df_auc.to_csv(f'{fd_out}/auc_{sample}.csv')















