import pandas as pd
import scanpy as sc
from pathlib import Path
from pyscenic.rss import regulon_specificity_scores

#-----------------------variable-----------------------------
f_lbl='./out/a09_meth-later_00_merge/clean_merged.h5ad'
f_auc='./out/a09_meth-later_02_scenic/auc.csv'
fd_out='./out/a09_meth-later_03_rss'

#---------------------setup----------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
df_auc=pd.read_csv(f_auc, index_col=0)
adata=sc.read(f_lbl)


######################################################################
#1. get l_anno
l_anno=adata.obs['anno'].fillna('Unknown')

#2. RSS
df_rss=regulon_specificity_scores(df_auc, l_anno).T
df_rss.to_csv(f'{fd_out}/rss_merged.csv')
