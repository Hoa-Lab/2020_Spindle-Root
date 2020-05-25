import pandas as pd
import scanpy as sc
from pathlib import Path
import pickle
from harmony import harmonize
import re

#------------------variable----------------------
nn=10
npc=30
n_genes=4000

fd_in='./out/a00_preprocess_00_pp'
fd_anno='./out/a02_modularity_01_anno'
fd_out='./out/a09_meth-later_00_merge'

#-----------------setup-------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)


###############################################################################
#-------------------merge----------------------------------
#1. load
ad1=sc.read(f'{fd_in}/clean_MethFix.h5ad')
ad2=sc.read(f'{fd_in}/clean_RNAlater.h5ad')
ad1_anno=sc.read(f'{fd_anno}/MethFix.h5ad')
ad2_anno=sc.read(f'{fd_anno}/RNAlater.h5ad')

#2. add anno
df1=ad1_anno.obs.loc[:, ['anno']]
df2=ad2_anno.obs.loc[:, ['anno']]

ad1.obs=ad1.obs.merge(df1, left_index=True, right_index=True)
ad2.obs=ad2.obs.merge(df2, left_index=True, right_index=True)

#3. get common genes
l_gene=[i for i in ad1.var.index if i in ad2.var.index]
l_gene=[i for i in l_gene if not re.match('^mt-*', i)]
l_gene=[i for i in l_gene if not re.match('^[A-Z][A-Z]+', i)]
l_gene=[i for i in l_gene if not ('Rik' in i)]
l_gene=[i for i in l_gene if not ('-' in i)]
l_gene=[i for i in l_gene if len(i)>1]
l_gene=[i for i in l_gene if not re.match('^Gm\d+', i)]
l_gene.sort()

#4. merge
ad1=ad1[:, l_gene].copy()
ad2=ad2[:, l_gene].copy()
adata=ad1.concatenate(ad2)
adata.write(f'{fd_out}/clean_merged.h5ad')

#5. Normalization, HVG
sc.pp.normalize_total(adata, exclude_highly_expressed=True)
sc.pp.log1p(adata)
adata.raw=adata

sc.pp.highly_variable_genes(adata, n_top_genes=n_genes)
adata=adata[:, adata.var['highly_variable']].copy()

#6. PCA
sc.pp.regress_out(adata, ['n_counts', 'perc_others'])
sc.tl.pca(adata, svd_solver='arpack')

#7. copy PCA adata for harmony
ad_pca=adata.copy()

#8. calculate neighbor
sc.pp.neighbors(adata, n_neighbors=nn, n_pcs=npc)

#9. embed- umap
sc.tl.umap(adata, n_components=2, random_state=42)
adata.write(f'{fd_out}/concat_merged.h5ad')


#-----------------------------harmony-----------------------------------
#1. rename ad_pca
adata=ad_pca.copy()

#2. harmony
Z=harmonize(adata.obsm['X_pca'], adata.obs, batch_key='sample')
adata.obsm['X_harmony']=Z

#3. calculate neighbor
sc.pp.neighbors(adata, n_neighbors=nn, n_pcs=npc, use_rep='X_harmony')

#4. embed- umap
sc.tl.umap(adata, n_components=2, random_state=42)

#5. save
adata.write(f'{fd_out}/harmony_merged.h5ad')




