# use normalized and log1p data
# concat adata
# MNN adata


import pandas as pd
import scanpy as sc
from pathlib import Path
import pickle
from harmony import harmonize

#---------------------variable----------------------------
fd_in='./out/a00_preprocess_00_pp'
fd_out='./out/a00_preprocess_01_harmony'

l_sample=['MethFix', 'RNAlater', 'Ctrl']  #merge meth and later first
nn=10
npc=30
n_genes=4000

#------------------setup-------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#-----------------function-----------------------------
def load_pkl(fname):
	with open(fname, 'rb') as f:
		pkl=pickle.load(f)
	return pkl

def dump_pkl(obj, fname):
	with open(fname, 'wb') as f:
		pickle.dump(obj, f)

####################################################################
#---------------------concat data------------------------------
#1. load normalized adata to dic_ad, and save all genes
dic_ad={}
l_gene=[]
for sample in l_sample:
	#1. load
	adata=sc.read(f'{fd_in}/norm_{sample}.h5ad')
	dic_ad[sample]=adata
	l_gene.extend(adata.var.index.tolist())

#2. get common genes
l_gene=list(set(l_gene))
for sample in l_sample:
	l_gene=[i for i in l_gene if i in dic_ad[sample].var.index.tolist()]
l_gene.sort()

#3. filter genes in adata, log1p and add to l_ad
l_ad=[]
for sample in l_sample:
	adata=dic_ad[sample][:, l_gene].copy()
	sc.pp.log1p(adata)
	l_ad.append(adata)

#4. concat adata without integration
adata=l_ad[0]
for adi in l_ad[1:]:
	adata=adata.concatenate(adi)

#5. HVG
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

#10. save
adata.write(f'{fd_out}/concat_all.h5ad')


##-----------------------------harmony-----------------------------------
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
adata.write(f'{fd_out}/harmony_all.h5ad')











