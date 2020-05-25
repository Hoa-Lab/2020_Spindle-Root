import scanpy as sc
from pathlib import Path

#------------------------variable------------------------------
nn=10
npc=30
n_genes=4000

fd_in='./out/a00_preprocess_00_pp'
f_gl='./raw/ref/gl_cell_cycle.txt'
fd_out='./out/a05_cc_00_pp'

l_sample=['Ctrl', 'MethFix', 'RNAlater']

#---------------------setup-----------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
l_gene=Path(f_gl).read_text().split('\n')


########################################################################
#1. main loop
for sample in l_sample:
	#1. setup
	adata=sc.read(f'{fd_in}/clean_{sample}.h5ad')

	#2. get valid gene names
	l_diss=[i for i in l_gene if i in adata.var_names]
	
	#3. normalization
	sc.pp.normalize_total(adata, exclude_highly_expressed=True)
	sc.pp.log1p(adata)
	adata.raw=adata
	sc.pp.scale(adata)

	#4. score
	sc.tl.score_genes(adata, gene_list=l_diss, score_name='cc')

	#5. HVG
	adata=adata[:,adata.X.sum(axis=0) > 0].copy()
	sc.pp.highly_variable_genes(adata, n_top_genes=n_genes)
	adata=adata[:, adata.var['highly_variable']]

	#6. PCA
	sc.pp.regress_out(adata, ['n_counts', 'perc_others', 'cc'])
	sc.tl.pca(adata, svd_solver='arpack')

	#7. calculate neighbor
	sc.pp.neighbors(adata, n_neighbors=nn, n_pcs=npc)

	#8. embed- umap
	sc.tl.umap(adata, n_components=2)
	
	adata.write(f'{fd_out}/{sample}.h5ad')










