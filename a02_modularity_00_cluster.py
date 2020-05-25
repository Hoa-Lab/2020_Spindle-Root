import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path

#--------------------variable--------------------------
fmt='tif'

fd_in='./out/a00_preprocess_00_pp'
fd_out='./out/a02_modularity_00_cluster'

l_sample=['Ctrl', 'MethFix', 'RNAlater']
dic_res={'Ctrl': 1, 'MethFix': 1, 'RNAlater': 1.5}

#----------------------setup--------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)



########################################################################
#1. main loop
for sample in l_sample:
	#1. load
	adata=sc.read(f'{fd_in}/embed_{sample}.h5ad')
	
	#2. cluster
	sc.tl.leiden(adata, resolution=dic_res[sample])
	adata.write(f'{fd_out}/{sample}.h5ad')
	
	#3. plot
	ax=sc.pl.umap(adata, color=['leiden'], legend_loc='on data', s=8, alpha=0.8, palette='tab20', show=False, frameon=False)

	ax.set_title(sample, fontsize=24, pad=20, weight='medium')
	plt.tight_layout()
	plt.savefig(f'{fd_out}/{sample}.{fmt}', dpi=300)
	plt.close()
	

