import scanpy as sc
import pandas as pd
from pathlib import Path
import re

#------------------------variable------------------------------
fd_in='./out/a00_preprocess_00_pp'
fd_out='./out/a08_combined-reg_00_pp'

l_sample=['Ctrl', 'MethFix', 'RNAlater']

#--------------------setup-----------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)


#####################################################################
#1. load
ad1=sc.read(f'{fd_in}/clean_Ctrl.h5ad')
ad2=sc.read(f'{fd_in}/clean_MethFix.h5ad')
ad3=sc.read(f'{fd_in}/clean_RNAlater.h5ad')

#2. get common genes
l_gene=ad1.var.index.tolist()
l_gene=[i for i in l_gene if i in ad2.var.index.tolist()]
l_gene=[i for i in l_gene if i in ad3.var.index.tolist()]

#3. clean gene
l_gene=[i for i in l_gene if not re.match('^mt-*', i)]
l_gene=[i for i in l_gene if not re.match('^[A-Z][A-Z]+', i)]
l_gene=[i for i in l_gene if not ('Rik' in i)]
l_gene=[i for i in l_gene if not ('-' in i)]
l_gene=[i for i in l_gene if len(i)>1]
l_gene=[i for i in l_gene if not re.match('^Gm\d+', i)]
l_gene.sort()

#4. make df
df1=pd.DataFrame(ad1.X.toarray(), index=ad1.obs.index, columns=ad1.var.index)
df1=df1.loc[:, l_gene]

df2=pd.DataFrame(ad2.X.toarray(), index=ad2.obs.index, columns=ad2.var.index)
df2=df2.loc[:, l_gene]

df3=pd.DataFrame(ad3.X.toarray(), index=ad3.obs.index, columns=ad3.var.index)
df3=df3.loc[:, l_gene]

#5. concat df
df=pd.concat([df1, df2, df3])

df.to_csv(f'{fd_out}/all_sample.csv')
print(df.shape)







