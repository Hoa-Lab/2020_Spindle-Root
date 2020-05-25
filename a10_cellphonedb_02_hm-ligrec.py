import pandas as pd
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

#----------------------variable-----------------------------
fmt='png'
s=0.3     #min score to keep
m=5     #max cbar

f_in='./out/a10_cellphonedb_00_clean/out/significant_means.txt'
fd_out='./out/a10_cellphonedb_02_hm-ligrec'

#----------------------------setup------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
df_pair=pd.read_csv(f_in, sep='\t', index_col=0).fillna(0)

#---------------------------function------------------------------
def extract_df(lig, rec, dfi=df_pair):
	'''
	extract non-zero pairs
	'''
	#1. setup
	df=dfi.copy()
	l_col=['gene_a', 'gene_b',  'receptor_a', 'receptor_b', f'{lig}|{rec}'] #, f'{rec}|{lig}']
	#2. clean df
	df=df.loc[:, l_col]
	df=df.loc[(df['gene_a']!=0) & (df['gene_b']!=0), :]
	df['gene_a']=df['gene_a'].str.capitalize()
	df['gene_b']=df['gene_b'].str.capitalize()
	#3. remove all 0
	df=df.loc[df[f'{lig}|{rec}']!=0, :]
	return df


###########################################################################
#-----------------------margina vs root--------------------------
c1='Marginal'
c2='Root'

#1. extract df
df1=extract_df(c1, c2)
df2=extract_df(c2, c1)

#2. concat
df1=df1.loc[df1[f'{c1}|{c2}']>s, ['gene_a', 'gene_b', f'{c1}|{c2}']]
df1.columns=['g1', 'g2', 'score']

df2=df2.loc[df2[f'{c2}|{c1}']>s, ['gene_b', 'gene_a', f'{c2}|{c1}']]
df2.columns=['g1', 'g2', 'score']

df=pd.concat([df1, df2])

#3. pivot table
df=pd.pivot_table(df, index=['g1'], columns=['g2'], values='score').fillna(0)

#4. heatmap
fig, ax=plt.subplots(figsize=(15,14))
ax=sns.heatmap(df, cmap='Purples', vmin=s, vmax=m)

plt.xlabel(c2, fontsize=34, weight='semibold', labelpad=10)
plt.ylabel(c1, fontsize=34, weight='semibold', labelpad=10)
plt.xticks(fontsize=16, rotation=45)#, weight='semibold')
plt.yticks(fontsize=22, rotation=0, weight='medium')

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=18)

plt.tight_layout()
plt.savefig(f'{fd_out}/{c1}_{c2}.{fmt}', dpi=300)
plt.close()



#-----------------------margina vs spindle--------------------------
c1='Marginal'
c2='Spindle'

#1. extract df
df1=extract_df(c1, c2)
df2=extract_df(c2, c1)

#2. concat
df1=df1.loc[df1[f'{c1}|{c2}']>s, ['gene_a', 'gene_b', f'{c1}|{c2}']]
df1.columns=['g1', 'g2', 'score']

df2=df2.loc[df2[f'{c2}|{c1}']>s, ['gene_b', 'gene_a', f'{c2}|{c1}']]
df2.columns=['g1', 'g2', 'score']

df=pd.concat([df1, df2])

#3. pivot table
df=pd.pivot_table(df, index=['g1'], columns=['g2'], values='score').fillna(0)

#4. heatmap
fig, ax=plt.subplots(figsize=(15,14))
ax=sns.heatmap(df, cmap='Purples', vmin=s, vmax=m)

plt.xlabel(c2, fontsize=34, weight='semibold', labelpad=10)
plt.ylabel(c1, fontsize=34, weight='semibold', labelpad=10)
plt.xticks(fontsize=16, rotation=45)#, weight='semibold')
plt.yticks(fontsize=22, rotation=0, weight='medium')

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=18)

plt.tight_layout()
plt.savefig(f'{fd_out}/{c1}_{c2}.{fmt}', dpi=300)
plt.close()



#-----------------------root vs spindle--------------------------
c1='Root'
c2='Spindle'

#1. extract df
df1=extract_df(c1, c2)
df2=extract_df(c2, c1)

#2. concat
df1=df1.loc[df1[f'{c1}|{c2}']>s, ['gene_a', 'gene_b', f'{c1}|{c2}']]
df1.columns=['g1', 'g2', 'score']

df2=df2.loc[df2[f'{c2}|{c1}']>s, ['gene_b', 'gene_a', f'{c2}|{c1}']]
df2.columns=['g1', 'g2', 'score']

df=pd.concat([df1, df2])

#3. pivot table
df=pd.pivot_table(df, index=['g1'], columns=['g2'], values='score').fillna(0)

#4. heatmap
fig, ax=plt.subplots(figsize=(15,14))
ax=sns.heatmap(df, cmap='Purples', vmin=s, vmax=m)

plt.xlabel(c2, fontsize=34, weight='semibold', labelpad=10)
plt.ylabel(c1, fontsize=34, weight='semibold', labelpad=10)
plt.xticks(fontsize=16, rotation=45)#, weight='semibold')
plt.yticks(fontsize=22, rotation=0, weight='medium')

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=18)

plt.tight_layout()
plt.savefig(f'{fd_out}/{c1}_{c2}.{fmt}', dpi=300)
plt.close()

