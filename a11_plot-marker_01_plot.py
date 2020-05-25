import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import ttest_ind

#--------------------------variable------------------------
fmt='png'

f_csv='./out/a11_plot-marker_00_pp/embed-exp_all.csv'
fd_out='./out/a11_plot-marker_01_plot'

cmap=['#c2eb7c', '#ebbfa4']   #root, spindle

#--------------------setup---------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
df_all=pd.read_csv(f_csv, index_col=0)
df_all['anno']=pd.Categorical(df_all['anno'], categories=['Spindle-Root-1', 'Spindle-Root-2'], ordered=True)


####################################################################
#-------------------root------------------------------
#1. prepare
f_out=f'{fd_out}/violin_root.{fmt}'
l_gene=['Lgr5', 'Epyc']

df=df_all.loc[df_all['gene'].isin(l_gene), :].copy()

#2. plot
sns.set()
fig, ax=plt.subplots(figsize=(10,6))
ax=sns.violinplot(x='gene', y='exp', hue='anno', data=df, palette=cmap, scale='count', dodge=True)

#3. adjust
ax.set_title('Spindle-Root-1 Markers', fontsize=27, pad=20, weight='semibold')
plt.xlabel('', fontsize=22, labelpad=10)
plt.xticks(fontsize=28, rotation=0, weight='semibold')
plt.ylabel('Normalized Counts', fontsize=22, labelpad=15, weight='semibold')
plt.yticks(fontsize=18, rotation=0)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, prop={'size': 15, 'weight':'semibold'})

#4. save
plt.tight_layout()
plt.savefig(f_out, dpi=300)
plt.close()


#-------------------spindle------------------------------
#1. prepare
f_out=f'{fd_out}/violin_spin.{fmt}'
l_gene=['Anxa1', 'Dpp10']

df=df_all.loc[df_all['gene'].isin(l_gene), :].copy()

#2. plot
sns.set()
fig, ax=plt.subplots(figsize=(10,6))
ax=sns.violinplot(x='gene', y='exp', hue='anno', data=df, palette=cmap, scale='count', dodge=True)

#3. adjust
ax.set_title('Spindle-Root-2 Markers', fontsize=27, pad=20, weight='semibold')
plt.xlabel('', fontsize=22, labelpad=10)
plt.xticks(fontsize=28, rotation=0, weight='semibold')
plt.ylabel('Normalized Counts', fontsize=22, labelpad=15, weight='semibold')
plt.yticks(fontsize=18, rotation=0)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, prop={'size': 15,  'weight':'semibold'})

#4. save
plt.tight_layout()
plt.savefig(f_out, dpi=300)
plt.close()
