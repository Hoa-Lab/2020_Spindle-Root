import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

#----------------------variable------------------------
f_reg='./out/a09_meth-later_09_reg-info/regulon.csv'
fd_reg='./out/a09_meth-later_06_barplot-rss'
fd_out='./out/a12_hoa_00_reg-go'

#l_root=['Sall2', 'Mta3', 'Isl1', 'Zfp850', 'Nr1h4',  'Irf6', 'Zfp712', 'Zfp2', 'Runx1', 'Bcl11a', 'Zkscan14', 'Prox1', 'Hes6']

#l_spin=['Zfp160', 'Ovol2', 'Bach2', 'Meis3', 'Asap3', 'Creb3l4', 'Tgif1', 'Atf6b', 'Barx2', 'Nr1h4', 'Foxg1', 'Irf6', 'Klf6']

#l_shared=['Grhl1', 'Zfp932', 'Zfp595', 'Isl2', 'Sox4', 'Rorb', 'Six1']

#-----------------------setup--------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
df=pd.read_csv(f_reg, index_col=0)

#-----------------------function----------------------
def get_reg(fname):
	l_reg=Path(fname).read_text().split('\n')
	return l_reg


####################################################################
#---------------------------------------------------------
#1. root
name='root'
l_reg=get_reg(f'{fd_reg}/root_reg.txt')

#2. make txt
l_gene=[]

for reg in l_reg:
	dfi=df.loc[reg, :].copy()
	l_gene.extend(dfi['Target'].tolist())
	
l_gene=list(set(l_gene))
Path(f'{fd_out}/{name}.txt').write_text('\n'.join(l_gene))

#-----------------------------------------------------
#1. spindle
name='spin'
l_reg=get_reg(f'{fd_reg}/spin_reg.txt')

#2. make txt
l_gene=[]

for reg in l_reg:
	dfi=df.loc[reg, :].copy()
	l_gene.extend(dfi['Target'].tolist())
	
l_gene=list(set(l_gene))
Path(f'{fd_out}/{name}.txt').write_text('\n'.join(l_gene))


#-------------------------------------------------------
#1. shared
name='shared'
l_reg=get_reg(f'{fd_reg}/shared_reg.txt')

#2. make txt
l_gene=[]

for reg in l_reg:
	dfi=df.loc[reg, :].copy()
	l_gene.extend(dfi['Target'].tolist())
	
l_gene=list(set(l_gene))
Path(f'{fd_out}/{name}.txt').write_text('\n'.join(l_gene))






