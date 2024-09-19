# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 23:21:15 2021

@author: congt
"""


import pandas as pd
path = '/path/to/file/'
df = pd.read_csv(path+'mRNA_seq_GAPDH_normalized_WT_1A5_21A4_29A3.bed', sep='\t')
df = df[['Gene', 'WT', '1A5', '21A4']]
df['KI/WT'] = df['21A4'] / df['WT']
df['KO/WT'] = df['1A5'] / df['WT']

#%%
'''
up- down-regluated in genes in KI cells
'''
df_select = df[df['KI/WT'] > 1.4]
df_select.to_csv(path+'upregulated_gene_in_KI_21A4.txt', sep='\t', index=False)

df_select = df[df['KI/WT'] < 1/1.4]
df_select.to_csv(path+'downregulated_gene_in_KI_21A4.txt', sep='\t', index=False)

#%%
'''
plot Gene set enrichment analysis
1. Submit gene list in upregulated/downregulated groups in the previous commands to Enrichr: https://maayanlab.cloud/Enrichr/
2. Download data in GO Biological Process tab
3. Plot top 10 processes (lowest p-val) value using barplot
'''
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

ax = plt.figure(figsize=(5,4))
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False

data = pd.read_csv(path+'GO_Biological_Process_2023_table_upregulated_genes_in_KI_(21A4).txt', sep='\t')
group = 'upregulation' #upregulation or downregulation
data = data[['Term', 'P-value']]
data.sort_values(['P-value'], ascending=True, inplace=True)
data['-Log10-P-value'] = -np.log10(data['P-value'])
data = data.nlargest(10, '-Log10-P-value')

ax = sns.barplot(data=data, x="-Log10-P-value",y='Term',color='blue')


ax.tick_params(axis='y', width = 1, length=8)
ax.tick_params(axis='x', width = 1, length=8)
plt.grid(axis = 'x', color = 'black', linestyle = '--', linewidth = 0.2)
ax.set_axisbelow(True)
plt.xlim(0,4)
plt.xticks([0,1,2,3,4])
plt.xlabel('')
plt.ylabel('')
plt.xticks(visible=False)
plt.yticks(visible=False)

plt.savefig(path+f'GO_biological_processes_{group}_KI_(21A4).png', dpi=150, bbox_inches='tight')
plt.show()  



