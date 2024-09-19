#Z-score analysis
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 23:21:15 2021

@author: congt
"""


import pandas as pd
path = '/path/to/file/'
df = pd.read_csv(path+'mRNA_seq_GAPDH_normalized_WT_1A5_21A4_29A3.bed', sep='\t')
df = df[['Gene', 'WT', '21A4']]

#%%
'''
Calculation of Z-score:
    1. Log10 transform expression value
    2. For each gene:
        2.1. Calculate the mean log10 transformed expression of all samples
        2.2. Calculate the standard deviation of the log10 transformed expression of all samples
        2.3. Calculate z-score in each sample: (Log10-expression - mean(log10-expression))  / std(log10-expression)
        
'''
import numpy as np
df['Log10-WT'] = np.log10(df['WT'])
df['Log10-E518K'] = np.log10(df['21A4'])

df['Mean-log10-expression'] = df[['Log10-WT', 'Log10-E518K']].mean(axis=1)
df['Standard_Deviation'] = df[['Log10-WT', 'Log10-E518K']].std(axis=1)

df = df[df['Standard_Deviation'] > 0]
df['Z-score-WT'] = (df['Log10-WT'] - df['Mean-log10-expression']) / df['Standard_Deviation']
df['Z-score-E518K'] = (df['Log10-E518K'] - df['Mean-log10-expression']) / df['Standard_Deviation']

#%%
'''
plot heatmap for Z-score
'''
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True

data = df[['Z-score-WT', 'Z-score-E518K']]
ax = plt.figure(figsize=(1.5,5))

ax = sns.heatmap(data=data, cmap='coolwarm', vmax=1, vmin=-1,
                 cbar_kws={"shrink": 0.75}, linewidths=0, linecolor='black', cbar=False)

for _, spine in ax.spines.items():
    spine.set_visible(True)
ax.tick_params(axis='y', width = 0, length=0)
ax.tick_params(axis='x', width = 0, length=0)
plt.axvline(x=1, color='black', linestyle='-', zorder=15, linewidth=1)

plt.xlabel('')
plt.ylabel('')
plt.xticks(visible=False, rotation=0)
plt.yticks(visible=False, rotation=0)

# plt.savefig(path+'Z-score-plotting-all-genes-WT-vs-E518K.png', dpi=150, bbox_inches='tight')
plt.show()

#%%
'''
plot Z-score for Ferroptosis gene
'''
df_ferroptosis = pd.read_excel(path+'Ferroptosis/Ferroptosis.xlsx')
df_ferroptosis = df_ferroptosis[['Gene', 'Function']]
df_ferroptosis.drop_duplicates(subset=['Gene'], keep='first', inplace=True)

data = df_ferroptosis.merge(df, on=['Gene'], how='inner')
data.sort_values(['Function'], ascending=True, inplace=True)

data = data[['Gene','Z-score-WT', 'Z-score-E518K']]
data.set_index(['Gene'], inplace=True)
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True


ax = plt.figure(figsize=(1.5,5))

ax = sns.heatmap(data=data, cmap='coolwarm', vmax=1, vmin=-1,
                 cbar_kws={"shrink": 0.75}, linewidths=0.005, linecolor='grey', cbar=False)
for _, spine in ax.spines.items():
    spine.set_visible(True)
ax.tick_params(axis='y', width = 0, length=0)
ax.tick_params(axis='x', width = 0, length=0)
plt.axhline(y=22, color='black', linestyle='-', zorder=15, linewidth=2)
plt.axvline(x=1, color='black', linestyle='-', zorder=15, linewidth=1)

plt.xlabel('')
plt.ylabel('')
plt.xticks(visible=False, rotation=0)
plt.yticks(visible=False, size=4, rotation=0)

plt.savefig(path+'Z-score-plotting-ferroptosis-genes-WT-vs-E518K.png', dpi=150, bbox_inches='tight')
plt.show()




