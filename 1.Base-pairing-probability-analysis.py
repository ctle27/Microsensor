# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 23:21:15 2021

@author: congt
"""

import matplotlib.pyplot as plt
import pandas as pd
path = 'path/to/file/'

df1 = pd.read_csv(path+'count/count_crossing_miR_1A5-rep1.bed',sep='\t')
df2 = pd.read_csv(path+'count/count_crossing_miR_1A5-rep2.bed',sep='\t')
df3 = pd.read_csv(path+'count/count_crossing_miR_1A5-rep3.bed',sep='\t')

df4 = pd.read_csv(path+'count/count_crossing_miR_21A4-rep1.bed',sep='\t')
df5 = pd.read_csv(path+'count/count_crossing_miR_21A4-rep2.bed',sep='\t')
df6 = pd.read_csv(path+'count/count_crossing_miR_21A4-rep3.bed',sep='\t')

df7 = pd.read_csv(path+'count/count_crossing_miR_293T-rep1.bed',sep='\t')
df8 = pd.read_csv(path+'count/count_crossing_miR_293T-rep2.bed',sep='\t')
df9 = pd.read_csv(path+'count/count_crossing_miR_293T-rep3.bed',sep='\t')
#%%remove DROSHA-independent pri-miRNAs
'''
list of DROSHA-independent pri-miRNAs was obtained from Atlas paper (Sup table S1)
'''
df_remove = pd.read_csv(path + 'Pri-miRNADROSHA-dependence_list_(from_Atlas_paper_refs).bed', sep='\t')
remove_list = df_remove['Pri-miRNA'].tolist()

#%%
def count_miRNA_read(df_input, sample):
    df = df_input.copy()
    df = df[df['Sequence'].str.len() >= 18]
    df['Count_hairpin_'+sample] = df.groupby(['hairpin'])['Count'].transform('sum')
    df.drop_duplicates(subset=['hairpin'], keep='first', inplace=True)
    df = df[df['Count_hairpin_'+sample] >= 50] #set minimum read count
    '''
    calculate sum reads of hsa-mir-320a and hsa-mir-320b-1 (DROSHA independent) to normalize
    '''
    # df_control = df[df['hairpin'].isin(['hsa-mir-320a','hsa-mir-320b-1'])]
    df_control = df[df['hairpin'].isin(remove_list)]
    count_control = df_control['Count_hairpin_'+sample].sum()
    df['Normalized_count_'+sample] = df['Count_hairpin_'+sample] / count_control
    df = df[['hairpin','Normalized_count_'+sample]]
    df.reset_index(inplace=True, drop=True)
    return (df)
df1 = count_miRNA_read(df1, '1A5_rep1')
df2 = count_miRNA_read(df2, '1A5_rep2')
df3 = count_miRNA_read(df3, '1A5_rep3')

df4 = count_miRNA_read(df4, '21A4_rep1')
df5 = count_miRNA_read(df5, '21A4_rep2')
df6 = count_miRNA_read(df6, '21A4_rep3')

df7 = count_miRNA_read(df7, '293T_rep1')
df8 = count_miRNA_read(df8, '293T_rep2')
df9 = count_miRNA_read(df9, '293T_rep3')
#%%merge samples
import numpy as np
df_merge = df1.merge(df2, on=['hairpin'], how='outer')
for df in [df3, df4, df5, df6, df7, df8, df9]:
    df_merge = df_merge.merge(df, on=['hairpin'], how='outer')
df_merge.fillna(0, inplace=True)

df_merge['Mean_1A5'] = df_merge[['Normalized_count_1A5_rep1', 'Normalized_count_1A5_rep2', 'Normalized_count_1A5_rep3']].mean(axis=1)
df_merge['Mean_21A4'] = df_merge[['Normalized_count_21A4_rep1', 'Normalized_count_21A4_rep2', 'Normalized_count_21A4_rep3']].mean(axis=1)
df_merge['Mean_293T'] = df_merge[['Normalized_count_293T_rep1', 'Normalized_count_293T_rep2', 'Normalized_count_293T_rep3']].mean(axis=1)

df_merge = df_merge[(df_merge['Mean_21A4'] != 0) & (df_merge['Mean_293T'] != 0)]
df_merge['293T_over_1A5'] = np.log2(df_merge['Mean_293T'] + 0.001) - np.log2(df_merge['Mean_1A5'] + 0.001)
df_merge['21A4_over_1A5'] = np.log2(df_merge['Mean_21A4'] + 0.001) - np.log2(df_merge['Mean_1A5'] + 0.001)
df_merge['293T_over_21A4'] = np.log2(df_merge['Mean_293T'] + 0.001) - np.log2(df_merge['Mean_21A4'] + 0.001)

df_merge.to_csv('HEK293TWT-DGCR8KO-DGCR8E518K-smallRNA-normalized-expression.tsv', sep='\t', index=False)

#%%remove DROSHA-independent pri-miRNAs
'''
list of DROSHA-independent pri-miRNAs was obtained from Atlas paper (Sup table S1)
'''
df_remove = pd.read_csv(path + 'Pri-miRNADROSHA-dependence_list_(from_Atlas_paper_refs).bed', sep='\t')
remove_list = df_remove['Pri-miRNA'].tolist()

df_merge = df_merge[~df_merge['hairpin'].isin(remove_list)]

#%%obtain sequence of pri-miRNAs
df_sequence = pd.read_csv(path + 'miRGeneDB_miRBase_primiRNA_miRNA_sequences.bed', sep='\t')
df_sequence = df_sequence[['MiRBase_ID', '30ntExt_premiRNA']]

df_merge = df_sequence.merge(df_merge, left_on='MiRBase_ID', right_on='hairpin', how='inner')
df_merge = df_merge.sort_values(by=['293T_over_21A4'], ascending = True)
df_merge = df_merge.drop_duplicates(subset=['MiRBase_ID'], keep='first')
df_merge['Pri_length'] = df_merge['30ntExt_premiRNA'].apply(len)
df_merge.reset_index(inplace=True, drop=True)

#%%
#pearson correlation analysis
import numpy as np
data = df_merge
sample = '21A4'
data1 = data['Normalized_count_'+sample+'_rep1'].to_numpy()
data2 = data['Normalized_count_'+sample+'_rep3'].to_numpy()
data1[np.isnan(data1)] = 0
data2[np.isnan(data2)] = 0
r = np.corrcoef(data1, data2)
print (r)


#%%
import seaborn as sns
ax = sns.ecdfplot(data=df_merge, x="293T_over_1A5",color='green')
ax = sns.ecdfplot(data=df_merge, x="21A4_over_1A5",color='red')

plt.xlim(-4,12)
plt.ylim(0,1)

plt.xlabel('Expression fold over pXG (log2)')
plt.ylabel('Cummulative fraction')

plt.xlabel('')
plt.ylabel('')

ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])

plt.tick_params(direction='out', length=5, width=1, colors='black')
plt.grid(color = 'black', linestyle = '--', linewidth = 0.4)
plt.savefig(path+'DGCR8_WT_E518K_cummulative_plot.png',dpi=300,bbox_inches='tight')
plt.show()


#%%
df_most = df_merge[-30:]
df_least = df_merge[:30]
import seaborn as sns

df_plot = df_most
title = 'most_affected'

ax = sns.ecdfplot(data=df_plot, x="293T_over_1A5",color='green')
ax = sns.ecdfplot(data=df_plot, x="21A4_over_1A5",color='red')

plt.xlim(-4,12)
plt.ylim(0,1)

plt.xlabel('Expression fold over pXG (log2)')
plt.ylabel('Cummulative fraction')

plt.xlabel('')
plt.ylabel('')

ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])

plt.tick_params(direction='out', length=5, width=1, colors='black')
plt.grid(color = 'black', linestyle = '--', linewidth = 0.4)
plt.savefig(path+f'DGCR8_WT_E518K_cummulative_plot_{title}_pri.png',dpi=300,bbox_inches='tight')
plt.show()

#%%         base-pairing probability of most and least affected pri-miRNAs
df_most = df_merge[-30:]
df_least = df_merge[:30]

#1 output fasta files for all pri-miRNA detected
f = open(path+'detected_pri-miRNA.fa','w+')
for i,val in enumerate(df_merge['MiRBase_ID']):
    f.write('>'+val+'\n')
    f.write(df_merge['30ntExt_premiRNA'][i]+'\n')
f.close()
#2 RNAfold -p input.fasta --> obtain files with bp probability

import pandas as pd
df = pd.DataFrame()
for i in range (150):
    df.loc[i,'Cummulative_prob'] = 0

#%%
from os import listdir
from os.path import isfile, join
def bp_probability(condition_list,len_list):
    df = pd.DataFrame()
    for i in range (150):
        df.loc[i,'Cummulative_prob_5p'] = 0
        df.loc[i,'Cummulative_prob_3p'] = 0
    path2 = 'path/to/files/'
    files_in_dir = [ f for f in listdir(path2) if isfile(join(path2,f)) ]
    for file in files_in_dir:
        check_name = file.strip('_dp.ps')
        if check_name in condition_list:
            pri_len = len_list[condition_list.index(check_name)]
            with open(path2+file, 'r') as f:
                for l in f:
                    line_split = []
                    if 'ubox' in l:
                        if '%' not in l:
                            if '/' not in l:
                                line_split = l.split(' ')
                                pos_5p = int(line_split[0]) - 1
                                pos_3p = pri_len - pos_5p - 1
                                prob1 = float(line_split[2])*float(line_split[2])
                                df.loc[pos_5p,'Cummulative_prob_5p'] += prob1
                                df.loc[pos_3p,'Cummulative_prob_3p'] += prob1
                                
                                pos_5p = int(line_split[1]) - 1
                                pos_3p = pri_len - pos_5p - 1
                                prob2 = float(line_split[2])*float(line_split[2])
                                df.loc[pos_5p,'Cummulative_prob_5p'] += prob2
                                df.loc[pos_3p,'Cummulative_prob_3p'] += prob2
    df.reset_index(inplace=True)
    for i,val in enumerate(df['index']):
        if val < 30:
            df.loc[i,'Position'] = val - 30
        if val >= 30:
            df.loc[i,'Position'] = val - 29

    df = df[df['Position'] < 26]
    df = df[df['Position'] > -11]
    return df
#%%
#list of all pri-miRNA detected
allpri_list = df_merge['MiRBase_ID'].tolist()
allpri_length_list = df_merge['Pri_length'].tolist()
#list of most affected pri-miRNAs
top_list = df_most['MiRBase_ID'].tolist()
top_length_list = df_most['Pri_length'].tolist()
#list of least affected pri-miRNAs
bot_list = df_least['MiRBase_ID'].tolist()
bot_length_list = df_least['Pri_length'].tolist()

#%%
df_all = bp_probability(allpri_list,allpri_length_list)
df_all['Cummulative_prob_5p'] = df_all['Cummulative_prob_5p']/len(allpri_list)
df_all['Cummulative_prob_3p'] = df_all['Cummulative_prob_3p']/len(allpri_list)
  
df_most = bp_probability(top_list,top_length_list)
df_most['Cummulative_prob_5p'] = df_most['Cummulative_prob_5p']/len(top_list)   
df_most['Cummulative_prob_3p'] = df_most['Cummulative_prob_3p']/len(top_list)  

df_least = bp_probability(bot_list,bot_length_list)
df_least['Cummulative_prob_5p'] = df_least['Cummulative_prob_5p']/len(bot_list)
df_least['Cummulative_prob_3p'] = df_least['Cummulative_prob_3p']/len(bot_list)

#%%bp-ring probability at 5p strand
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

ax = plt.figure(figsize=(10,2))
ax = sns.lineplot(x='Position',y='Cummulative_prob_5p',data=df_most,color='red',marker="o")
ax = sns.lineplot(x='Position',y='Cummulative_prob_5p',data=df_least,color='green',marker="o")
ax = sns.lineplot(x='Position',y='Cummulative_prob_5p',data=df_all,color='purple',marker="o")

plt.tick_params(direction='out', length=5, width=1, colors='black')
plt.grid(color = 'black', linestyle = '--', linewidth = 0.3)
plt.ylim(-0.1,1.1)
plt.yticks([0,0.25,0.5,0.75,1])
plt.xlim(-10.5,25.5)
xticks = [i for i in range(-10, 26)]
plt.xticks(xticks)
plt.gca().xaxis.set_major_locator(MultipleLocator(5))
plt.gca().xaxis.set_minor_locator(AutoMinorLocator(5))

ax.set(xticklabels=[])
ax.set(yticklabels=[])
ax.set(xlabel=None)
ax.set(ylabel=None)
plt.savefig(path+'bp-ing_probability_5p.png',dpi=300,bbox_inches='tight')
plt.show()

#%%bp-ring probability at 3p strand
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

ax = plt.figure(figsize=(10,2))
ax = sns.lineplot(x='Position',y='Cummulative_prob_3p',data=df_most,color='red',marker="o")
ax = sns.lineplot(x='Position',y='Cummulative_prob_3p',data=df_least,color='green',marker="o")
ax = sns.lineplot(x='Position',y='Cummulative_prob_3p',data=df_all,color='purple',marker="o")

plt.tick_params(direction='out', length=5, width=1, colors='black')
plt.grid(color = 'black', linestyle = '--', linewidth = 0.3)
plt.ylim(-0.1,1.1)
plt.yticks([0,0.25,0.5,0.75,1])
plt.xlim(-10.5,25.5)
xticks = [i for i in range(-10, 26)]
plt.xticks(xticks)
plt.gca().xaxis.set_major_locator(MultipleLocator(5))
plt.gca().xaxis.set_minor_locator(AutoMinorLocator(5))

ax.set(xticklabels=[])
ax.set(yticklabels=[])
ax.set(xlabel=None)
ax.set(ylabel=None)
plt.savefig(path+'bp-ing_probability_3p.png',dpi=300,bbox_inches='tight')
plt.show()








