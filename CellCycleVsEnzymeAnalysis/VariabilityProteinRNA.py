#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
import os
from scipy import stats
import json
import ast
from collections import Counter
from functools import reduce
from sklearn.metrics import r2_score
from statannotations.Annotator import Annotator
from mycolorpy import colorlist as mcp

# %%

### user input
directory = ''
date = ''
hpa = pd.read_csv('20231012_HPAv23_scvmultiloccomp.csv', engine='python') ###HPA v23
enzymes = pd.read_csv('20230111_Human1_ProteinMapping.csv', delimiter=',') ### list of enzymes and pathway annotations from Human1

rnaccd = pd.read_csv('../RNA_CCD_Annotationfrompaper.csv', delimiter=',') ### from Mahdessian D et al, 2021 
protccd = pd.read_csv('../Prot_CCD_Annotationfrompaper.csv', delimiter=',') ### from Mahdessian D et al, 2021 

variance = pd.read_csv("VarianceRNAProtein.csv", delimiter=',') ### from Mahdessian D et al, 2021 

# %%
enzlist = enzymes['ENSG'].unique()
hpa['enzyme'] = False
for i in hpa.index:
    if hpa.at[i, 'Gene'] in enzlist:
        hpa.at[i, 'enzyme'] = True
    
# %%
### now plot the properties of enzymes vs ccd/nonccd proteins
ccdlist = protccd[protccd['CCD_COMP'] == True]['ENSG'].unique()
nonccdlist = protccd[protccd['CCD_COMP'] == False]['ENSG'].unique()

hpalist = hpa['Gene'].unique()
enzlist = hpa[hpa['enzyme'] == True]['Gene'].unique()
scvenzlist = hpa[(hpa['enzyme'] == True) & (hpa['SCV'] == True)]['Gene'].unique()
nonscvenzlist = hpa[(hpa['enzyme'] == True) & (hpa['SCV'] == False)]['Gene'].unique()
nonscvprotlist = hpa[(hpa['SCV'] == False)]['Gene'].unique()
scvprotlist = hpa[(hpa['SCV'] == True)]['Gene'].unique()

nonccdenz = set(enzlist).intersection(nonccdlist)
ccdenz = set(enzlist).intersection(ccdlist)



def boxplot_variability_multiplecond(datalists, ticklabellist, var, rnaorprot):
    data = []
    if rnaorprot == 'rna':
        col = var + '_rna'

    else: 
        col = var + '_cell_prot'
    for dlist in datalists:
        dseries = variance[variance['ENSG'].isin(dlist)][col].dropna()
        data.append(dseries)

    fig = plt.figure(figsize =(7, 7)) 
    # Creating axes instance 
    ax = fig.add_axes([0, 0, 1, 1]) 
    # Creating plot 
    bp = ax.boxplot(data, patch_artist = True, showfliers = False) 
    colors=mcp.gen_color(cmap="viridis",n=len(data))
    #colors = ['#39568CFF', '#73D055FF', '#FDE725FF', 'blue'] 
    for patch, color in zip(bp['boxes'], colors): 
        patch.set_facecolor(color) 
    for median in bp['medians']: 
        median.set(color ='black', linewidth = 1) 

    ax.set_xticklabels(ticklabellist)
    ax.set_ylabel(col)
    title = col  
    plt.title(title)
    plt.ylim(ymin = 0)
    filename_plot = date + '_VariabilityProteinRNAEnzymes_' + col +'.pdf'
    filename_stats = date + '_VariabilityProteinRNAEnzymes_' + col +'_pvaluesKruskal.csv'
    
    pvaldf = pd.DataFrame(columns=['SampleSet'] + ticklabellist)
    i_idx = 0
    for i in data:
        currentdataset = ticklabellist[i_idx]
        pvallist = [currentdataset]
        statlist = [currentdataset]
        for g in data:
            stat, pval = stats.kruskal(i, g)
            pvallist.append(pval)
            statlist.append(stat)
        pvaldf.loc[len(pvaldf)] = pvallist
        i_idx += 1
    pvaldf.to_csv(filename_stats, index=False)

    plt.show
    plt.savefig(filename_plot, bbox_inches="tight")

for var in ['gini', 'cv']:
    ticklist = []
    for protrna in ['rna', 'prot']:
        #boxplot_prot_prop_multiplecond([hpalist, enzlist, ccdlist, nonccdlist], ['AllProteins', 'AllEnzymes', 'CCD Proteins', 'NonCCD Proteins'], cond)
        print('Now processing: ', var, protrna)
        boxplot_variability_multiplecond([hpalist, enzlist, nonccdenz, ccdlist], ['AllFucciProteins', 'AllfucciEnzymes', 'NonCCD enzymes', 'CCD Proteins'], var, protrna)
        #boxplot_variability_multiplecond([hpalist, enzlist, ccdlist, ccdenz, nonccdlist, nonccdenz], ['AllFucciProteins', 'AllFucciEnzymes', 'CCD Proteins', 'CCD Enzymes', 'NonCCD Proteins', 'NonCCD Enzymes'], var, protrna)
# %%
