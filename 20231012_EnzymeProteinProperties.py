#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from mycolorpy import colorlist as mcp
from statannotations.Annotator import Annotator


#%%
date = ''
pathway = pd.read_csv("EnzymeLists/20230111_Human1_ProteinMapping.csv", engine='python') ## Human1 export of all enzymes and pathway annotations
properties = pd.read_csv("../EnzymeLists/PropertiesTable.csv", engine='python')
hpa = pd.read_csv("EnzymeLists/20231012_HPAv23_scvmultiloccomp.csv", engine='python') ### HPA v23
hpa = hpa.rename(columns={'multilocalizing': 'Multilocalizing'})
fucci = pd.read_csv("Prot_CCD_Annotationfrompaper.csv", engine='python') ### results from Mahdessian D et al, 2021
#hpa['SCV'] = (hpa['SCV'].astype(str)).str.upper()

###add column "enzyme" to hpa dataframe
hpa['enzyme'] = 'False'
for i in hpa.index:
    if hpa.Gene_name[i] in pathway['Gene_name'].tolist():
       hpa.enzyme[i] = 'True'
    else:
        hpa.enzyme[i] = 'False'


#%%
### merge enzyme/hpa data with properties
enzprop = pd.merge(pathway, properties, left_on='Gene_name',right_on='Protein', how='inner')
hpaprop = pd.merge(hpa, properties, left_on='Gene_name',right_on='Protein', how='inner')

hpaprop['Multilocalizing'] = hpaprop['Multilocalizing'].map({True: 'True', False: 'False'})
hpaprop['SCV'] = hpaprop['SCV'].map({True: 'True', False: 'False'})
enzprop['Multilocalizing'] = enzprop['Multilocalizing'].map({True: 'True', False: 'False'})
enzprop['SCV'] = enzprop['SCV'].map({True: 'True', False: 'False'})

enzhpaprop = hpaprop[hpaprop['enzyme'] == 'True']

#%%
### function to plot the protein properties in a boxplot
def propertiesboxplot(df, condition, testing, dataset):
    plt.figure()
    sns.set_style("whitegrid")
    sns.boxplot(x = condition, y = testing, data = df)
    plt.show
    filename = date + '_' + dataset + '_' + condition + '_' + testing + '.pdf'
    plt.show
    plt.savefig('../ProteinProperties/' + filename)
    plt.close

def propertiesboxplot_with_test(df, condition, testing, dataset, cond_list, stat_boxpairs):   
    #df[condition] = df[condition].map({True: 'True', False: 'False'})
    fig, ax = plt.subplots(figsize=(5, 5))
    #hue_order=['0uM', '1uM', '10uM']
    order = cond_list
    bp = sns.boxplot(y=testing, x=condition, data=df, palette="colorblind", order=order, showfliers = False)
    title = testing + '_' + dataset
    plt.title(title)

    annotator = Annotator(ax=bp, pairs=stat_boxpairs,y=testing, x=condition, data=df, order=order)
    annotator.configure(test='Mann-Whitney', verbose=True, text_format='full').apply_and_annotate()
    filename = date + '_' + dataset + '_' + condition + '_' + testing + '.pdf'
    plt.savefig('../ProteinProperties/' + filename)
    plt.close
    
#%%
pvalues = pd.DataFrame(columns = ['Dataset','Condition', 'ProteinProperty', 'Pvalue', 'comment'])
i = 0
for df in [hpaprop, enzhpaprop]:
    if i == 0:
        dataset = 'hpa_proteins'
    else:
        dataset='Enzymes'
    for cond in ['Multilocalizing', 'SCV', 'compartment']:
        for test in ['Melting Temperature (deg C)', 'Log2 Length', 'Disordered Residue Fraction']:
            if cond != 'compartment':
                propertiesboxplot_with_test(df,cond,test, dataset, ['False', 'True'], [('False', 'True')])
            if cond == 'compartment':
                propertiesboxplot_with_test(df,cond,test, dataset, ['cell', 'cyto', 'nuc'], [('cell', 'cyto'), ('cell', 'nuc'), ('cyto', 'nuc')])
                nucdf = df[df[cond]=='nuc'][test].dropna()
                cytodf = df[df[cond]=='cyto'][test].dropna()
                celldf = df[df[cond]=='cell'][test].dropna()
                pval = stats.ttest_ind(nucdf,cytodf)
                #pval = stats.kruskal(nucdf,cytodf)
                pvalues = pvalues.append({'Dataset':dataset,'Condition':cond, 'ProteinProperty':test,
                                        'Pvalue':pval, 'comment':'nuc_vs_cyto'}, ignore_index=True)
                pval = stats.ttest_ind(nucdf,celldf)
                #pval = stats.kruskal(nucdf,celldf)
                pvalues = pvalues.append({'Dataset':dataset,'Condition':cond, 'ProteinProperty':test,
                                        'Pvalue':pval, 'comment':'nuc_vs_cell'}, ignore_index=True)
                pval = stats.ttest_ind(cytodf,celldf)
                #pval = stats.kruskal(cytodf,celldf)
                pvalues = pvalues.append({'Dataset':dataset,'Condition':cond, 'ProteinProperty':test,
                                        'Pvalue':pval, 'comment':'cyto_vs_cell'}, ignore_index=True)
            else:      
                truedf = df[df[cond]==True][test].dropna()
                falsedf = df[df[cond]==False][test].dropna()
                pval = stats.ttest_ind(truedf,falsedf)
                #pval = stats.kruskal(truedf,falsedf)
                pvalues = pvalues.append({'Dataset':dataset,'Condition':cond, 'ProteinProperty':test,
                                        'Pvalue':pval, 'comment':''}, ignore_index=True)
    if dataset == 'hpa_proteins':
        for test in ['Melting Temperature (deg C)', 'Log2 Length', 'Disordered Residue Fraction']:
            propertiesboxplot_with_test(df,'enzyme',test, dataset, ['False', 'True'], [('False', 'True')])
            truedf = df[df['enzyme']=='True'][test].dropna()
            falsedf = df[df['enzyme']=='False'][test].dropna()
            pval = stats.ttest_ind(truedf,falsedf)
            #pval = stats.kruskal(truedf,falsedf)
            pvalues = pvalues.append({'Dataset':dataset,'Condition':'enzyme', 'ProteinProperty':test,
                                    'Pvalue':pval, 'comment':''}, ignore_index=True)
    i+=1
pvalues.to_csv('../ProteinProperties/' + date + '_ProteinPropertiesPvalues_HPAandEnzymes.csv', index=False)

#%%
### now plot the properties of enzymes vs ccd/nonccd proteins
ccdlist = fucci[fucci['CCD_COMP'] == True]['ENSG'].unique()
nonccdlist = fucci[fucci['CCD_COMP'] == False]['ENSG'].unique()
hpalist = hpa['Gene'].unique()
enzlist = hpa[hpa['enzyme'] == 'True']['Gene'].unique()
scvenzlist = hpa[(hpa['enzyme'] == 'True') & (hpa['SCV'] == True)]['Gene'].unique()
nonscvenzlist = hpa[(hpa['enzyme'] == 'True') & (hpa['SCV'] == False)]['Gene'].unique()
nonscvprotlist = hpa[(hpa['SCV'] == False)]['Gene'].unique()
scvprotlist = hpa[(hpa['SCV'] == True)]['Gene'].unique()

def boxplot_prot_prop_multiplecond(datalists, ticklabellist, cond, stat_boxpairs):
    df_list = []
    i = 0
    for dlist in datalists:
        df = hpaprop[hpaprop['Gene'].isin(dlist)]
        df['label'] = ticklabellist[i]
        df_list.append(df)
        i += 1
    data = pd.concat(df_list)
    data = data[data[cond].notna()]
    print(type(data))
    fig, ax = plt.subplots(figsize=(5, 5))
    bp = sns.boxplot(y=cond, x='label', data=data, palette="colorblind", order=ticklabellist, showfliers = False)
    title = cond + '_' + dataset
    plt.title(title)
    annotator = Annotator(ax=bp, pairs=stat_boxpairs,y=cond, x='label', data=data, order=ticklabellist)
    annotator.configure(test='Mann-Whitney', verbose=True, text_format='full').apply_and_annotate()
    filename = date + '_ProteinProperties_' + cond +'.pdf'
    #plt.show
    plt.savefig('../ProteinProperties/' +filename, bbox_inches="tight")



for cond in ['Disordered Residue Fraction', 'Melting Temperature (deg C)', 'Log2 Length']:
    #boxplot_prot_prop_multiplecond([hpalist, enzlist, ccdlist, nonccdlist], ['AllProteins', 'AllEnzymes', 'CCD Proteins', 'NonCCD Proteins'], cond)
    stat_boxpairs = [('AllProteins', 'AllEnzymes'), ('AllProteins', 'SCV Enzymes'), ('AllProteins', 'CCD Proteins'), 
                        ('AllEnzymes', 'SCV Enzymes'), ('AllEnzymes', 'CCD Proteins'), ('SCV Enzymes', 'CCD Proteins')]
    check = boxplot_prot_prop_multiplecond([hpalist, enzlist, scvenzlist, ccdlist], ['AllProteins', 'AllEnzymes', 'SCV Enzymes', 'CCD Proteins'], cond, stat_boxpairs)
    #boxplot_prot_prop_multiplecond([hpalist, enzlist, scvprotlist, scvenzlist, nonscvprotlist, nonscvenzlist, ccdlist], ['AllProteins', 'AllEnzymes', 'SCV Proteins', 'SCV Enzymes', 'NonSCV Proteins', 'NonSCV Enzymes', 'CCD Proteins'], cond)

# %%
