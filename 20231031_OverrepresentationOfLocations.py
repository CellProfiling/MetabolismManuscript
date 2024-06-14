### author christian Gnann
### this script contains function to calculate the overrepresentation of locations across different pathways groups and bitweeen different sets of proteins (e.g. scv vs stably expressed)
### it generates dataframes with the significance calcuated from binomial test (BH corrected) as well as plots the results in heatmaps and/or barplots

#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from functools import reduce
from collections import Counter
import seaborn as sns
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from statsmodels.stats.multitest import multipletests


#%%
### load files 
date = ''
hpa = pd.read_csv("20231012_HPAv23_scvmultiloccomp.csv", engine='python') ### HPA v23
pathway = pd.read_csv("20230111_Human1_ProteinMapping.csv", engine='python') ## file with all enzymes and pathways from Human1
enz_hpa = hpa[hpa['Gene'].isin(pathway['ENSG'].unique())]

hpaproteins = hpa['Gene_name'].unique().tolist()

#%% stuff needed further down --> grouped location into structures and compartments
locations = ['Actin filaments', 'Aggresome', 'Cell Junctions', 'Centriolar satellite', 'Centrosome',
                'Cleavage furrow', 'Cytokinetic bridge', 'Cytoplasmic bodies', 'Cytosol', 
                'Endoplasmic reticulum', 'Endosomes', 'Focal adhesion sites', 'Golgi apparatus', 
                'Intermediate filaments', 'Kinetochore', 'Lipid droplet', 'Lysosomes', 'Microtubule ends',
                'Microtubules', 'Midbody', 'Midbody ring', 'Mitochondria', 'Mitotic chromosome', 'Mitotic spindle',
                'Nuclear bodies', 'Nuclear membrane', 'Nuclear speckles', 'Nucleoli', 'Nucleoli fibrillar center',
                'Nucleoli rim', 'Nucleoplasm','Peroxisomes', 'Plasma membrane', 'Rods & Rings', 'Vesicles']

csk_comp = ['Actin filaments', 'Centriolar satellite', 'Centrosome',
                'Cleavage furrow', 'Cytokinetic bridge', 'Focal adhesion sites', 'Intermediate filaments', 'Microtubule ends',
                'Microtubules', 'Midbody', 'Midbody ring', 'Mitotic spindle']

nuc_comp = ['Kinetochore', 'Mitotic chromosome', 'Nuclear bodies', 'Nuclear membrane', 'Nuclear speckles', 'Nucleoli', 'Nucleoli fibrillar center',
                'Nucleoli rim', 'Nucleoplasm']

sec_comp = ['Cell Junctions', 'Endoplasmic reticulum', 'Endosomes', 'Golgi apparatus', 'Lipid droplet', 'Lysosomes', 
            'Peroxisomes', 'Plasma membrane', 'Vesicles']

cyt_comp = ['Actin filaments', 'Aggresome', 'Cell Junctions', 'Cytoplasmic bodies', 'Cytosol', 
                'Endoplasmic reticulum', 'Endosomes', 'Focal adhesion sites', 'Golgi apparatus', 
                'Intermediate filaments', 'Lipid droplet', 'Lysosomes', 'Microtubule ends',
                'Microtubules', 'Mitochondria', 'Peroxisomes', 'Plasma membrane', 'Rods & Rings', 'Vesicles']
mitosis = ['Cleavage furrow', 'Cytokinetic bridge','Kinetochore', 'Midbody', 'Midbody ring', 
                'Mitotic chromosome', 'Mitotic spindle']
non_mitosis = list((Counter(locations) - Counter(mitosis)).elements())

groupedlocations = ['Actin filaments', 'Centrosome', 'Cytoplasmic Aggregates', 'Cytosol', 'Endoplasmic reticulum', 'Golgi apparatus', 
                    'Intermediate filaments', 'Microtubules', 'Mitochondria',
                    'Nuclear bodies', 'Nuclear membrane', 'Nuclear speckles', 'Nucleoli', 
                    'Nucleoplasm', 'Plasma membrane', 'Vesicles']

#%%

def loc_count(df, dataset, loclist, allorgrouped):
    '''
    this function counts the number of genes in a certain location for a subset of proteins
    input: df --> subsetted hpa dataframe
    dataset: string describing the dataset; only used for printing
    '''
    loc_count_df = pd.DataFrame(columns = ['Localization', 'Count', 'WhichList'])
    for location in loclist:
        #print('Now processing: ',location)
        if allorgrouped == 'all':
            locseries = df['All_Locations']
        elif allorgrouped == 'grouped':
            locseries = df['All_GroupedLocations']
        count = int((locseries.str.count(location)).sum())
        #print('The count is: ', count)
        loc_count_df = loc_count_df.append({'Localization': location, 'Count': count,
                                            'WhichList': dataset}, ignore_index=True)
    return loc_count_df

hpa_loc_count = loc_count(hpa, 'all hpa', locations, 'all')
hpa_loc_count_grouped = loc_count(hpa, 'all hpa - grouped', groupedlocations, 'grouped')
enz_loc_count_grouped = loc_count(enz_hpa, 'all enz - grouped', groupedlocations, 'grouped')
multilocenz_loc_count_grouped = loc_count(enz_hpa[enz_hpa['multilocalizing']==True], 'multilocalizing enz - grouped', groupedlocations, 'grouped')
scvenz_loc_count_grouped = loc_count(enz_hpa[enz_hpa['SCV']==True], 'multilocalizing enz - grouped', groupedlocations, 'grouped')

#%%
def LocationOverrep(LocCount_bg , bg, pathway_df, loclist, allorgrouped):
    '''
    this function ccalculates whether a protein is enriched in certain stuctures compared to a background list
    input: Loc count df (output from loc_count function; subsetted )
    dataset: string describing the dataset; only used for printing
    it will do so in a for loop and go through all metabolic pathways
    '''    
    loc_overrep_df = pd.DataFrame(columns = ['Localization', 'Count', 'percentage', 'enrichment', 'pvalue_underrepresentation', 
                                            'pvalue_overrepresentation', 'background_probability', 'pathway'])
    len_bg = len(bg)
    pathway_list = pathway['PathwayGroup'].unique()
        
    for path in pathway_list:
        print(path)
        test_list = pathway_df[pathway_df['PathwayGroup'] == path]['ENSG'].unique()
        test = bg[bg['Gene'].isin(test_list)]    # subset hpa with all genes for pathway
        len_test = len(test)
        print(bg.shape, test.shape)
        
        LocCount_test = loc_count(test, path, loclist, allorgrouped)
        #print(LocCount_test.head())
        for i in LocCount_test.index:
            location = LocCount_test.Localization[i]
            count = int(LocCount_test.Count[i])
            test_probability = count/len_test
            print(location, ': ', count)
            
            bg_count = (LocCount_bg.loc[(LocCount_bg.Localization == location)].Count).item()
            bg_probability = bg_count/len_bg
        
            pval_under = stats.binom_test(count, n=len_test, p=bg_probability, alternative='less')
            pval_over = stats.binom_test(count, n=len_test, p=bg_probability, alternative='greater')
            print(test_probability, pval_under, pval_over)
            ratio = test_probability/bg_probability
            
            if ratio <= 1:
                if test_probability != 0:
                    ratio = bg_probability/test_probability*(-1)
                else:
                    ratio = 0
                print(ratio)
            
            loc_overrep_df = loc_overrep_df.append({'Localization': location, 'Count': count, 'percentage':test_probability, 'enrichment':ratio,
                                                    'pvalue_underrepresentation': pval_under, 'pvalue_overrepresentation':pval_over,
                                                    'background_probability':bg_probability, 'pathway':path, }, ignore_index=True)
    ### Adjust p values
    reject, p_values_corr_under, _, _ = multipletests(loc_overrep_df['pvalue_underrepresentation'], alpha=0.5, method='fdr_bh')
    loc_overrep_df['p_val_under_bonf'] = p_values_corr_under
    reject, p_values_corr_over, _, _ = multipletests(loc_overrep_df['pvalue_overrepresentation'], alpha=0.5, method='fdr_bh')
    # Add the corrected p-values to the DataFrame
    loc_overrep_df['p_val_over_bonf'] = p_values_corr_over
    return loc_overrep_df

overrep_df = LocationOverrep(hpa_loc_count, hpa, pathway, locations, 'all')

overrep_df_grouped = LocationOverrep(hpa_loc_count_grouped, hpa, pathway, groupedlocations, 'grouped')

enzoverrep_df_grouped = LocationOverrep(enz_loc_count_grouped, enz_hpa, pathway, groupedlocations, 'grouped')

# %%
### heatmap of manually selected locations

heatmap_loc = ['Mitochondria', 'Cytosol', 'Endoplasmic reticulum', 'Golgi apparatus', 'Vesicles', 'Plasma membrane',
                    'Nuclear membrane', 'Nucleoplasm', 'Nucleoli', 'Nuclear speckles', 'Nuclear bodies', 'Actin filaments', 
                    'Intermediate filaments', 'Centrosome'  ]

overrep_df_heatmap = overrep_df_grouped[overrep_df_grouped['Localization'].isin(heatmap_loc)]

### --> set the enrichment for the non significant ones to 0
for i in overrep_df_heatmap.index:
    if overrep_df_heatmap.p_val_under_bonf[i] >= 0.05 and overrep_df_heatmap.p_val_over_bonf[i] >= 0.05:
        overrep_df_heatmap.enrichment[i] = 0

overrep_df2 = overrep_df_heatmap.pivot("pathway", "Localization", "enrichment")
#ax = sns.heatmap(overrep_df2, cmap="vlag")
overrep_df2 = overrep_df2[heatmap_loc]
sns.heatmap(overrep_df2, cmap="seismic", vmin = -10, vmax = 10)
plt.savefig('20231031_Heatmap_GroupedLocationOverrep_perPathway.pdf', bbox_inches='tight')


# %%
### now let's see whether there is something interesting when looking at enzymes only as the background

heatmap_loc = ['Mitochondria', 'Cytosol', 'Endoplasmic reticulum', 'Golgi apparatus', 'Vesicles', 'Plasma membrane',
                    'Nuclear membrane', 'Nucleoplasm', 'Nucleoli', 'Nuclear speckles', 'Nuclear bodies', 'Actin filaments', 
                    'Intermediate filaments', 'Centrosome']

overrep_df_heatmap = enzoverrep_df_grouped[enzoverrep_df_grouped['Localization'].isin(heatmap_loc)]

### --> set the enrichment for the non significant ones to 0
for i in overrep_df_heatmap.index:
    if overrep_df_heatmap.p_val_under_bonf[i] >= 0.05 and overrep_df_heatmap.p_val_over_bonf[i] >= 0.05:
        overrep_df_heatmap.enrichment[i] = 0
    if overrep_df_heatmap.pvalue_underrepresentation[i] >= 0.05 and overrep_df_heatmap.pvalue_overrepresentation[i] >= 0.05:
        overrep_df_heatmap.enrichment[i] = 0

overrep_df2 = overrep_df_heatmap.pivot("pathway", "Localization", "enrichment")
#ax = sns.heatmap(overrep_df2, cmap="vlag")
overrep_df2 = overrep_df2[heatmap_loc]
sns.heatmap(overrep_df2, cmap="seismic", vmin = -11, vmax = 11)
plt.savefig('20231031_Heatmap_GroupedLocationOverrepVsEnzymes_perPathway_bhcorrected.pdf', bbox_inches='tight')
plt.show()
# %%
### now let's look into big scale location overrepresentation (scv vs nonscv; enz vs non enz, ...)
def LocationOverrep_nosubsetting(LocCount_test, test, LocCount_bg , bg):
    loc_overrep_df = pd.DataFrame(columns = ['Localization', 'Count', 'percentage', 'enrichment', 'pvalue_underrepresentation', 
                                            'pvalue_overrepresentation', 'background_probability'])
    len_test = len(test)
    len_bg = len(bg)
    #print(len_test, len_background)
    for i in LocCount_test.index:
        location = LocCount_test.Localization[i]
        count = int(LocCount_test.Count[i])
        test_probability = count/len_test
        print(location, ': ', count)
        bg_count = (LocCount_bg.loc[(LocCount_bg.Localization == location)].Count).item()
        bg_probability = bg_count/len_bg
        pval_under = stats.binom_test(count, n=len_test, p=bg_probability, alternative='less')
        pval_over = stats.binom_test(count, n=len_test, p=bg_probability, alternative='greater')
        print(test_probability, pval_under, pval_over)
        ratio = test_probability/bg_probability
        if ratio <= 1:
            ratio = (bg_probability/test_probability)*-1
        
        loc_overrep_df = loc_overrep_df.append({'Localization': location, 'Count': count, 'percentage':test_probability,'enrichment':ratio,
                                                'pvalue_underrepresentation': pval_under, 'pvalue_overrepresentation':pval_over,
                                                'background_probability':bg_probability, }, ignore_index=True)
    
    reject, p_values_corr_under, _, _ = multipletests(loc_overrep_df['pvalue_underrepresentation'], alpha=0.5, method='fdr_bh')
    loc_overrep_df['p_val_under_bonf'] = p_values_corr_under
    reject, p_values_corr_over, _, _ = multipletests(loc_overrep_df['pvalue_overrepresentation'], alpha=0.5, method='fdr_bh')
    loc_overrep_df['p_val_over_bonf'] = p_values_corr_over
    return loc_overrep_df

enz_loc_overrep = LocationOverrep_nosubsetting(enz_loc_count_grouped, enz_hpa, hpa_loc_count_grouped, hpa)
multilocenz_loc_overrep = LocationOverrep_nosubsetting(multilocenz_loc_count_grouped, enz_hpa[enz_hpa['multilocalizing']==True], enz_loc_count_grouped, enz_hpa)
scvenz_loc_overrep = LocationOverrep_nosubsetting(scvenz_loc_count_grouped, enz_hpa[enz_hpa['SCV']==True], enz_loc_count_grouped, enz_hpa)

enz_loc_overrep.to_csv(date + '_LocationOverrepresentation_Enzymes.csv', index=False)
multilocenz_loc_overrep.to_csv(date + '_LocationOverrepresentation_MulitlocEnzymes.csv', index=False)
scvenz_loc_overrep.to_csv(date + '_LocationOverrepresentation_ScvEnzymes.csv', index=False)
#%%
def barplot_loc_overrep_enzymes(df, dataset, bh):
    df = df[df['Localization'] != 'Centrosome']
    df = df.sort_values(by=['enrichment']).reset_index()
    #df.plot(kind='bar', stacked=True, color=['orange', 'lightgrey'])
    
    objects = df['Localization'].to_list()
    y_pos = np.arange(len(objects))
    enrichment = df['enrichment'].to_list()
    for i in range(len(enrichment)):
        if enrichment[i] >= 1:
            enrichment[i] -= 1
        else:
            enrichment[i] += 1
    if bh == True:
        pvalsunder = df['p_val_under_bonf'].to_list()
        pvalsover = df['p_val_over_bonf'].to_list()
    else:
        pvalsunder = df['pvalue_underrepresentation'].to_list()
        pvalsover = df['pvalue_overrepresentation'].to_list()
    bluecmap = mpl.colors.LinearSegmentedColormap.from_list('custom blue', [mpl.colors.to_hex(plt.cm.get_cmap('seismic')(0)),mpl.colors.to_hex(plt.cm.get_cmap('seismic')(0.47))], N=256)
    redcmap = plt.cm.get_cmap('Reds_r')
    redcmap = mpl.colors.LinearSegmentedColormap.from_list('custom red', [mpl.colors.to_hex(plt.cm.get_cmap('seismic')(0.99999)),mpl.colors.to_hex(plt.cm.get_cmap('seismic')(0.53))], N=256)
    maxpvalup = min(pvalsover)
    maxpvaldown = min(pvalsunder)
    maxpval = min([maxpvalup, maxpvaldown])
    colorlist = []
    for location in objects:
        i = objects.index(location)
        if enrichment[i] <= 0:
            if pvalsunder[i] <= 0.05:
                colorvalue = abs(np.log10(pvalsunder[i]))
                color_cmap = abs(colorvalue-abs(np.log10(maxpval)))
                color_cmap = (abs(np.log10(maxpval)) - colorvalue)/abs(np.log10(maxpval))
                colorlist.append(bluecmap(color_cmap)) #
            else:
                colorlist.append('white')
        else:
            if pvalsover[i] <= 0.05:
                colorvalue = abs(np.log10(pvalsover[i]))
                color_cmap = (abs(np.log10(maxpval)) - colorvalue)/abs(np.log10(maxpval))
                colorlist.append(redcmap(color_cmap))
            else:
                colorlist.append('white')

    plt.bar(y_pos, enrichment, align='center', alpha=1, color = colorlist, edgecolor = "black")
    ax = plt.gca()
    # remove the existing ticklabels
    ax.set_xticklabels([])
    # remove the extra tick on the negative bar
    ax.set_xticks([idx for (idx, x) in enumerate(enrichment) if x < 0])
    ax.spines["bottom"].set_position(("data", 0))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    label_offset = 0.5
    for location, (x_position, y_position) in zip(objects, enumerate(enrichment)):
        if y_position > 0:
            label_y = -label_offset
        else:
            label_y = y_position - label_offset
        ax.text(x_position, label_y, location, ha="center", va="top",rotation='vertical')
    ax.set_title('Enrichment of locations ' + dataset)
    plt.ylabel('fold enrichment')
    divider = make_axes_locatable(plt.gca())
    ax_cb = divider.new_horizontal(size="5%", pad=0.05)    
    cb1 = mpl.colorbar.ColorbarBase(ax_cb, cmap=mpl.cm.seismic, orientation='vertical', ticks=[0, 0.5, 1])
    cbar_label = '-' + str(round(np.log10(maxpval)))
    cb1.ax.set_yticklabels([cbar_label, '0', cbar_label])
    cb1.ax.set_ylabel('log10 p value binom test')
    plt.gcf().add_axes(ax_cb)
    if bh == True:
        plt.savefig(date + '_Overrepresentation_Annotations_' + dataset + '_bhcorrected.pdf', bbox_inches = 'tight')
    else:
        plt.savefig(date + '_Overrepresentation_Annotations_' + dataset + '.pdf', bbox_inches = 'tight')
    plt.show()

for bh_test in [True, False]:
    barplot_loc_overrep_enzymes(enz_loc_overrep, 'enzymes', bh_test)
    barplot_loc_overrep_enzymes(multilocenz_loc_overrep, 'multilocenzymes', bh_test)
    barplot_loc_overrep_enzymes(scvenz_loc_overrep, 'scvenzymes', bh_test)

