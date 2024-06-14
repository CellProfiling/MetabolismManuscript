#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anthony.cesnik
"""
    
import pandas as pd
import numpy as np
import statsmodels.stats.multitest as smm
import scipy.stats as stats

def proportion_test(phosHuman, names_hpamapped, names_subset, names_alt=[]):
    '''
    Test the proportions of kinase families with Fisher exact test
    NB: Kinases lacking a family, i.e. "Other," are dropped here
    '''
    subset = np.isin(phosHuman["SUBSTRATE"], names_subset)
    mapped_minus = np.isin(phosHuman["SUBSTRATE"], names_hpamapped if len(names_alt) == 0 else names_alt) & ~subset # mutually exclusive mapped proteome
    counts_subset = phosHuman["KINASE_FAMILY"][subset].value_counts().drop(labels='Other', errors='ignore')
    counts_mappedminus = phosHuman["KINASE_FAMILY"][mapped_minus].value_counts().drop(labels='Other', errors='ignore')
    labels_subset = counts_subset.index
    labels_mappedminus = counts_mappedminus.index
    proportion_subset = np.array([phosHuman["KINASE_FAMILY"][subset].value_counts()[idx] / sum(phosHuman["KINASE_FAMILY"][subset].value_counts()) for idx in np.arange(len(labels_subset))])
    proportion_mappedminus = np.array([counts_mappedminus.iloc[idx] / sum(phosHuman["KINASE_FAMILY"].iloc[mapped_minus].value_counts()) for idx in np.arange(len(labels_mappedminus))])    

    labels_comb = np.unique(np.concatenate((labels_subset, labels_mappedminus)))
    counts_comb_subset = np.array([counts_subset[np.array(labels_subset) == x][0] if x in labels_subset else 0 for x in labels_comb])
    counts_comb_mappedminus = np.array([counts_mappedminus[np.array(labels_mappedminus) == x][0] if x in labels_mappedminus else 0 for x in labels_comb])
    values_comb_subset = np.array([proportion_subset[np.array(labels_subset) == x][0] if x in labels_subset else 0 for x in labels_comb])
    values_comb_mappedminus = np.array([proportion_mappedminus[np.array(labels_mappedminus) == x][0] if x in labels_mappedminus else 0 for x in labels_comb])
    fisher_comb_subset = np.array([stats.fisher_exact([
        [counts_comb_subset[idx], sum(counts_comb_subset) - counts_comb_subset[idx]], 
        [counts_comb_mappedminus[idx], sum(counts_comb_mappedminus) - counts_comb_mappedminus[idx]]], "greater") for idx in np.arange(len(labels_comb))])
    return labels_comb, counts_comb_subset, values_comb_subset, counts_comb_mappedminus, values_comb_mappedminus, fisher_comb_subset
    

def kinase_families():
    '''Investigate whether there are differences in upstream kinases from mapped proteins'''
    print("Running kinase family analysis")
    kinaseFams = pd.read_csv("input/KinHubKinaseFamilies.csv")
    kinFamDict = dict([(row[7], row[4]) for idx, row in kinaseFams.iterrows()])
    phosphositeplus = pd.read_csv("input/Kinase_Substrate_Dataset_240517_withoutheader.txt.gz", sep="\t")
    phosphositeplus["KINASE_FAMILY"] = [kinFamDict[x] if x in kinFamDict else "Other" for x in phosphositeplus["KIN_ACC_ID"]]
    phosHuman = phosphositeplus[phosphositeplus["SUB_ORGANISM"] == "human"]

    # Get HPA information
    hpamapped = pd.read_csv("input/HPAexport_ProteinProperties_ForAnthony.csv")
    hpa_all = pd.read_csv("input/proteinatlas_v23.tsv.gz", sep="\t")
    names_hpamapped = hpa_all["Gene"][hpa_all['Subcellular location'].notna() & hpa_all['Subcellular location'].str.strip().astype(bool)]
    print(f"{len(names_hpamapped)}: including this many proteins in background")
    names_ccdprotein = hpamapped[hpamapped["IsCCD"]]["Gene_name"]
    names_enzymes = hpamapped[hpamapped["enzyme"]]["Gene_name"]
    names_scv = hpamapped[hpamapped["SCV"] & hpamapped["enzyme"]]["Gene_name"]
    names_nonccd = hpamapped[hpamapped["IsNonCCD"]]["Gene_name"]

    # Are there any overrepresented kinase families upstream phosphosites on CCD or non-CCD proteins?
    labels_comb_ccd, counts_comb_ccd, values_comb_ccd, counts_comb_mappedminus_ccd, values_comb_mappedminus_ccd, fisher_comb_ccd = proportion_test(
        phosHuman, names_hpamapped, names_ccdprotein)
    labels_comb_enzymes, counts_comb_enzymes, values_comb_enzymes, counts_comb_mappedminus_enzymes, values_comb_mappedminus_enzymes, fisher_comb_enzymes = proportion_test(
        phosHuman, names_hpamapped, names_enzymes)
    labels_comb_scv, counts_comb_scv, values_comb_scv, counts_comb_mappedminus_scv, values_comb_mappedminus_scv, fisher_comb_scv = proportion_test(
        phosHuman, names_hpamapped, names_scv)
    labels_comb_nonccd, counts_comb_nonccd, values_comb_nonccd, counts_comb_mappedminus_nonccd, values_comb_mappedminus_nonccd, fisher_comb_nonccd = proportion_test(
        phosHuman, names_hpamapped, names_nonccd)
    labels_comb_scvenzVsEnz, counts_comb_scvenzVsEnz, values_comb_scvenzVsEnz, counts_comb_mappedminus_scvenzVsEnz, values_comb_mappedminus_scvenzVsEnz, fisher_comb_scvenzVsEnz = proportion_test(
        phosHuman, names_enzymes, names_scv)
    
    allfisher = np.array(np.concatenate((fisher_comb_ccd[:,1], fisher_comb_enzymes[:,1], 
                                        # fisher_comb_scv[:,1], 
                                        # fisher_comb_nonccd[:,1],
                                        fisher_comb_scvenzVsEnz[:,1] 
                                        )), dtype=float)
    reject_bonf, pvals_corrected_bonf, _,_ = smm.multipletests(allfisher, alpha=0.05, method="bonferroni")
    pd.DataFrame({
        "Group" : np.concatenate((["ccd"] * len(labels_comb_ccd), 
                                    ["enzymes"] * len(labels_comb_enzymes),
                                    # ["scv-enzymes"] * len(labels_comb_scv),
                                    # ["nonccd"] * len(labels_comb_nonccd),
                                    ["scv-enzymes-vs-nonscv-enzymes"] * len(labels_comb_scvenzVsEnz))),
        "KinaseFamily" : np.concatenate((labels_comb_ccd, labels_comb_enzymes, 
                                        #  labels_comb_scv,
                                        # labels_comb_nonccd, 
                                        labels_comb_scvenzVsEnz)),
        "PhosphositeCount_MappedProteome" : np.concatenate((counts_comb_mappedminus_ccd, counts_comb_mappedminus_enzymes, 
                                                # counts_comb_mappedminus_scv, 
                                                # counts_comb_mappedminus_nonccd, 
                                                counts_comb_mappedminus_scvenzVsEnz)),
        "FractionPhosphositesDownstreamOfKinase_MappedProteome" : np.concatenate((values_comb_mappedminus_ccd, values_comb_mappedminus_enzymes, 
                                            # values_comb_mappedminus_scv, 
                                            # values_comb_mappedminus_nonccd,
                                            values_comb_mappedminus_scvenzVsEnz)),
        "PhosphositeCount" : np.concatenate((counts_comb_ccd, counts_comb_enzymes, 
                                                    # counts_comb_scv, 
                                                    # counts_comb_nonccd, 
                                                    counts_comb_scvenzVsEnz)),
        "FractionPhosphositesDownstreamOfKinase" : np.concatenate((values_comb_ccd, values_comb_enzymes, 
                                            # values_comb_scv, 
                                            # values_comb_nonccd, 
                                            values_comb_scvenzVsEnz)),
        "FisherPValue" : np.concatenate((fisher_comb_ccd[:,1], fisher_comb_enzymes[:,1], 
                                        # fisher_comb_scv[:,1], 
                                        # fisher_comb_nonccd[:,1], 
                                        fisher_comb_scvenzVsEnz[:,1])),
        "FisherPValue_BonfCorrected" : pvals_corrected_bonf,
        "FisherPValue_BonfAlpha0.05Pass" : reject_bonf,
        }).to_csv("output/upstreamKinaseResults.csv", index = False)
        
def main():
   kinase_families() 

if __name__ == "__main__":
    main()