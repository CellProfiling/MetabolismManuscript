### author: Christian Gnann
### map UniProt locations to HPA locations and filter for locations with experimental evidence - output is a csv file with those locations

#%%
import pandas as pd

#%% 
noGOid = pd.read_csv("Uniprot_subcellular_location without_GO-id.csv", engine='python',encoding='utf-8')
#hpago = pd.read_csv("KnownLocationsFromUniprot/GOmapping_HPA.csv", engine='python')
hpago = pd.read_csv("GOmapping_HPA_updated.csv", engine='python', encoding='python')
uniprot = pd.read_csv("20230222_UniprotLocations.csv", engine='python', encoding='python')
date = ''

#%%    
### add ECO code and localizations together to allow for filtering
### first filter rows with no localization or evidence data
uniprot = uniprot[uniprot['location'] != '']
uniprot = uniprot[uniprot['evidence'] != '']
uniprot = uniprot.dropna()

uniprot = uniprot.reset_index(drop=True)

### now add the localization together with the ECO code
uniprot['loc_eco'] = ''
uniprot['loc_eco_comment'] = ''

for i in uniprot.index:
    loclist = uniprot.location[i].split(';')
    ecolist = uniprot.evidence[i].split(';')
    if len(loclist) == len(ecolist):
        loceco = ''
        for x in range(0, len(loclist)):
            loceco = loceco + ';' + loclist[x] + '_' + ecolist[x]
        loceco = loceco.strip(';')
        uniprot.loc_eco[i] = loceco
    else: ### those are examples where there are multiple locations per evidence or the other way round
        uniprot.loc_eco_comment[i] = 'needs special treatment'
        id = uniprot.uniprot_id[i]
        #print(id, ecolist, loclist)

#%%
### filter for localizations with experimental evidence. (ECO code ECO:0000269) --> https://www.uniprot.org/help/evidences
### first get rid of genes with no experimental evidence at all
#uniprot = uniprot[uniprot['loc_eco_comment'] == 'NO EXPERIMENTAL EVIDENCE']
uniprot_filt = uniprot[uniprot['loc_eco'] != '']
print('exclude inconclusive ones: ',uniprot_filt.shape)
uniprot_filt = uniprot_filt[uniprot_filt['evidence'].str.contains('ECO:0000269')]
print('Proteins with any experimental evidence : ', uniprot_filt.shape)

uniprot_filt['loc_exp_evidence'] = ''
for i in uniprot_filt.index:
    allloceco = uniprot_filt.loc_eco[i].split(';')
    exp = ''
    for loceco in allloceco:
        if 'ECO:0000269' in loceco:
            exploc = loceco.split('_')[0]
            #print(exploc)
            exp = exp + ';' + exploc
    exp = exp.strip(';')
    #print(exp)
    uniprot_filt.loc_exp_evidence[i] = exp
uniprot_filt = uniprot_filt[uniprot_filt['loc_exp_evidence'] != '']
print('Double check filtering of experimental evidence: ', uniprot_filt.shape)

#%%        
### deal with GO term hierarchy and extract everything
idlist, namelist, namespacelist, is_a_list, rel_list = [], [], [], [], []
with open("KnownLocationsFromUniprot/GOhierarchy.txt") as openfile:
    wholetext = openfile.read()
    for go in wholetext.split('[Term]'):
        id, name, namespace, is_a, relationship = '', '', '', '', ''
        #print(len(wholetext))
        #print(go)
        for line in go.split('\n'):
            if line.startswith('id:'):
                id = line[4:].strip()
            elif line.startswith('name:'):
                name = line[6:].strip()
            elif line.startswith('namespace:'):
                namespace = line[11:].strip()
            elif line.startswith('is_a:'):
                is_a = is_a + ',' + line[6:16]
            elif line.startswith('relationship: part_of'):
                relationship = relationship + ',' + line[22:32]
        idlist.append(id)
        namelist.append(name)
        namespacelist.append(namespace)
        is_a_list.append(is_a[1:])
        rel_list.append(relationship[1:])

print(len(idlist), len(namelist), len(namespacelist), len(is_a_list))
Dict = {'GOid':idlist, 'location':namelist, 'category':namespacelist, 'family':is_a_list, 'relationship':rel_list}
golevels = pd.DataFrame(Dict)

#%%
### walk up the GO hierarchy for "is_a" and "relationsship" until there is a match with one of the HPA location GO terms
golevels = golevels.set_index('GOid')
hpaterms = set(hpago['GO_term'].tolist())

### parentlevel is an integer, indicating whether we are at direct parent or grandparents; 1 for first generation, 2 for grandparents
def find_parent(original_go, current_go): 
    #print('Function has been called again for the original GO term', original_go)
    ### relationsship (part_of:) is srtonger than family (is_a) so search for this one first
    parent_isa = golevels.family[current_go]
    parent_rel = golevels.relationship[current_go]
    ### make sure to get all parents --> but "part of" comes first in the list
    parent_go_list = (parent_rel + ',' + parent_isa).strip(',').split(',')

    if set(parent_go_list).intersection(hpaterms):
        hpaparent = ",".join(str(s) for s in set(parent_go_list).intersection(hpaterms))
    elif parent_go_list != ['']:
        grandparentlist = []
        for parent in parent_go_list: #### work with relationsship first
            #print('Now checking for new parent', parent)
            preparent = find_parent(go, parent)
            grandparentlist.append(preparent)
        #print('The grandparentset: ',grandparentset)
        grandparentstring = ','.join(str(s) for s in grandparentlist)
        grandparentunique = list(set(grandparentstring.split(',')))
        #print('On to the next one _________________________\n', grandparentlist, '\n', grandparentunique)
        hpaparent = ','.join(str(s) for s in grandparentunique).strip(',')
    else:
        hpaparent = ''
    #print('the hpa parent: ',hpaparent)
    return hpaparent

golevels = golevels[golevels['category'] == 'cellular_component']
golevels['hpa_go_parent'] = ''

for go in golevels.index:
    if go in hpaterms:
        golevels.hpa_go_parent[go] = go
    else:
        golevels.hpa_go_parent[go] = find_parent(go, go)

#%%
### now we'll have to map the GO term to the corresponding HPA location
golevels['hpa_go_parent_final'] = ''
go_dict_hierarchy = pd.Series(hpago.hierarchical_level.values,index=hpago.GO_term).to_dict()
for i in golevels.index:
    allgotermlist = list(filter(None, (golevels.hpa_go_parent[i]).split(',')))
    #print(allgotermlist)
    if len(allgotermlist) == 1:
        """
        pick the current GO term as hpa parent term
        """       
        golevels.hpa_go_parent_final[i] = golevels.hpa_go_parent[i]
    elif len(allgotermlist) >= 2:
        """
        pick the GO term with the lower hierarchy (based on ranking from excel file) as hpa parent term
        """       
        hierarchy_list = []
        for go in allgotermlist:
            level = go_dict_hierarchy.get(go)
            hierarchy_list.append(level)        
        lowest_hierarchy_index = hierarchy_list.index(max(hierarchy_list))
        #hierarchy_set = set(hierarchy_list)
        final_go = allgotermlist[lowest_hierarchy_index]
        golevels.hpa_go_parent_final[i] = final_go
        ###if len(hierarchy_set) == 1:
        ###    print(i, hierarchy_list)
    else:
        """
        those are examples for which there were no HPA parent terms found
        """       
        continue

godict = pd.Series(hpago.location.values,index=hpago.GO_term).to_dict()
golevels['hpa_location'] = golevels['hpa_go_parent_final']
golevels['hpa_location'] = golevels['hpa_location'].map(godict)

#%%
#### first we need to map the UniProt locations to the GO terms used in the HPA
#### download from: ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/external2go/
unicodelist, uniloclist, gocodelist, goloclist = [], [], [], []
#with open("KnownLocationsFromUniprot/uniprot_GO_mapping_subcelllocation.txt") as openfile:
with open("KnownLocationsFromUniprot/20230222_uniprot_GO_mapping_subcelllocation.txt") as openfile:
    lines = openfile.readlines()
    for line in lines:
        if line.startswith('UniProtKB-SubCell:'):
            unicode = line[18:25].strip()
            uniloc = line.split('>')[0][26:].strip()
            gocode = line[-11:].strip('\n')
            goloc = (line.split('>')[1])[1:-14]
            unicodelist.append(unicode)
            uniloclist.append(uniloc)
            gocodelist.append(gocode)
            goloclist.append(goloc)

print(len(unicodelist), len(uniloclist), len(gocodelist), len(goloclist))
Dict = {'UniProt_SL_code':unicodelist, 'UniProt_location_name':uniloclist, 'GOterm':gocodelist, 'GO_location_name':goloclist}
gomap_uniprot = pd.DataFrame(Dict)

#%%
#### map the correct GO term to the uniProt locations:
uni_go_dict = pd.Series(gomap_uniprot.GOterm.values,index=gomap_uniprot.UniProt_location_name).to_dict()
loclistuniprot = set(gomap_uniprot['UniProt_location_name'].unique())
### filter special cases for only human locations --> this is done on the list provided by Fredric
noGOid = noGOid[noGOid['species']=='human']
loclistspecial = set(noGOid['location'].unique())
special_go_dict = pd.Series(noGOid.GO_term.values,index=noGOid.location).to_dict()

### further exceptions are covered here --> manual assignment to the mother HPA GO term by CG
noGOid_manual = pd.read_csv("KnownLocationsFromUniprot/20210611_ManualGOtermAssignment_notinspecialcases.csv", engine='python')
noGOid_manual['location'] = noGOid_manual['location'].str.capitalize()
loclistmanual = set(noGOid_manual['location'].unique())
manualspecial_go_dict = pd.Series(noGOid_manual.parentGOterm.values,index=noGOid_manual.location).to_dict()

uniprot_filt['UniProt_final_locations'] = ''
uniprot_filt['UniProt_Corresponding_GOTerm'] = ''
checklocations = []
for i in uniprot_filt.index:
    allexploc = uniprot_filt.loc_exp_evidence[i].split(';')
    finalloc = ''
    allgo = ''
    for loc_hierarchy in allexploc:
        loc = loc_hierarchy.split(',')[-1]
        #print(loc)
        if loc.capitalize() in loclistuniprot:
            go = uni_go_dict[loc.capitalize()]
            finalloc = finalloc + ';' + loc
            #print(location, '______', go)
        elif loc.capitalize() in loclistspecial:
            #print(location)
            go = special_go_dict[loc.capitalize()]
            finalloc = finalloc + ';' + loc
            #checklocations.append(location)
        elif loc.capitalize() in loclistmanual:
            go = manualspecial_go_dict[loc.capitalize()]
            finalloc = finalloc + ';' + loc
        else:
            continue
            #checklocations.append(location)
        
        allgo = allgo + ';' + go
    finalloc = finalloc.strip(';')
    allgo = allgo.strip(';')
    uniprot_filt.UniProt_final_locations[i] = finalloc
    uniprot_filt.UniProt_Corresponding_GOTerm[i] = allgo



#%%
### for each of the "UniProt GOTerms walk up the GO term hierarchy to find the HPA GO terms
### then generate a new column containing corresponding HPA localization data

### make a dictionary for original GO and parent GO; first filter golevels df to things that are mapped to a parent
allgodict_df = golevels[golevels['hpa_go_parent_final'] != '']
og_parent_dict = pd.Series(allgodict_df.hpa_go_parent_final.values,index=allgodict_df.index).to_dict()
goterms = set(allgodict_df.index.to_list())
### use the previous dictionary to find the parent GO term --> then find the corresponding HPA location via another dictionary
### write to cell
### first remove rows without corresponding GOterm; e.g. location = flagellum axoneme
uniprot_filt = uniprot_filt[uniprot_filt['UniProt_Corresponding_GOTerm'] != '']
uniprot_filt['CorrespondingGOTerm_HPA'] = ''
uniprot_filt['Corresponding_HPA_location'] = ''

checkgo = []
#golevels = golevels.set_index('location')
for i in uniprot_filt.index:
    og_go_list = uniprot_filt.UniProt_Corresponding_GOTerm[i].split(';')
    allparents = []
    hpalocations = []
    for go in og_go_list:
        #print(go)
        if go in goterms:
            parentgo = og_parent_dict[go]
            #print(go, parentgo)
        else: 
            """
            there are three terms that are not contained in the dictionary
            1. GO:0000777 --> old term; replaced by GO:0000776
            2. GO:0031514 --> flagellum --> just delete it
            3. GO:0005737 --> cytoplasm; not sure why this one was skipped
            """
            if go == 'GO:0000777':
                parentgo = 'GO:0000776'
            elif go == 'GO:0031514':
                continue
            elif go == 'GO:0005737':
                parentgo = 'GO:0005829'
        allparents.append(parentgo)
    allparents = list(set(allparents))
    ### now find the corresponding HPA locations
    for go in allparents:
        hpaloc = godict[go]
        hpalocations.append(hpaloc)
    allhpaloc = list(set(hpalocations))
    gohpamapped = ';'.join([str(elem) for elem in allparents])
    lochpamapped = ';'.join([str(elem) for elem in allhpaloc])
    ### add to dataframe
    uniprot_filt.CorrespondingGOTerm_HPA[i] = gohpamapped
    uniprot_filt.Corresponding_HPA_location[i] = lochpamapped

        
        
uniprot_filt.to_csv('KnownLocationsFromUniprot/' + date + '_UniprotKnownLocations_mappedtohpa.csv', index=False)        
        
