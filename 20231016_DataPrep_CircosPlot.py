#%%
import pandas as pd
import itertools


#%%
date = '20231016'

hpa = pd.read_csv("", engine='python') ### hpa v23
pathway = pd.read_csv("", engine='python') ## file with all reactome ids and kegg ids in human genes --> filter for metabolic terms 


#%%
###add column "enzyme" to hpa dataframe and use it to filter for the enzyme data
hpa['enzyme'] = ''
for i in hpa.index:
    if hpa.Gene_name[i] in pathway['Gene_name'].tolist():
       hpa.enzyme[i] = 'TRUE'
    else:
        hpa.enzyme[i] = 'FALSE'


#%%
### merge substructures into organelles (fibrillar center and nucleoi and nuci rim, ...)
locations = ['Actin filaments', 'Aggresome', 'Cell Junctions', 'Centriolar satellite', 'Centrosome',
                'Cleavage furrow', 'Cytokinetic bridge', 'Cytoplasmic bodies', 'Cytosol', 
                'Endoplasmic reticulum', 'Endosomes', 'Focal adhesion sites', 'Golgi apparatus', 
                'Intermediate filaments', 'Kinetochore', 'Lipid droplets', 'Lysosomes', 'Microtubule ends',
                'Microtubules', 'Midbody', 'Midbody ring', 'Mitochondria', 'Mitotic chromosome', 'Mitotic spindle',
                'Nuclear bodies', 'Nuclear membrane', 'Nuclear speckles', 'Nucleoli', 'Nucleoli fibrillar center',
                'Nucleoli rim', 'Nucleoplasm','Peroxisomes', 'Plasma membrane', 'Rods & Rings', 'Vesicles']
groupedlocations = ['Actin filaments','Centrosome','Cytoplasmic Aggregates', 'Cytosol','Endoplasmic reticulum',
                    'Golgi apparatus','Intermediate filaments','Microtubules','Mitochondria','Nuclear bodies',
                    'Nuclear membrane','Nuclear speckles','Nucleoli','Nucleoplasm','Plasma membrane', 'Vesicles']


#%%
def circosplotinput_allloc(df,datasetstring):
    input_df = pd.DataFrame(columns = ['Loc1', 'Loc2','bothloc', 'flow'])
    for i in groupedlocations:
        for g in groupedlocations:
            locpair = i + ',' + g
            input_df = input_df.append({'Loc1':i,'Loc2':g,'bothloc':locpair,
                                        'flow':0}, ignore_index=True)
    for i in input_df.index:
        locpairlist = sorted(input_df.bothloc[i].split(','))
        input_df.bothloc[i] = ';'.join(locpairlist)
    input_df.drop_duplicates(subset ="bothloc", keep = 'first', inplace = True) 
    for i in df.index:
        loc_list = sorted(df.All_GroupedLocations[i].split(';'))
        locations = ';'.join(loc_list)   ###alphabetically sorted to match other stuff
        #print(loc_list, locations)
        if len(loc_list) == 1:
            locations = locations + ';' + locations
            idx = input_df.index[input_df['bothloc'] == locations].tolist()[0]
            input_df.flow[idx] = input_df.flow[idx] + 1
            ### +1 where loc list[0] = locatiosn and loc list[1] = locations
        elif len(loc_list) == 2:
            idx = input_df.index[input_df['bothloc'] == locations].tolist()[0]
            input_df.flow[idx] = input_df.flow[idx] + 1
        else:
            #### for all double combinations in location --> run the same as for the elif statement
            #print('New gene with more than 2 locations: ',loc_list)
            for locpair in itertools.combinations(loc_list, 2):
                locations = ';'.join(locpair)
                idx = input_df.index[input_df['bothloc'] == locations].tolist()[0]
                input_df.flow[idx] = input_df.flow[idx] + 1
    ### add color to dataframe
    ### first define subcompartments
    nuc = ['Nuclear bodies', 'Nuclear membrane', 'Nuclear speckles', 'Nucleoli', 'Nucleoplasm']
    cyto = ['Actin filaments','Centrosome','Cytoplasmic Aggregates', 'Cytosol', 'Intermediate filaments', 'Microtubules', 'Mitochondria']
    sec = ['Endoplasmic reticulum', 'Golgi apparatus', 'Plasma membrane', 'Vesicles']
    input_df['color'] = ''
    for i in input_df.index:
        loc1 = input_df.Loc1[i]
        loc2 = input_df.Loc2[i]
        if loc1 == loc2:
            col = 'white'
            '''
            ### this parts reduces the amount of nodes by two to represent the real ratio in the circos plot later
            ### not needed anymore though
            no_nodes = input_df.flow[i]
            input_df.flow[i] = math.ceil(no_nodes/2)
            '''
        elif loc1 in nuc and loc2 in nuc:
            col = 'blue'
        elif loc1 in nuc and loc2 in cyto:
            col = 'magenta'
        elif loc1 in nuc and loc2 in sec:
            col = 'cyan'
        elif loc1 in cyto and loc2 in nuc:
            col = 'magenta'
        elif loc1 in cyto and loc2 in cyto:
            col = 'red'
        elif loc1 in cyto and loc2 in sec:
            col = 'yellow'
        elif loc1 in sec and loc2 in nuc:
            col = 'cyan'
        elif loc1 in sec and loc2 in cyto:
            col = 'yellow'
        elif loc1 in sec and loc2 in sec:
            col = 'green'
        input_df.color[i] = col
       
    filename = date + '_CirclizeDataInput_FullDataset_'+ datasetstring + '.csv'
    input_df.to_csv(filename,index=False)
    
# %%
### generate subgroups for which we should prepare circos plots
enz_hpa = hpa[hpa['enzyme']=='TRUE']

#%%
### actual file generation
def run_circosplotprep(df, dataset):
    circosplotinput_allloc(df, dataset)

run_circosplotprep(enz_hpa, 'enzymes')

# %%
### generate pathway specific files
for path in pathway['PathwayGroup'].unique().tolist():
    subset = pathway[pathway['PathwayGroup'] == path]
    genelist = subset['ENSG'].unique()
    path = path.replace('/Biological Oxidations', '')   
    path = path.replace(' ', '_')   
    hpa_path = hpa[hpa['Gene'].isin(genelist)]
    print(path, hpa_path.shape) 
    run_circosplotprep(hpa_path, path)


