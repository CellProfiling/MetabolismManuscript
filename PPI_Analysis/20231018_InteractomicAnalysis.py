### author Christian Gnann
### generate PPI interaction network from a custom database (OpenCell, BioPlex, Hein et al)

#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from itertools import combinations

### for network plots
import networkx as nx
import requests
from matplotlib.lines import Line2D
from mycolorpy import colorlist as mcp

# %%
date = ''
hpa = pd.read_csv('20231012_HPAv23_scvmultiloccomp.csv', engine='python')
pathway = pd.read_csv("../EnzymeLists/20230111_Human1_ProteinMapping.csv", engine='python') ## file with all reactome ids and kegg ids in human genes --> filter for metabolic terms 
enzymelist = pathway['Gene_name'].unique()
proteome_merged_db = pd.read_csv('20231012_ProteomeInteractors.csv', engine='python') ### contains all interactors from the databases (Hein 2015, OpenCell and BioPlex)

#%%
hpa['enzyme'] = False
for i in hpa.index:
    if hpa.Gene_name[i] in set(enzymelist):
        hpa.enzyme[i] = True

#%%
### first define locations list so we're able to assign compartments
hpa_loc = ['Actin filaments', 'Aggresome', 'Cell Junctions', 'Centriolar satellite', 'Centrosome',
                'Cleavage furrow', 'Cytokinetic bridge', 'Cytoplasmic bodies', 'Cytosol', 
                'Endoplasmic reticulum', 'Endosomes', 'Focal adhesion sites', 'Golgi apparatus', 
                'Intermediate filaments', 'Kinetochore', 'Lipid droplet', 'Lysosomes', 'Microtubule ends',
                'Microtubules', 'Midbody', 'Midbody ring', 'Mitochondria', 'Mitotic chromosome', 'Mitotic spindle',
                'Nuclear bodies', 'Nuclear membrane', 'Nuclear speckles', 'Nucleoli', 'Nucleoli fibrillar center',
                'Nucleoli rim', 'Nucleoplasm','Peroxisomes', 'Plasma membrane', 'Rods & Rings', 'Vesicles']
hpa_cyto_list = ['Actin filaments', 'Aggresome', 'Cell Junctions', 'Cytoplasmic bodies', 'Cytosol', 
                'Endoplasmic reticulum', 'Endosomes', 'Focal adhesion sites', 'Golgi apparatus', 
                'Intermediate filaments', 'Lipid droplet', 'Lysosomes', 'Microtubule ends',
                'Microtubules', 'Mitochondria', 'Peroxisomes', 'Plasma membrane', 'Rods & Rings', 'Vesicles']
hpa_nuc_list = ['Kinetochore', 'Mitotic chromosome', 'Nuclear bodies', 'Nuclear membrane', 'Nuclear speckles', 'Nucleoli', 'Nucleoli fibrillar center',
                'Nucleoli rim', 'Nucleoplasm']

### Step 3: for each interactor, check location/compartment
'''
make dictionary of compartment for all open cell and HPA data that can be called in the following function
'''
hpa_compdict = dict(zip(hpa.Gene_name, hpa.compartment))


#%%
### generate PPI network plots
def pull_string_data(protein_list):
    '''
    pulls interactions for a list of proteins from string
    the interaction dataframe can be filtered for strength of interaction, ...
    returns a dataframe with the two interactors as well as the interaction score 
    '''
    proteins = '%0d'.join(protein_list)
    url = 'https://string-db.org/api/tsv/network?identifiers=' + proteins + '&species=9606'
    r = requests.get(url)

    lines = r.text.split('\n') # pull the text from the response object and split based on new lines
    data = [l.split('\t') for l in lines] # split each line into its components based on tabs
    # convert to dataframe using the first row as the column names; drop empty, final row
    df = pd.DataFrame(data[1:-1], columns = data[0])
    interactions = df[['preferredName_A', 'preferredName_B', 'score']]
    #for i in protein_list:

    return interactions
        

def rescale(l,newmin,newmax):
    '''
    this function is needed if you wnt to adjust the size of the nodes, thickness of the edges and so on
    '''
    arr = list(l)
    if not max(arr) == min(arr):
        return [(x-min(arr))/(max(arr)-min(arr))*(newmax-newmin)+newmin for x in arr]
    else:
        return [(x-min(arr))/(max(arr) + 0.01 - min(arr))*(newmax-newmin)+newmin for x in arr]

def network_plot(int_df, gene, dataset, isbaitorprey):
    G=nx.Graph(name='Protein Interaction Graph')
    interactions = np.array(int_df)
    for i in range(len(interactions)):
        interaction = interactions[i]
        a = interaction[0] # protein a node
        b = interaction[1] # protein b node
        w = float(interaction[2]) # score as weighted edge where high scores = low weight
        G.add_weighted_edges_from([(a,b,w)]) # add weighted edge to graph
    
    print(G.degree(v) for v in G)
    # use the matplotlib viridis colormap
    #graph_colormap = cm.get_cmap('viridis', 12)
    # node color varies with Degree
    #c = rescale([G.degree(v) for v in G],0.01,0.9)
    #c = [graph_colormap(i) for i in c]
    compcolordict = {'nuc':'#084594', 'cyto':'#FFA500', 'cell':'#008E89'}
    c=[]
    for node in G:
        if node == gene:
            c.append('#900603')
        elif hpa_compdict.get(node) in ['nuc', 'cyto', 'cell']:
            c.append(compcolordict.get(hpa_compdict.get(node)))
        else:
            c.append('grey')
    
    #c = [compcolordict.get(hpa_compdict.get(node)) if hpa_compdict.get(node) in ['nuc', 'cyto', 'cell']  else 'grey' for node in G]

    # node size varies with betweeness centrality - map to range [10,100] 
    bc = nx.betweenness_centrality(G) # betweeness centrality
    s =  rescale([v for v in bc.values()],1500,2500)
    # edge width shows 1-weight to convert cost back to strength of interaction 
    ew = rescale([float(G[u][v]['weight']) for u,v in G.edges],0.1,4)
    # edge color also shows weight
    #ec = rescale([float(G[u][v]['weight']) for u,v in G.edges],0.1,1)
    #ec = [graph_colormap(i) for i in w]

    pos = nx.spring_layout(G)
    fig = plt.figure(figsize=(10,10),facecolor=[0.0,0.0,0.0,0.0])
    nx.draw_networkx(G, pos=pos, with_labels=True, node_color=c, node_size=s,width=ew,
                    font_color='white',font_weight='bold',font_size='9') ### edge_color= ec --> to color the nodes as well
    plt.axis('off')
    

    #cbar = plt.colorbar(heatmap)
    #cbar.ax.set_ylabel('edge weight',labelpad=15,rotation=270)


    legendhandles = []
    for key in compcolordict:
        patch = Line2D([0], [0], marker='o', color=compcolordict.get(key), label=key, 
                        markerfacecolor=compcolordict.get(key), markersize=12)
        legendhandles.append(patch)
    legendhandles.append(Line2D([0], [0], marker='o', color='#900603', label='ProteinOfInterest', 
                        markerfacecolor='#900603', markersize=12))
    legendhandles.append(Line2D([0], [0], marker='o', color='grey', label='NoLocationDataHPA', 
                        markerfacecolor='grey', markersize=12))
    plt.title('Interactors ' + gene + '_' + dataset, fontsize = 30) #, fontdict=20
    plt.legend(handles=legendhandles)
    ### save as pdf and png
    plt.savefig('figures/' + date + '_InteractionNetwork/' + dataset + '/' + date + '_' + dataset + '_' + isbaitorprey + '_' + gene + '_interactors.pdf' )
    plt.savefig('figures/' + date + '_InteractionNetwork/' + dataset + '/' + date + '_' + dataset + '_' + isbaitorprey + '_' + gene + '_interactors.png' )
    #plt.show()
    plt.close(fig)    
    
def pairwise_interactor_check(protein_list, db_df, string_df): #, db_df
    '''
    takes an input of a protein list and checks whether there is data for those proteins in the interactomics databases(s)
    this has to happen with the original data though cause not all interactors/baits are covered in the processed dataframes (filtered for enzymes)

    '''
    allgenes = set(db_df['Gene_name'].unique())
    prot_pairs = list(combinations(sorted(protein_list), 2))
    #print(prot_pairs)
    nameA_list = []
    nameB_list = []
    scorelist = []
    for pair in prot_pairs:
        #print('\nnext pair: ',pair)
        #print(pair[0], pair[1])
        if pair[0] in allgenes and pair[1] in allgenes:
            intlist = set(db_df[db_df['Gene_name']==pair[0]].reset_index().at[0,'allinteractors'].split(';'))
            #print(intlist)
            if pair[1] in intlist:
                '''
                first check whether this pair is already in the string dataframe
                
                '''
                #print('\nnext pair: ',pair)
                #print('do something')
                pairstring = pair[0] + ';' + pair[1]
                if pairstring in string_df['interactor_pair'].unique():
                    #print('will have to chnage the score')
                    continue
                elif pairstring in string_df['interactor_pair_reverse'].unique():
                    #print('will have to chnage the score in the interaction dataframe')
                    continue
                else:
                    #print('New interaction found')
                    nameA_list.extend([pair[0], pair[0]])
                    nameB_list.extend([pair[1], pair[1]])
                    scorelist.extend([1.0, 1.0])
            else:
                #print('No additional interaction has been found')
                continue
        else:
            #print('we are exiting')
            continue
    #print(type(string_df))
    string_df = string_df[['preferredName_A', 'preferredName_B', 'score']]
    if not len(nameA_list) == 0:
        novel_int = pd.DataFrame(list(zip(nameA_list, nameB_list, scorelist)), columns =['preferredName_A', 'preferredName_B', 'score'])
        string_df = pd.concat([string_df, novel_int])
    return string_df

def generate_network_plot(protein_list, gene, dataset, baitorprey):
    interaction_df = pull_string_data(protein_list)
    '''
    now add the interaction between gene and protein list manually to make sure they are included in the network
    1. for this we first remove all interactions with the gene from the list
    2. then add the interactions back for the entire gene list (2x) with a score of 0.9 --> high confidence
    '''
    int_df = interaction_df[(interaction_df['preferredName_A'] != gene) & (interaction_df['preferredName_B'] != gene)]
    nameA_list = []
    nameB_list = []
    scorelist = []
    for prey in protein_list:
        if not prey == gene:
            nameA_list.extend([gene, gene])
            nameB_list.extend([prey, prey])
            scorelist.extend([1.0, 1.0])
    updated_prey_int = pd.DataFrame(list(zip(nameA_list, nameB_list, scorelist)), columns =['preferredName_A', 'preferredName_B', 'score'])
    int_df = pd.concat([int_df, updated_prey_int])
    int_df['interactor_pair'] = int_df["preferredName_A"] + int_df["preferredName_B"]
    int_df['interactor_pair_reverse'] = int_df["preferredName_B"] + ';' + int_df["preferredName_A"]
    #return int_df
    if not int_df.empty:
        print('Now processing: ',gene)
        updated_int_df = pairwise_interactor_check(protein_list, proteome_merged_db, int_df)
        network_plot(updated_int_df,  gene, dataset, baitorprey) #dataset

def run_generate_network_plt_df(df, dataset):
    for i in df.index:
        gl = (df.at[i,'allinteractors']).split(';')
        gene = df.at[i, 'Gene_name']
        if df.at[i, 'IsBait'] == 'Y':
            baitprey = 'bait'
        else:
            baitprey = 'prey'
        if len(gl) >= 2:
            generate_network_plot(gl, gene, dataset, baitprey)
            #output = generate_network_plot(gl, gene, dataset)
            #return output

run_generate_network_plt_df(proteome_merged_db, 'MergedData_Exp_evidence')