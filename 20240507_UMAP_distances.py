#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator
import statsmodels.api as sm


#%%
### user input
date = ''
hpa = pd.read_csv('20231012_HPAv23_scvmultiloccomp.csv', engine='python') ###HPA v23
enzymes = pd.read_csv('20230111_Human1_ProteinMapping.csv', delimiter=',') ### list of enzymes and pathway annotations from Human1
dist = pd.read_csv('20230111_Human1_ProteinMapping_Pathway_avgcosinedistance.csv', delimiter=',') ### cosine distances per pathway calculated from HPAv23 UMAP
dist_group = pd.read_csv('20230111_Human1_ProteinMapping_PathwayGroup_avgcosinedistance.csv', delimiter=',') ### cosine distances per pathway group calculated from HPAv23 UMAP

# %%
enzlist = enzymes['ENSG'].unique()
hpa['enzyme'] = False
for i in hpa.index:
    if hpa.at[i, 'Gene'] in enzlist:
        hpa.at[i, 'enzyme'] = True
    
#%%
### add the number of proteins per pathway to the dist dataframe
dist['number_enzymes'] = 0
for i in dist.index:
    currentpath = dist.Pathway[i]
    enzlist = enzymes[enzymes['Pathway'] == currentpath]['ENSG']
    dist.number_enzymes[i] = len(hpa[hpa['Gene'].isin(enzlist)])

# %%
### boxplot of cosine distance of pathways and pathway groups vs random
distlist = dist['distance_cosine']
distgrouplist = dist_group['distance_cosine']
randomdistlist = dist['distance_random_mean']

def boxplot_dist(dist_path, dist_group, cutoff=1):
    dist_path = dist_path[dist_path['number_enzymes']>=cutoff].reset_index()
    df_path = dist_path[['Pathway', 'distance_cosine']]
    df_random = dist_path[['Pathway', 'distance_random_mean']]
    df_pg = dist_group[['PathwayGroup', 'distance_cosine']]

    df_path['label'] = 'Pathway'
    df_random['label'] = 'Random'
    df_pg['label'] = 'PathwayGroup'

    df_random = df_random.rename(columns={'distance_random_mean':'distance_cosine'})
    
    df = pd.concat([df_path, df_pg, df_random])
    
    labelorder = ['Pathway', 'PathwayGroup', 'Random']
    fig = plt.figure(figsize =(7, 7)) 
    # Creating axes instance 
    ax = fig.add_axes([0, 0, 1, 1]) 
    bp = sns.boxplot(y='distance_cosine', x='label', data=df, order = labelorder, palette="colorblind", showfliers = False)
    bp.set_ylim(0,1)
    stat_boxpairs = [('Pathway', 'PathwayGroup'), ('Pathway', 'Random'), ('PathwayGroup', 'Random')]
    annotator = Annotator(ax=bp, pairs=stat_boxpairs,y='distance_cosine', x='label', data=df, order=labelorder)
    annotator.configure(test='Mann-Whitney', verbose=True, text_format='full').apply_and_annotate()
    
    filename = f"{date}_CosineDistance_cutoff{str(cutoff)}.pdf"
    plt.savefig(filename, bbox_inches="tight")
    plt.show()

for c in [1,2,4,8]:
    boxplot_dist(dist, dist_group, cutoff=c)

# %%
### checking for outliers in a regression model
### first process the dataframe properly
dist['number_enzymes_log'] = np.log2(dist['number_enzymes'])
dist.replace([np.inf, -np.inf], np.nan, inplace=True)
# drop rows with NaN values
dist = dist.dropna().reset_index()

#%%
# Fit the linear regression model
#model = sm.OLS(dist['distance_euclidean'], dist['number_enzymes_log']).fit()
# Fit the linear regression model
def linear_regression_model_scatter(full_df, whichcol, cutoff, dataset = ''):
    '''
    cutoff --> number of enzymes per pathway
    linear regression model of cosine distance over log2 number of enzymes
    scatterplot colored by significant outliers and a second plot colored by selected pathways
    '''
    print(f"cutoff used for number of enzymes per pathway: {cutoff}")
    df = full_df[full_df['number_enzymes']>=cutoff].reset_index()
    X = sm.add_constant(df["number_enzymes_log"])
    y = df[whichcol]
    model = sm.OLS(y, X).fit()

    # Calculate the predicted values and confidence intervals
    predict = model.get_prediction(X)
    df['pred_mean'] = predict.predicted_mean
    ci = predict.conf_int(alpha=0.05)

    # Calculate the standardized residuals
    residuals = model.resid
    std_resid = abs(residuals / np.std(residuals))

    # Identify potential outliers using a threshold (e.g., 2 or 3 standard deviations)
    threshold = 2
    outlier_indices = np.where(std_resid > threshold)[0]
    df['outlier'] = 'N'
    df.loc[outlier_indices, 'outlier'] = 'Y'
    plt.figure(figsize=(6, 6))
    # Plot the original data with the linear regression line and confidence intervals
    s = sns.scatterplot(x=df["number_enzymes_log"], y=df[whichcol], hue = df["outlier"], palette = {'Y':'orange', 'N':'gray'})
    plt.plot(df["number_enzymes_log"], df["pred_mean"], color="black", label="Linear regression line")
    plt.fill_between(df["number_enzymes_log"], ci[:, 0], ci[:, 1], alpha=0.2, color="lightgray", label="95% Confidence Interval")
    plt.xlabel("Number of enzymes (log2)")
    plt.ylabel("Distance (cosine)")
    plt.ylim(0,1)
    plt.title("Linear regression with confidence intervals")
    plt.legend()
    filename = f"{date}_ScatterDistanceOverNumberEnzymes{dataset}_regression_cutoff{str(cutoff)}.pdf"
    plt.savefig(filename, bbox_inches="tight")
    print(df.head())
    for line in range(0,df.shape[0]):
        if not df['outlier'][line] == 'N':
            df['Pathway'][line]
            s.text(df['number_enzymes_log'][line]+0.25, df[whichcol][line], 
                df['Pathway'][line], horizontalalignment='left', 
                size='small', color='black')
    filename = f"{date}_ScatterDistanceOverNumberEnzymes{dataset}_regression_labelledoutliers_cutoff{str(cutoff)}.pdf"
    plt.savefig(filename, bbox_inches="tight")
    print('pval: ',model.pvalues[1])
    print('R: ', model.rsquared)
    plt.close()
    #### now plot but highlight manually
    # Plot the original data with the linear regression line and confidence intervals
    df['highlight']= 'N'
    for i in df.index:
        if df.at[i, 'Pathway'] in ['Oxidative phosphorylation', 'Cholesterol metabolism', 'Aminoacyl-tRNA biosynthesis']:
            df.at[i,'highlight']='Y'
    plt.figure(figsize=(6, 6))
    s = sns.scatterplot(x=df["number_enzymes_log"], y=df[whichcol], hue = df["highlight"], palette = {'Y':'orange', 'N':'gray'}, s=50)
    plt.plot(df["number_enzymes_log"], df["pred_mean"], color="black", label="Linear regression line")
    plt.fill_between(df["number_enzymes_log"], ci[:, 0], ci[:, 1], alpha=0.2, color="lightgray", label="95% Confidence Interval")
    plt.xlabel("Number of enzymes (log2)")
    plt.ylabel("Distance (Cosine)")
    plt.legend()
    plt.ylim(0,1)
    for line in range(0,df.shape[0]):
        if not df['highlight'][line] == 'N':
            s.text(df['number_enzymes_log'][line]+0.25, dist['distance_cosine'][line], 
                df['Pathway'][line], horizontalalignment='left', 
                size='small', color='black')
    filename = f"{date}_ScatterDistanceOverNumberEnzymes{dataset}_regression_selectedpathways_cutoff{str(cutoff)}.pdf"
    plt.savefig(filename, bbox_inches="tight")

for c in [2]: #[1,2,4,8]
    linear_regression_model_scatter(dist, 'distance_cosine', c)



# %%
