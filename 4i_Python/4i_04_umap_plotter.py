from sklearn.preprocessing import StandardScaler # v 1.4.1.post1
import matplotlib.pyplot as plt # v 3.8.2
import pandas as pd # v 2.1.4
import umap # v 0.1.1
import seaborn as sns # v 0.13.2
import os
from scipy.stats import zscore # v 1.11.4
import re
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score



pd.set_option('display.max_columns', None)

path = "X:\CELL_PROFILING\\research\Alina\\4i\\4i_U2OS_set01\\"
files = os.listdir(path)
reducer = umap.UMAP(random_state=42)
data_path = [f for f in files if ("cell by cell" in f)&(".csv" in f)]


df_raw = pd.concat(pd.read_csv(path+csv) for csv in data_path)
df_raw["Well"] = df_raw["Cell ID"].str.extract(r'(r\d\dc\d\d)')

df = df_raw.loc[(df_raw["Is complete"] == True) & (df_raw["Nucleus size"] > 1500)]
wells = set(df["Well"].values)

df=df.dropna()

tubulin_correction_features = [tcf for tcf in df.columns if (("intensity" in tcf) or ("avg" in tcf)) and ("Average intensity of tubulin" not in tcf)]
for tcf in tubulin_correction_features:
    r = re.search(r'\d\d_',tcf).group()
    print(r)
    print(tcf)
    df[tcf] = df[tcf]/df[r+"Average intensity of tubulin"]


for well in wells:
    print(well)
    well_mask = df["Well"].isin([well])
    df.loc[well_mask, df.columns[~df.columns.isin(["Cell ID", "Is complete", "Well"])]] = df.loc[well_mask, df.columns[~df.columns.isin(["Cell ID", "Is complete", "Well"])]].apply(zscore)


features = [f for f in df.columns if "marker" in f]


df_data = df[["Cell size", "Nucleus size", "NCR"]+features].values


scaled_df_data = StandardScaler().fit_transform(df_data)
embedding = reducer.fit_transform(scaled_df_data)


plt.scatter(embedding[:,0],
                embedding[:,1],
                alpha=0.6,
                s=15)
plt.gca().set_aspect('equal','datalim')
plt.show()

dapi_tub = [f for f in df.columns if ("DAPI" in f) or ("tubulin" in f)]

parameters = df.drop(["Cell ID", "Is complete", "Well"], axis=1).columns
hpa_parameters = [p for p in parameters if "HPA" in p]
shared_parameters = [p for p in parameters if "HPA" not in p]
print(shared_parameters)

df[shared_parameters]=df[shared_parameters].apply(zscore)
for w in wells:
    df.loc[df["Well"]==w, hpa_parameters] = df.loc[df["Well"]==w, hpa_parameters].apply(zscore)


i = 0

for w in wells:
    for parameter in parameters:
        measurements = df[parameter]
        gradient = sns.color_palette("coolwarm", as_cmap=True)
        if "marker" in parameter:
            plt.scatter(embedding[:,0],
                        embedding[:,1],
                        c = "silver",
                        cmap=gradient,
                       vmin=-1,
                       vmax=1,
                        alpha=[1 if (w not in x) else 0 for x in df["Well"]],
                        s=15)
        if "marker" in parameter:
            plt.scatter(embedding[:,0],
                        embedding[:,1],
                        c = [measurements],
                        cmap=gradient,
                       vmin=-1,
                       vmax=1,
                        alpha=[0.6 if (w in x) else 0 for x in df["Well"]],
                        s=15)
            plt.colorbar()
            plt.gca().set_aspect('equal','datalim')
            plt.title(parameter+"_"+w)
            plt.show()
    i=1

dfp = pd.DataFrame(dict(x=embedding[:,0], y=embedding[:,1], z=df.index))

silhouettes = []
for i in range(3,25):
    kmeans = KMeans(n_clusters=i, random_state=42)
    cluster_labels = kmeans.fit_predict(embedding).astype(str)
    silhouette = silhouette_score(embedding, cluster_labels, random_state=42)
    print(str(i) + ": " + str(silhouette))
    silhouettes.append(silhouette)

#kmeans = KMeans(n_clusters=silhouettes.index(max(silhouettes))+3, random_state=42)
kmeans = KMeans(n_clusters=18, random_state=42)
cluster_labels = kmeans.fit_predict(embedding).astype(str)

df["Cluster"]=cluster_labels
df.to_excel(path+"clustering results.xlsx")
