from sklearn.preprocessing import StandardScaler # v 1.4.1.post1
import matplotlib.pyplot as plt # v 3.8.2
import pandas as pd # v 2.1.4
import umap # v 0.1.1
import os
from scipy.stats import zscore # v 1.11.4
import numpy as np # v 1.26.3
import re



pd.set_option('display.max_columns', None)

path = "X:\CELL_PROFILING\\research\Alina\\4i\\4i_U2OS_set01\\"
files = os.listdir(path)


reducer = umap.UMAP(random_state=42)

data_path = [f for f in files if ("cell by cell" in f)&(".csv" in f)]


df_raw = pd.concat(pd.read_csv(path+csv) for csv in data_path)
df_raw["Well"] = df_raw["Cell ID"].str.extract(r'(r\d\dc\d\d)')

df = df_raw.loc[(df_raw["Is complete"] == True) & (df_raw["Nucleus size"] > 1500)]
df = df.drop([c for c in df.columns if "MDA" in c], axis=1)
df = df.drop([c for c in df.columns if "Morph" in c], axis=1)
df=df.dropna()

wells = set(df["Well"].values)

tubulin_correction_features = [tcf for tcf in df.columns if (("intensity" in tcf) or ("avg" in tcf)) and ("Average intensity of tubulin" not in tcf)]
for tcf in tubulin_correction_features:
    r = re.search(r'\d\d_',tcf).group()
    df[tcf] = df[tcf]/df[r+"Average intensity of tubulin"]

for well in wells:
    well_mask = df["Well"].isin([well])
    df.loc[well_mask, df.columns[~df.columns.isin(["Cell ID", "Is complete", "Well"])]] = df.loc[well_mask, df.columns[~df.columns.isin(["Cell ID", "Is complete", "Well"])]].apply(zscore)

features = [f for f in df.columns if "marker" in f]

df_data = df[["Cell size", "Nucleus size", "NCR"]+features].values


scaled_df_data = StandardScaler().fit_transform(df_data)
embedding = reducer.fit_transform(scaled_df_data)


dfp = pd.DataFrame(dict(x=embedding[:,0], y=embedding[:,1], z=df.index))

df["Cluster"] = np.nan


df.loc[(df["01_Average intensity of marker"] >=0)&
(df["02_Average intensity of marker"] >=0)&
(df["03_Average intensity of marker"] >=0), "Cluster"]=1

df.loc[(df["01_Average intensity of marker"]>=0)&(df["02_Average intensity of marker"]<0)&(df["03_Average intensity of marker"]>=0), "Cluster"]=2
df.loc[(df["01_Average intensity of marker"]>=0)&(df["02_Average intensity of marker"]<0)&(df["03_Average intensity of marker"]<0), "Cluster"]=3
df.loc[(df["01_Average intensity of marker"]>=0)&(df["02_Average intensity of marker"]>=0)&(df["03_Average intensity of marker"]<0), "Cluster"]=4
df.loc[(df["01_Average intensity of marker"]<0)&(df["02_Average intensity of marker"]>=0)&(df["03_Average intensity of marker"]>=0), "Cluster"]=5
df.loc[(df["01_Average intensity of marker"]<0)&(df["02_Average intensity of marker"]<0)&(df["03_Average intensity of marker"]>=0), "Cluster"]=6
df.loc[(df["01_Average intensity of marker"]<0)&(df["02_Average intensity of marker"]>=0)&(df["03_Average intensity of marker"]<0), "Cluster"]=7
df.loc[(df["01_Average intensity of marker"]<0)&(df["02_Average intensity of marker"]<0)&(df["03_Average intensity of marker"]<0), "Cluster"]=8



df["Embedding_X"] = embedding[:, 0]
df["Embedding_Y"] = embedding[:,1]
df.to_csv(path+"clustering results.csv")

df["Cluster_color"] = ["red" if l==1 else "silver" for l in df["Cluster"]]

plt.scatter(df.loc[df["Cluster"].isin([1, 2, 5, 6]), "Embedding_X"],
                df.loc[df["Cluster"].isin([1, 2, 5, 6]), "Embedding_Y"],
                c="silver",
                alpha=1,
                s=6)
plt.scatter(df.loc[df["Cluster"].isin([3, 4, 7, 8]), "Embedding_X"],
                df.loc[df["Cluster"].isin([3, 4, 7, 8]), "Embedding_Y"],
                c="blue",
                alpha=0.6,
                s=6)
plt.gca().set_aspect('equal', 'datalim')
plt.title("SOD2-high")
plt.show()
