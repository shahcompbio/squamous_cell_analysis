---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

```python
!which python3
#!pip install module sacnpy==1.9.1
#%pip install module scanpy
import scanpy as sc
import seaborn as sns
import numpy as np
import scrublet as scr
import scipy.io
```

```python
#Squamous analysis
```

```python
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.manifold import TSNE
from sklearn.preprocessing import MinMaxScaler
import seaborn as sns
import matplotlib.pyplot as plt
import pandas
import umap
import tqdm
import scanpy as sc
import networkx as nx
from sklearn.metrics import confusion_matrix
from scipy.special import softmax
from scipy.spatial import distance
import numpy
from scipy.sparse import csr_matrix, csc_matrix, find
import numpy as np
import operator
import collections
import os  
import pandas as pd
import gc
def normalized_marker_expression(self, normalized_matrix, genes, cells, markers):
    normalized_expression = collections.defaultdict(dict)
    normalized_matrix.eliminate_zeros()
    row_indices, column_indices = normalized_matrix.nonzero()
    nonzero_values = normalized_matrix.data
    entries = list(zip(nonzero_values, row_indices, column_indices))
    for value, i, j in tqdm.tqdm(entries):
        if value > 0 and genes[j].upper() in self.embed.embeddings and genes[j] in markers:
            normalized_expression[cells[i]][genes[j]] = value
    return normalized_expression
class bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
def phenotype_probability(self, adata, phenotype_markers, mark_unknown=True, expression_weighted=True, target_col="genevector"):
    for x in adata.obs.columns:
        if "Pseudo-probability" in x:
            del adata.obs[x]
    mapped_components = dict(zip(list(self.data.keys()),self.matrix))
    genes = adata.var.index.to_list()
    cells = adata.obs.index.to_list()
    matrix = csr_matrix(adata.X)
    embedding = csr_matrix(cembed.matrix)
    all_markers = []
    for _, markers in phenotype_markers.items():
        all_markers += markers
    all_markers = list(set(all_markers))
    normalized_expression = normalized_marker_expression(self, matrix, genes, cells, all_markers)
    probs = dict()
    def generate_weighted_vector(self, genes, weights):
        vector = []
        for gene, vec in zip(self.genes, self.vector):
            if gene in genes and gene in weights:
                vector.append(weights[gene] * numpy.array(vec))
        if numpy.sum(vector) == 0:
            return None
        else:
            return list(numpy.mean(vector, axis=0))
    matrix = matrix.todense()
    for pheno, markers in phenotype_markers.items():
        dists = []
        print(bcolors.OKBLUE+"Computing similarities for {}".format(pheno)+bcolors.ENDC)
        print(bcolors.OKGREEN+"Markers: {}".format(", ".join(markers))+bcolors.ENDC)
        odists = []
        for x in tqdm.tqdm(adata.obs.index):
            weights = normalized_expression[x]
            vector = generate_weighted_vector(self.embed, markers, weights)
            if vector != None:
                dist = 1. - distance.cosine(mapped_components[x], numpy.array(vector))
                odists.append(dist)
            else:
                odists.append(0.)
        probs[pheno] = odists
    distribution = []
    celltypes = []
    for k, v in probs.items():
        distribution.append(v)
        celltypes.append(k)
    distribution = list(zip(*distribution))
    probabilities = softmax(numpy.array(distribution),axis=1)
    res = {"distances":distribution, "order":celltypes, "probabilities":probabilities}
    barcode_to_label = dict(zip(list(self.data.keys()), res["probabilities"]))
    ct = []
    probs = collections.defaultdict(list)
    for x in adata.obs.index:
        ctx = res["order"][numpy.argmax(barcode_to_label[x])]
        ct.append(ctx)
        for ph, pb in zip(res["order"],barcode_to_label[x]):
            probs[ph].append(pb)
    adata.obs[target_col] = ct
    def load_predictions(adata,probs):
        for ph in probs.keys():
            adata.obs[ph+" Pseudo-probability"] = probs[ph]
        return adata
    adata = load_predictions(adata,probs)
    prob_cols = [x for x in adata.obs.columns if "Pseudo" in x]
    cts = adata.obs[target_col].tolist()
    probs = adata.obs[prob_cols].to_numpy()
    adj_cts = []
    for ct,p in zip(cts,probs):
        if len(set(p)) == 1:
            adj_cts.append("Unknown")
        else:
            adj_cts.append(ct)
    adata.obs[target_col] = adj_cts
    return adata
```

```python
adata = sc.read("/work/shah/users/Mobina/Project_Squamous/Squamousdata/Final_data.h5ad")
```

```python
adata.obs
```

```python
adata.obs["orig.ident"].value_counts()
```

```python
5253928-909277
```

```python
adata=adata[adata.obs["orig.ident"]!="3590_YL-1540_1200a_IGO_12437_AB_11",]
```

```python
adata.shape
```

```python
sc.pp.filter_cells(adata, min_genes=200)
```

```python
adata.var['mt'] = adata.var_names.str.startswith('MT-') 
```

```python
adata.shape
```

```python
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
```

```python
adata.obs["index"]=adata.obs.index.tolist()
```

```python
Data_doublet=[]
df_doublet=[]
for ph in np.unique(adata.obs["orig.ident"]):
 Data=adata[adata.obs["orig.ident"]==ph]    
 scrub = scr.Scrublet(Data.X)
 doublet_scores, predicted_doublets = scrub.scrub_doublets()
 scrub.plot_histogram()
 print(ph)
 df = pd.DataFrame({'score':doublet_scores, 'doublet':predicted_doublets})
 df["orig.ident"]=Data.obs["orig.ident"].tolist()
 df["0.4threshold"]= df["score"]<0.4
 df_doublet.append(df)  
 print(np.unique(df["0.4threshold"], return_counts=True))
 Data=Data[df["0.4threshold"],:]
 Data_doublet.append(Data.obs)

```

```python
Data_doublet_concat= pd.concat(Data_doublet, ignore_index=True)
df_doublet_concat=pd.concat(df_doublet, ignore_index=True)
```

```python
df_doublet_concat
```

```python
Data_doublet_concat.index=Data_doublet_concat["index"]
#del Data_doublet_concat['index']
```

```python
#len(Data_doublet_concat["index"])
```

```python
adata=adata[Data_doublet_concat.index.tolist(),]
```

```python
adata.shape
```

```python
#Solve the problem for adata.write (ref:https://github.com/theislab/scvelo/issues/255)
#adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
```

```python
#adata.write("/work/shah/users/Mobina/Project_Squamous/Squamousdata/Final_data_1200_cells_doublet_removed.h5ad")
```

```python
#del adata_2.obs['index']
```

```python
sc.pl.violin(adata, ['n_genes_by_counts'],groupby="orig.ident",
             jitter=False, multi_panel=True,rotation=90)

```

```python
sc.pl.violin(adata, ['total_counts'],groupby="orig.ident",
             jitter=False, multi_panel=True,rotation=90)

```

```python
sc.pl.violin(adata, ['pct_counts_mt'],groupby="orig.ident",
             jitter=False, multi_panel=True,rotation=90)

```

```python
#set(adata.obs["orig.ident"])
```

```python
#adata.obs.groupby('orig.ident').plot(kind='scatter', x='total_counts', y='n_genes_by_counts',subplots=True,s=1,cmap= "tab10",c="genevector")
```

```python
adata.obs.groupby('orig.ident').plot(kind='scatter', x='n_genes_by_counts', y='pct_counts_mt',subplots=True)
```

```python
adata.obs
```

```python
#adata.obs.groupby('orig.ident').plot(kind='scatter', x='total_counts', y='n_genes_by_counts',subplots=True,s=1,cmap= "tab10",c="genevector")
```

```python
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
```

```python
ct = []
cells = []
for x in set(adata.obs["orig.ident"]):
    cells.append(len(adata[adata.obs["orig.ident"]==x].obs.index))
    ct.append(x)
df = pandas.DataFrame.from_dict({"Patient":cells,"Cell Count":ct})
import seaborn as sns
sns.barplot(data=df,x="Patient",y="Cell Count")
plt.xticks(rotation=90)
plt.show()
```

```python
import torch as t
from genevector.data import GeneVectorDataset
from genevector.model import GeneVector
from genevector.embedding import GeneEmbedding, CellEmbedding
```

```python
data_wide = df_doublet_concat.pivot(columns='orig.ident',
                     values='score')
data_wide.head()
data_wide.plot.density(figsize = (7, 7),
                       linewidth = 1)
```

```python
df_mt=adata.obs[["orig.ident","pct_counts_mt"]]
```

```python
data_wide = df_mt.pivot(columns='orig.ident',
                     values='pct_counts_mt')
data_wide.head()
data_wide.plot.density(figsize = (7, 7),
                       linewidth = 1)
```

```python
df_doublet_concat
```

```python
adata = adata[adata.obs.pct_counts_mt <35, :]
```

```python
adata.var_names_make_unique()
```

```python
print(adata)
```

```python
genes = [x for x in adata.var.index if "." not in x]
genes = [x for x in genes if "-" not in x or "HLA" in x]
genes = [x for x in genes if not x.startswith("MT-")]
genes = [x for x in genes if not x.startswith("RPL")]
genes = [x for x in genes if not x.startswith("RPS")]
adata = adata[:,genes]
print(adata)
```

```python
adata_QC=adata.copy()
```

```python
sc.pp.normalize_total(adata, target_sum=1e4)
```

```python
sc.pp.log1p(adata)
```

```python
sc.pp.calculate_qc_metrics(adata,inplace=True)
```

```python
# adata[:,adata.var.index.isin(["CD3D","LYZ","CD3E","CD3G","CD4"])].var

```

```python
#flavor="seurat_v3"
#subset=True
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000,batch_key="orig.ident")
```

```python
adata.obs
```

```python
for x in adata.var.index.tolist():
    if "CD3E" in x:
        print(x)
```

```python

```

```python
#adata=adata[:, adata.var['mt'] == False]
```

```python
"CD3E" in adata.var.index.tolist()
```

```python
VariableGenes=adata.var.highly_variable
```

```python
VariableGenes
```

```python
"CD3E" in VariableGenes
```

```python
#adata = sc.read("../Project_Squamous/Squamousdata/seurat_data_withQC_1200removed.h5ad")
```

```python
#adata.shape
```

```python
#adata.var_names_make_unique()
```

```python
#genes = [x for x in adata.var.index if "." not in x]
#genes = [x for x in genes if "-" not in x or "HLA" in x]
#genes = [x for x in genes if not x.startswith("MT-")]
#genes = [x for x in genes if not x.startswith("RPL")]
#genes = [x for x in genes if not x.startswith("RPS")]
#adata = adata[:,genes]
#print(adata)
```

```python
adata = adata_QC[:, VariableGenes]
```

```python
adata.shape
```

```python
"CD1C" in adata.var.index.tolist()
```

```python
adata_QC
```

```python
adata
```

```python
#adata.X = adata.X.todense()
adata_original=adata.copy()
```

```python
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
#adata.write("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_genevector_input.h5ad")
adata.write("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_genevector_input_new.h5ad")
```

```python
adata_QC.__dict__['_raw'].__dict__['_var'] = adata_QC.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
#adata_QC.write("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_QC_genevector_input.h5ad")
adata_QC.write("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_QC_genevector_input_new.h5ad")
```

```python

```
