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
import scanpy.external as sce
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
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.rc_file_defaults()
sns.set(font_scale=1.)
def celltype_composition(adata_tmp, ax, key, cat, order,title="", legend=True):
    sizes = adata_tmp.obs.groupby([cat, key]).size()
    props = sizes.groupby(level=1).apply(lambda x: 1. * x / x.sum()).reset_index() 
    props = props.pivot(columns=key, index=cat).T
    props.index = props.index.droplevel(0)
    props.fillna(0, inplace=True)
    props.reindex(order).plot(
        kind="bar", 
        stacked=True, 
        ax=ax, 
        legend=None,
        colormap=plt.get_cmap("tab20"),
    )
    ax.legend(bbox_to_anchor=(1.01, 1), frameon=False, title=title)
    sns.despine(fig, ax)
    ax.tick_params(axis="x", rotation=90)
    ax.set_xlabel("")
    ax.set_ylabel("")
```

```python
def barplot(adata,x_lab,y_name):
    #,color,order
    ct = []
    cells = []
    for x in set(adata.obs[x_lab]):
        cells.append(len(adata[adata.obs[x_lab]==x].obs.index))
        ct.append(x)
        df = pandas.DataFrame.from_dict({"Count":cells,y_name:ct})
    import seaborn as sns
    sns.barplot(data=df,x="Count",y= y_name)
    #,palette=color,order=order
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
#Data_myeloid_liftover
```

```python
adata_myeloid_liftover = sc.read("/work/shah/users/Mobina/Project_Squamous/Squamousdata/Data_myeloid_liftover.h5ad")
```

```python
adata_myeloid_liftover.X[0,:].todense().sum()
```

```python
adata_myeloid_liftover
```

```python
fig, axes = plt.subplots(1, 1, figsize=(5,3), gridspec_kw={'wspace':0.2})
celltype_composition(adata_myeloid_liftover, axes,"predicted.ann_level_3","genevector_updated", 
                     list(set(adata_myeloid_liftover.obs["predicted.ann_level_3"])), legend=False)
```

```python
fig, axes = plt.subplots(1, 1, figsize=(5,3), gridspec_kw={'wspace':0.2})
celltype_composition(adata_myeloid_liftover, axes,"predicted.ann_level_3","orig.ident", 
                     list(set(adata_myeloid_liftover.obs["predicted.ann_level_3"])), legend=False)
```

```python
barplot(adata_myeloid_liftover,"predicted.ann_level_3","Cell Types")
```

```python
genes=["CSF1R","FCGR1A","FCGR1B","MAFB","CD68","ITGAM","CD14","MRC1","TREM2","APOE","VSIG4","LILRB1","LILRB2","LILRB3","LILRB4","LILRB5","ITGAX","LAMP3","CCR7","XCR1","IRF8","CD1C","CD163","CD3D"]
```

```python
adata_myeloid_liftover_monocyte=adata_myeloid_liftover[(adata_myeloid_liftover.obs["predicted.ann_level_3"]=="Monocytes")]
```

```python
adata_myeloid_liftover_macrophage=adata_myeloid_liftover[(adata_myeloid_liftover.obs["predicted.ann_level_3"]=="Macrophages")]
```

```python
adata_myeloid_liftover_DC=adata_myeloid_liftover[(adata_myeloid_liftover.obs["predicted.ann_level_3"]=="Dendritic cells")]
```

```python
set(adata_myeloid_liftover_DC.obs["predicted.ann_level_3"])
```

```python
#Monocytes
```

```python
sc.tl.pca(adata_myeloid_liftover_monocyte, svd_solver='arpack')
sc.pp.neighbors(adata_myeloid_liftover_monocyte, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_myeloid_liftover_monocyte)
```

```python
sc.pl.umap(adata_myeloid_liftover_monocyte,color=["predicted.ann_level_3","genevector_updated"],use_raw=False,ncols=2)
```

```python
barplot(adata_myeloid_liftover_monocyte,"genevector_updated","Cell Types")
```

```python
adata_myeloid_liftover_monocyte=adata_myeloid_liftover_monocyte[adata_myeloid_liftover_monocyte.obs["genevector_updated"].isin(["Monocytes","T_Cell","Epithelial"]),]
```

```python
set(adata_myeloid_liftover_monocyte.obs["genevector_updated"])
```

```python
adata_myeloid_liftover_monocyte.uns['log1p']["base"] = None
sc.tl.rank_genes_groups(adata_myeloid_liftover_monocyte, 'genevector_updated', method='t-test', use_raw=False)
sc.pl.rank_genes_groups_dotplot(adata_myeloid_liftover_monocyte,use_raw=False,min_logfoldchange=0.25,title="Myeloids")
```

```python
genes=["CD3D","EPCAM","S100A8","TREM2","APOE","TYROBP","LYZ","FCER1G","AIF1"]
for i in range(0,23,3):
    sc.pl.violin(adata_myeloid_liftover_monocyte,[genes[i],genes[(i+1)],genes[(i+2)]],groupby="genevector_updated",use_raw=False,rotation=90) 
#sc.pl.umap(adata_myeloid_liftover_monocyte,color=["CD3D","EPCAM","S100A8","TREM2","APOE","TYROBP","LYZ","FCER1G","AIF1"],use_raw=False,ncols=3)
```

```python
#sc.pl.scatter(adata_myeloid_liftover_monocyte, x='total_counts', y='pct_counts_mt')
#sc.pl.scatter(adata_myeloid_liftover_monocyte, x='total_counts', y='n_genes_by_counts')
adata_myeloid_liftover_monocyte.obs.plot(kind='scatter', x='total_counts', y='n_genes_by_counts',subplots=True,s=1,cmap="coolwarm",c="genevector_updated")
```

```python
#Macrophage

sc.tl.pca(adata_myeloid_liftover_macrophage, svd_solver='arpack')
sc.pp.neighbors(adata_myeloid_liftover_macrophage, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_myeloid_liftover_macrophage)
sc.pl.umap(adata_myeloid_liftover_macrophage,color=["predicted.ann_level_3","genevector_updated"],use_raw=False,ncols=1)
```

```python
barplot(adata_myeloid_liftover_macrophage,"genevector_updated","Cell Types")
```

```python
adata_myeloid_liftover_macrophage=adata_myeloid_liftover_macrophage[adata_myeloid_liftover_macrophage.obs["genevector_updated"].isin(["Macrophages","T_Cell","Epithelial","Cancer","Unknown"]),]
```

```python
set(adata_myeloid_liftover_macrophage.obs["genevector_updated"])
```

```python
adata_myeloid_liftover_macrophage.uns['log1p']["base"] = None
sc.tl.rank_genes_groups(adata_myeloid_liftover_macrophage, 'genevector_updated', method='t-test', use_raw=False)
sc.pl.rank_genes_groups_dotplot(adata_myeloid_liftover_macrophage,use_raw=False,min_logfoldchange=0.25,title="Myeloids")
```

```python
genes=["CD3D","EPCAM","S100A8","TREM2","APOE","TYROBP","LYZ","FCER1G","AIF1"]
for i in range(0,23,3):
    sc.pl.violin(adata_myeloid_liftover_macrophage,[genes[i],genes[(i+1)],genes[(i+2)]],groupby="genevector_updated",use_raw=False,rotation=90) 
#sc.pl.umap(adata_myeloid_liftover_monocyte,color=["CD3D","EPCAM","S100A8","TREM2","APOE","TYROBP","LYZ","FCER1G","AIF1"],use_raw=False,ncols=3)
```

```python
#sc.pl.scatter(adata_myeloid_liftover_monocyte, x='total_counts', y='pct_counts_mt')
#sc.pl.scatter(adata_myeloid_liftover_monocyte, x='total_counts', y='n_genes_by_counts')
adata_myeloid_liftover_macrophage.obs.plot(kind='scatter', x='total_counts', y='n_genes_by_counts',subplots=True,s=1,cmap="coolwarm",c="genevector_updated")
```

```python
#DC

sc.tl.pca(adata_myeloid_liftover_DC, svd_solver='arpack')
sc.pp.neighbors(adata_myeloid_liftover_DC, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_myeloid_liftover_DC)
sc.pl.umap(adata_myeloid_liftover_DC,color=["predicted.ann_level_3","genevector_updated"],use_raw=False,ncols=2)
```

```python
barplot(adata_myeloid_liftover_DC,"genevector_updated","Cell Types")
```

```python
adata_myeloid_liftover_DC=adata_myeloid_liftover_DC[adata_myeloid_liftover_DC.obs["genevector_updated"].isin(["Dendritic_cells","T_Cell","Epithelial","Unknown","Plasmacytoid_DCs"]),]
```

```python
set(adata_myeloid_liftover_DC.obs["genevector_updated"])
```

```python
adata_myeloid_liftover_DC.uns['log1p']["base"] = None
sc.tl.rank_genes_groups(adata_myeloid_liftover_DC, 'genevector_updated', method='t-test', use_raw=False)
sc.pl.rank_genes_groups_dotplot(adata_myeloid_liftover_DC,use_raw=False,min_logfoldchange=0.25,title="Myeloids")
```

```python
genes=["CD3D","EPCAM","S100A8","TREM2","APOE","TYROBP","LYZ","FCER1G","AIF1"]
for i in range(0,23,3):
    sc.pl.violin(adata_myeloid_liftover_DC,[genes[i],genes[(i+1)],genes[(i+2)]],groupby="genevector_updated",use_raw=False,rotation=90) 
#sc.pl.umap(adata_myeloid_liftover_monocyte,color=["CD3D","EPCAM","S100A8","TREM2","APOE","TYROBP","LYZ","FCER1G","AIF1"],use_raw=False,ncols=3)
```

```python
#sc.pl.scatter(adata_myeloid_liftover_monocyte, x='total_counts', y='pct_counts_mt')
#sc.pl.scatter(adata_myeloid_liftover_monocyte, x='total_counts', y='n_genes_by_counts')
adata_myeloid_liftover_DC.obs.plot(kind='scatter', x='total_counts', y='n_genes_by_counts',subplots=True,s=1,cmap="coolwarm",c="genevector_updated")
```

```python

```
