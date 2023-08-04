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
import scanpy as sc
from genevector.data import GeneVectorDataset
from genevector.model import GeneVector
from genevector.embedding import GeneEmbedding, CellEmbedding
import os
import gdown

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=100, facecolor='white')
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
import torch as t
from genevector.data import GeneVectorDataset
from genevector.model import GeneVector
from genevector.embedding import GeneEmbedding, CellEmbedding
```

```python
#!pip install gdown
```

```python
adata = sc.read("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_Full_genevector_from_top.h5ad")
```

```python
adata = adata[adata.obs["genevector"].isin(["Myeloid+PDC"])]
adata
```

```python
set(adata.obs["genevector"])
```

```python
sc.pp.highly_variable_genes(adata,n_top_genes=2000,flavor="seurat_v3",subset=True)
adata
```

```python
sc.pl.umap(adata,color="genevector")
```

```python
dataset = GeneVectorDataset(adata)
```

```python
#cmps = GeneVector(dataset,
#                  output_file="myeloid.vec",
#                  batch_size=100000000,
#                  c=100.,
#                  emb_dimension=100)
```

```python
#cmps.train(100)
```

```python
embed = GeneEmbedding("myeloid.vec", dataset, vector="average")
```

```python
cembed = CellEmbedding(dataset, embed)
adata = cembed.get_adata()

```

```python
# sc.pl.umap(adata, color="CD14", wspace=0.3, size=5)
```

```python
adata.X[0,:].todense().sum()
```

```python
sc.tl.leiden(adata,resolution=0.1)
```

```python
sc.pl.umap(adata,color="leiden")
```

```python
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
```

```python
sc.tl.rank_genes_groups(adata,"leiden",use_raw=False)
```

```python
sc.tl.dendrogram(adata,"leiden")
sc.pl.rank_genes_groups_dotplot(adata,use_raw=False,groups=["0"], n_genes=20)
```

```python
markers = dict()
markers["pDC"] = ["GZMB","CLIC3","LILRA4"]
markers["Macrophage"] = ["APOE","APOC1"]
markers["Dendritic"] = ["CD1C","CD1A","CD1E","CLEC10A","FCER1A","S100B","CST7","IRF8"]
markers["Monocyte"] = ["VCAN"]
markers["Plasma"] = ["MZB1","JCHAIN"]
markers["B"] = ["CD79B","MS4A1"]
markers["Mast"] = ["KIT","CLU","GATA2","CPA3"]
markers["NK"] = ["GNLY","NCAM1","GZMA","KLRD1"]
```

```python
sc.pl.umap(adata,color=markers["B"],use_raw=False)
```

```python
#sc.pl.umap(full,color="IRF8",use_raw=False)
```

```python

#sc.pl.umap(full,color="leiden",use_raw=False)
```

```python
adata = phenotype_probability(cembed, adata, markers)
```

```python
prob_cols = [x for x in adata.obs.columns.tolist() if "Pseudo-probability" in x]
```

```python
prob_cols
```

```python
sc.pl.umap(adata,color=["genevector"],add_outline=True,s=10,wspace=0.5)
```

```python
# xdata = full[full.obs["leiden"]=="0"]
# sc.pp.highly_variable_genes(xdata,n_top_genes=3000,subset=True,flavor="seurat_v3")
# sc.pp.normalize_total(xdata)
# sc.pp.log1p(xdata)
# sc.tl.pca(xdata)
# sc.pp.neighbors(xdata)
# sc.tl.umap(xdata)
# sc.tl.leiden(xdata,resolution=0.3)
#xdata.obsm["X_umap"] = adata.obsm["X_umap"]
adata
```

```python
#sc.tl.rank_genes_groups(adata,"genevector",use_raw=False)

#dfs = []
#for x in set(adata.obs["genevector"]):
#    df = sc.get.rank_genes_groups_df(adata,x)
#    dfs.append(df)
#df = pandas.concat(dfs)
#df.to_csv("/Users/ceglian/coarse_myeloid_subtypes.csv")
```

```python
#sc.tl.rank_genes_groups(xdata,"leiden",use_raw=False)

#dfs = []
#for x in set(xdata.obs["leiden"]):
#    df = sc.get.rank_genes_groups_df(xdata,x)
#    df["cluster"] = x
#    dfs.append(df)
#df = pandas.concat(dfs)
#df.to_csv("/Users/ceglian/leiden_myeloid.csv")
```

```python
sc.pl.umap(adata,color=["leiden","genevector"])
sc.pl.rank_genes_groups_dotplot(adata,use_raw=False,min_logfoldchange=0.25,swap_axes=True,n_genes=30)
```

```python
sc.pl.dotplot(adata,markers,groupby="genevector",use_raw=False)
```

```python
emarkers = dict()
emarkers["pDC"] = ["GZMB","CLIC3","IRF4","LILRA4"]
emarkers["Macrophage"] = ["APOE","APOC1","C1QA","C1QB","C1QC","CD14","TREM2"]
emarkers["Dendritic"] = ["CD1C","CD1A","CD1E","CLEC10A","FCER1A","S100B","CST7","CCR7"]
emarkers["Monocyte"] = ["VCAN","CD14","S100A8","CD163","S100A9","FCN1"]
emarkers["Plasma"] = ["MZB1","JCHAIN"]
emarkers["B"] = ["CD79B","MS4A1","CD79B","MS4A1"]
emarkers["Mast"] = ["KIT","CLU","GATA2","CPA3"]
emarkers["NK"] = ["GNLY","NCAM1","GZMA","KLRD1"]
sc.pl.dotplot(adata,emarkers,groupby="genevector",use_raw=False)
```

```python
ct = []
for x in adata.obs["genevector"]:
    if x == "Unknown":
        ct.append("Myeloid Undefined")
    else:
        ct.append(x)
adata.obs["genevector"] = ct
```

```python
sc.pl.umap(adata,color=prob_cols)
```

```python
adata.shape
```

```python
#df = adata.obs[["genevector"]]
#df.to_csv("/Users/ceglian/myeloid_subpopulations.csv")
```

```python
full = sc.read("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_Full_genevector_from_top.h5ad")
sp = dict(zip(adata.obs.index,adata.obs["genevector"]))
ct = []
for x,y in zip(full.obs.index,full.obs["genevector"]):
    if x in sp:
        ct.append(sp[x])
    else:
        ct.append(y)
full.obs["genevector_l2"] = ct
```

```python
full.obs
```

```python
set(full.obs["genevector_l2"])
```

```python
di= {"NK": "NK_Cell","B":"B_Cell","Mast":"Mast_Cell"}

full.obs=full.obs.replace({"genevector_l2": di})
```

```python
set(full.obs["genevector_l2"])
```

```python
ct = []
cells = []
for x in set(full.obs["genevector_l2"]):
    cells.append(len(full[full.obs["genevector_l2"]==x].obs.index))
    ct.append(x)
df = pandas.DataFrame.from_dict({"Cells":cells,"Cell Type":ct})
import seaborn as sns
sns.barplot(data=df,x="Cell Type",y="Cells")
plt.xticks(rotation=90)
```

```python
full.write("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_Full_genevector_from_top_with_subpops.h5ad")
```

```python
set(full.obs["genevector_l2"])
```
