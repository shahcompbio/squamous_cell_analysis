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
adata = sc.read("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_genevector_output.h5ad")
adata_QC = sc.read("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_QC_genevector_input.h5ad")
```

```python
adata_QC.obs
```

```python
adata.obs
```

```python
adata_Full = adata_QC.copy()
sc.pp.normalize_total(adata_Full, target_sum=1e4)
sc.pp.log1p(adata_Full)
sc.pp.calculate_qc_metrics(adata_Full,inplace=True)
adata_Full.obs=adata.obs
adata_Full.obsm=adata.obsm
```

```python
adata_Full.obs
```

```python
genes = [x for x in adata_Full.var.index if "." not in x]
genes = [x for x in genes if "-" not in x or "HLA" in x]
genes = [x for x in genes if not x.startswith("MT-")]
genes = [x for x in genes if not x.startswith("RPL")]
genes = [x for x in genes if not x.startswith("RPS")]
adata_Full = adata_Full[:,genes]
```

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

```python
adata
```

```python
adata_Full
```

```python
adata_QC
```

```python
markers = dict()
markers["T Cell"] = ["CD3D","IL32"]
markers["B Cell"] = ["CD79A","CD19"]
markers["Plasma"] = ["MZB1","JCHAIN"]
markers["Myeloid"] = ["LYZ","FCER1G","AIF1"]
markers["Cancer"] = ["KRT19","KRT8","KRT18"]
markers["Fibroblast"] = ["COL1A1","COL3A1","SPARC"]
#markers["Mast Cell"] = ["KIT","CLU"]
```

```python
markers.keys()
```

```python
sc.pl.umap(adata_Full,color=["EPCAM","PTPRC","genevector"],use_raw=False)
```

```python
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.rc_file_defaults()
sns.set(font_scale=1.)
def celltype_composition(adata_tmp, ax, key, cat, order, title="", legend=True):
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
adata.obs["index"]=adata.obs.index.tolist()
```

```python
adata_Full.obs["index"]=adata_Full.obs.index.tolist()
```

```python
adata_Full.obs
```

```python
import scanpy.external as sce
```

```python
Data = adata_Full[adata_Full.obs["genevector"]=="Myeloid"]
sc.pp.highly_variable_genes(Data, flavor="seurat_v3",subset=True, n_top_genes=2500,batch_key="orig.ident")
sc.tl.pca(Data, svd_solver='arpack')
sce.pp.harmony_integrate(Data, "orig.ident")
sc.pp.neighbors(Data, n_neighbors=10, n_pcs=40)
sc.tl.umap(Data)
sc.tl.leiden(Data,resolution=0.1)
sc.pl.umap(Data,color=["leiden"],use_raw=False,title="Myeloid")

```

```python
sc.pl.umap(Data,color=["CD79A","CD79B","CD3E"],use_raw=False,ncols=3,vmax=1)
```

```python
Data.uns["leiden_colors"]
```

```python
sc.pl.embedding(Data, basis='X_genevector')
```

```python
#sc.pl.umap(Data,color=["leiden","orig.ident"],use_raw=False,title=ph)
```

```python
#same analysis for all the celltypes
Data_celltype=[]
Data_celltype_adata=[]
for ph,genes in markers.items():
  Data = adata_Full[adata_Full.obs["genevector"]==ph]
  sc.pp.highly_variable_genes(Data, flavor="seurat_v3",subset=True, n_top_genes=2500,batch_key="orig.ident")
  sc.tl.pca(Data, svd_solver='arpack')
  sce.pp.harmony_integrate(Data, "orig.ident")
  sc.pp.neighbors(Data, n_neighbors=10, n_pcs=40)
  sc.tl.umap(Data)
  sc.tl.leiden(Data,resolution=0.1)
  sc.tl.rank_genes_groups(Data, 'leiden', method='t-test', use_raw=False)
  sc.pl.rank_genes_groups_dotplot(Data,use_raw=False,min_logfoldchange=1.,title=ph)
  sc.pl.umap(Data,color=["leiden","orig.ident"],use_raw=False,title=ph)
  fig, axes = plt.subplots(1, 1, figsize=(5,3), gridspec_kw={'wspace':0.2})
  celltype_composition(Data, axes, "orig.ident","leiden", 
                     list(set(Data.obs["orig.ident"])), legend=False)
  #sc.pl.umap(Data,color=["orig.ident"],use_raw=False,title=ph)  
  Data_celltype.append(Data.obs)
  Data_celltype_adata.append(Data)
  print(ph)
```

```python
Data = adata_Full[adata_Full.obs["genevector"]=="Unknown"]
sc.pp.highly_variable_genes(Data, flavor="seurat_v3",subset=True, n_top_genes=2500,batch_key="orig.ident")
sc.tl.pca(Data, svd_solver='arpack')
sce.pp.harmony_integrate(Data, "orig.ident")
sc.pp.neighbors(Data, n_neighbors=10, n_pcs=40)
sc.tl.umap(Data)
sc.tl.leiden(Data,resolution=0.1)
sc.tl.rank_genes_groups(Data, 'leiden', method='t-test', use_raw=False)
sc.pl.rank_genes_groups_dotplot(Data,use_raw=False,min_logfoldchange=1.,title="Unknown")
sc.pl.umap(Data,color=["leiden","orig.ident"],use_raw=False,title="Unknown")
fig, axes = plt.subplots(1, 1, figsize=(5,3), gridspec_kw={'wspace':0.2})
celltype_composition(Data, axes, "orig.ident","leiden", 
                     list(set(Data.obs["orig.ident"])), legend=False)
Data_celltype.append(Data.obs)
Data_celltype_adata.append(Data)
```

```python
sc.pl.umap(adata,color=["genevector"],use_raw=False)
```

```python
Data_celltype[4]
```

```python
#sns.violinplot(x="orig.ident",y="pct_counts_mt",data=adata_Full.obs)
#plt.xticks(rotation=90)
#plt.show()
```

```python
#sns.violinplot(x="orig.ident",y="pct_counts_mt",data=Data_celltype[0].obs)
#plt.xticks(rotation=90)
#plt.show()

```

```python
for i in range(7):
    print(set(Data_celltype[i]["leiden"]))
```

```python
markers.keys()
```

```python
di_Tcell= {"0": "0_Tcell", "1": "1_Tcell", "2": "2_Tcell"}
di_Bcell= {"0": "0_Bcell" ,"1": "1_Bcell", "2": "2_Bcell"}
di_Plasma = {"0": "0_Plasma", "1": "1_Plasma"}
di_Myeloid= {"0":"0_Myeloid","1": "1_Myeloid", "2": "2_Myeloid", "3":"3_Myeloid","4":"4_Myeloid","5":"5_Myeloid"}
di_Fibroblast= {"0":"0_Fibroblast","1": "1_Fibroblast", "2": "2_Fibroblast", "3":"3_Fibroblast"}
di_Unknown= {"0":"0_Unknown","1": "1_Unknown", "2": "2_Unknown"}
```

```python
Data_celltype[0]=Data_celltype[0].replace({"leiden": di_Tcell})
Data_celltype[1]=Data_celltype[1].replace({"leiden": di_Bcell})
Data_celltype[2]=Data_celltype[2].replace({"leiden": di_Plasma})
Data_celltype[3]=Data_celltype[3].replace({"leiden": di_Myeloid})
Data_celltype[5]=Data_celltype[5].replace({"leiden": di_Fibroblast})
Data_celltype[6]=Data_celltype[6].replace({"leiden": di_Unknown})
```

```python
Data_celltype[0]
```

```python
Data_celltype_concat= pd.concat(Data_celltype, ignore_index=True)
```

```python
Data_celltype_concat.index=Data_celltype_concat["index"]
```

```python
Data_celltype_concat
```

```python
adata_Full.obs=Data_celltype_concat
```

```python
adata_Full
```

```python
#adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
adata_Full.write("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_Full_subtype_output.h5ad")
```

```python
#adata_Full = sc.read("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_Full_subtype_output.h5ad")
```

```python

```

```python

```

```python

```

```python

```

```python

```

```python

```

```python

```

```python

```

```python

```

```python

```

```python

```

```python

```

```python

```

```python

```
