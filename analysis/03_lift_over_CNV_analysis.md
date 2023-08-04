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
def barplot(adata,x_lab,y_name,color,order):
    ct = []
    cells = []
    for x in set(adata.obs[x_lab]):
        cells.append(len(adata[adata.obs[x_lab]==x].obs.index))
        ct.append(x)
        df = pandas.DataFrame.from_dict({"Count":cells,y_name:ct})
    import seaborn as sns
    sns.barplot(data=df,x="Count",y= y_name,palette=color,order=order)
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
#adata_Full = sc.read("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_Full_liftover_output.h5ad")
adata_Full = sc.read("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_Full_subtype_output_filtered_with_Unknown.h5ad")

```

```python
#sc.pp.log1p(adata_Full) it's already log normalized
adata_Full.uns['log1p']["base"] = None
```

```python
adata_Full.obs
```

```python
adata_Full.obs.columns.tolist()
```

```python
sc.pl.umap(adata_Full,color=["genevector"],use_raw=False)
```

```python
sc.pl.umap(adata_Full,color=["genevector"],use_raw=False,ncols=1)
sc.pl.umap(adata_Full,color=["genevector_updated"],use_raw=False,ncols=1)
sc.pl.umap(adata_Full,color=["predicted.ann_level_3"],use_raw=False,ncols=1)
sc.pl.umap(adata_Full,color=["predicted.ann_level_4"],use_raw=False,ncols=1)
sc.pl.umap(adata_Full,color=["leiden"],use_raw=False,ncols=1)
```

```python
a=adata_Full[(adata_Full.obs["genevector"]=="Cancer")]
sc.pl.umap(a,color=["genevector_updated"],use_raw=False,ncols=1)
sc.pl.umap(a,color=["predicted.ann_level_3"],use_raw=False,ncols=1)
```

```python
a=adata_Full[(adata_Full.obs["leiden"]=="Epithelial")]
```

```python
sc.pl.umap(a,color=["leiden"],use_raw=False)
```

```python
a=adata_Full[(adata_Full.obs["genevector"]=="Cancer")]
sc.pl.umap(a,color=["leiden"],use_raw=False)
a=adata_Full[(adata_Full.obs["genevector"]=="Myeloid")]
sc.pl.umap(a,color=["leiden"],use_raw=False)
```

```python
sc.pl.umap(a,color=["predicted.ann_level_3"],use_raw=False)
```

```python
Unknown=adata_Full[(adata_Full.obs["genevector"]=="Unknown")]
```

```python
sc.pl.umap(Unknown,color=["leiden"],use_raw=False,ncols=1,palette="tab10")
sc.pl.umap(Unknown,color=["predicted.ann_level_3"],use_raw=False,ncols=1,palette="tab10")
sc.pl.umap(Unknown,color=["leiden","predicted.ann_level_2","predicted.ann_level_3","predicted.ann_level_4"],use_raw=False,ncols=1,palette="tab10")
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
#Myeloid group
```

```python
myeloid_genevector=adata_Full[(adata_Full.obs["genevector"]=="Myeloid")]
myeloid_liftover=adata_Full[(adata_Full.obs["predicted.ann_level_2"]=="Myeloid")]
```

```python
myeloid_genevector
```

```python
myeloid_liftover
```

```python
sc.pl.umap(myeloid_genevector,color=["predicted.ann_level_3"],use_raw=False,ncols=1,palette="tab10")
sc.tl.pca(myeloid_genevector, svd_solver='arpack')
sce.pp.harmony_integrate(myeloid_genevector, "orig.ident")
sc.pp.neighbors(myeloid_genevector, n_neighbors=10, n_pcs=40)
sc.tl.umap(myeloid_genevector)
sc.pl.umap(myeloid_genevector,color=["predicted.ann_level_3"],use_raw=False)
```

```python
sc.pl.umap(myeloid_genevector,color=["predicted.ann_level_3"],use_raw=False,palette=myeloid_genevector.uns["predicted.ann_level_3_colors"])
```

```python
set(myeloid_genevector.obs["predicted.ann_level_3"])
```

```python
order=['B cell lineage',
 'Dendritic cells',
 'Fibroblasts',
 'Innate lymphoid cell NK',
 'Macrophages',
 'Mast cells',
 'Monocytes',
 'T cell lineage']
```

```python
barplot(myeloid_genevector,"predicted.ann_level_3","lift_over_subgroups",myeloid_genevector.uns["predicted.ann_level_3_colors"],order=order)
```

```python
list(set(myeloid_genevector.obs["leiden"]))
```

```python
fig, axes = plt.subplots(1, 1, figsize=(5,3), gridspec_kw={'wspace':0.2})
celltype_composition(myeloid_genevector , axes, "leiden","predicted.ann_level_3", 
                     list(set(myeloid_genevector.obs["leiden"])),legend=False)
```

```python
fig, axes = plt.subplots(1, 1, figsize=(5,3), gridspec_kw={'wspace':0.2})
celltype_composition(myeloid_genevector , axes, "leiden","predicted.ann_level_4", 
                     list(set(myeloid_genevector.obs["leiden"])), legend=False)
```

```python
fig, axes = plt.subplots(1, 1, figsize=(5,3), gridspec_kw={'wspace':0.2})
celltype_composition(myeloid_genevector , axes, "leiden","predicted.ann_level_2", 
                     list(set(myeloid_genevector.obs["leiden"])), legend=False)
```

```python
# look at Plasmacytoid DCs
```

```python
Plasmacytoid_DCs=adata_Full[(adata_Full.obs["predicted.ann_level_4"]=="Plasmacytoid DCs")]
```

```python
sc.tl.pca(Plasmacytoid_DCs, svd_solver='arpack')
sc.pp.neighbors(Plasmacytoid_DCs, n_neighbors=10, n_pcs=40)
sc.tl.umap(Plasmacytoid_DCs)
sc.pl.umap(Plasmacytoid_DCs,color=["leiden"],use_raw=False,palette="tab10")
#barplot(Plasmacytoid_DCs,"leiden","leiden_subgroups")
```

```python
adata_Full = sc.read("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_Full_subtype_output_filtered.h5ad")
```

```python
adata_Full.uns['log1p']["base"] = None
```

```python
adata_Full
```

```python
set(adata_Full.obs["genevector_updated"])
```

```python
sc.pl.umap(adata_Full,color=["genevector_updated"],use_raw=False,ncols=1,palette="tab20")
```

```python
adata_Full
```

```python
sc.tl.pca(adata_Full, svd_solver='arpack')
sce.pp.harmony_integrate(adata_Full, "orig.ident")
sc.pp.neighbors(adata_Full, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_Full)
sc.pl.umap(adata_Full,color=["predicted.ann_level_3"],use_raw=False)
```

```python
sc.pl.umap(adata_Full,color=["genevector_updated"],use_raw=False)
```

```python
set(adata_Full.obs["genevector_updated"])
```

```python
order=['B_Cell',
 'Cancer',
 'Dendritic_cells',
 'Epithelial',
 'Fibroblast',
 'Innate_lymphoid_cell_NK',
 'Macrophages',
 'Mast_Cell',
 'Monocytes',
 'Plasma',
 'Plasmacytoid_DCs',
 'T_Cell']
```

```python
barplot(adata_Full,"genevector_updated","Cell Types",adata_Full.uns["genevector_updated_colors"],order=order)
```

```python
#sc.tl.dendrogram(adata_Full, groupby='leiden')
sc.tl.rank_genes_groups(adata_Full,'genevector_updated', method='t-test', use_raw=False)
```

```python
for i in set(adata_Full.obs["genevector_updated"]):
    sc.pl.rank_genes_groups_dotplot(adata_Full,use_raw=False,min_logfoldchange=1.,groupby="genevector_updated",groups=[i])
```

```python
sc.pl.umap(adata_Full,color=["genevector_updated","KLRD1","CD3D","ITGAM","PSAP","JCHAIN","IGKC","KRT19","GZMB"],use_raw=False,ncols=3,vmax=1)
```

```python
"ITGAM" in adata_Full.var.index.tolist()
```

```python
sc.pl.umap(adata_Full,color=["genevector_updated","CSF1R","FCGR1A","FCGR1B","MAFB","CD68","ITGAM","CD14","MRC1","TREM2","APOE","VSIG4","LILRB1","LILRB2","LILRB3","LILRB4","LILRB5","ITGAX","LAMP3","CCR7","XCR1","IRF8","CD1C","CD163"],use_raw=False,ncols=3,vmax=1)
```

```python
genes=["CSF1R","FCGR1A","FCGR1B","MAFB","CD68","ITGAM","CD14","MRC1","TREM2","APOE","VSIG4","LILRB1","LILRB2","LILRB3","LILRB4","LILRB5","ITGAX","LAMP3","CCR7","XCR1","IRF8","CD1C","CD163","CD3D"]
```

```python
for i in range(0,23,2):
    sc.pl.violin(adata_Full,[genes[i],genes[(i+1)]],groupby="genevector_updated",use_raw=False,rotation=90) 
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
#Cellphonedb
adata_Full.obs["celltype"]=adata_Full.obs["genevector_updated"]
adata_Full.write("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_Full_cellphonedb.h5ad")
```

```python
#df_meta = pd.DataFrame(data={'Cell':list(adata_Full.obs.index),
#                             'cell_type':[ i for i in adata_Full.obs['celltype']]
#                            })
```

```python
#df_meta.set_index('Cell', inplace=True)
```

```python
#df_meta.to_csv('/work/shah/users/Mobina/Project_Squamous/Squamousdata/squamous_metadata.tsv', sep = '\t')
```

```python
#df_meta
```

```python

```
