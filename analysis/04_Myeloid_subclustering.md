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
#Myeloid group
```

```python
adata_Full = sc.read("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_Full_subtype_output_filtered.h5ad")
```

```python
adata_Full.X[0,:].todense().sum()
```

```python
adata_Full
```

```python
sc.pl.violin(adata_Full,["IGKC","IGHG2","IGHG1","MZB1","JCHAIN","GZMB"],groupby="genevector_updated",use_raw=False,rotation=90) 
```

```python
sc.pl.violin(adata_Full,["MS4A1","XBP1","IL3RA"],groupby="genevector_updated",use_raw=False,rotation=90) 
```

```python
sc.pl.violin(adata_Full,["IRF7","IRF8","CD4","CLEC4C","PTGDS","LILRA4"],groupby="genevector_updated",use_raw=False,rotation=90) 
```

```python
sc.pl.violin(adata_Full,["FCER1A","GPX1","RGS1","CLEC10A","S100B","LYZ","AIF1","FCER1G","IL3RA","TRBV2","CAP2","FOXH1","CD1C"],groupby="genevector_updated",use_raw=False,rotation=90) 
```

```python
sc.pl.violin(adata_Full,["FCER1A","TREM2","CD163","FCER1G","LYZ","AIF1"],groupby="genevector_updated",use_raw=False,rotation=90) 
```

```python
ct = []
cells = []
for x in set(adata_Full.obs["genevector_updated"]):
    cells.append(len(adata_Full[adata_Full.obs["genevector_updated"]==x].obs.index))
    ct.append(x)
df = pandas.DataFrame.from_dict({"Cells":cells,"Cell Type":ct})
import seaborn as sns
sns.barplot(data=df,x="Cell Type",y="Cells")
plt.xticks(rotation=90)
```

```python
for x in adata_Full.var.index.tolist():
    if "TRBV" in x:
        print(x)
```

```python
genes=["CSF1R","FCGR1A","FCGR1B","MAFB","CD68","ITGAM","CD14","MRC1","TREM2","APOE","VSIG4","LILRB1","LILRB2","LILRB3","LILRB4","LILRB5","ITGAX","LAMP3","CCR7","XCR1","IRF8","CD1C","CD163","CD3D"]
```

```python
adata_Full.uns['log1p']["base"] = None
```

```python
sc.pp.highly_variable_genes(adata_Full, flavor="seurat_v3",subset=True, n_top_genes=3000)
sc.tl.pca(adata_Full, svd_solver='arpack')
#sce.pp.harmony_integrate(adata_myeloid, "orig.ident")
sc.pp.neighbors(adata_Full, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata_Full)
sc.tl.rank_genes_groups(adata_Full, 'genevector_updated', method='t-test', use_raw=False)
sc.pl.rank_genes_groups_dotplot(adata_Full,use_raw=False,min_logfoldchange=0.25,title="Myeloids")
sc.pl.umap(adata_Full,color=["genevector_updated"],use_raw=False,title="Myeloids")
```

```python
adata_cancer_epithelial=adata_Full[(adata_Full.obs["genevector_updated"]=="Epithelial")|(adata_Full.obs["genevector_updated"]=="Cancer")]
```

```python
set(adata_cancer_epithelial.obs["genevector_updated"])
```

```python
adata_cancer_epithelial.uns['log1p']["base"] = None
```

```python
sc.pp.highly_variable_genes(adata_cancer_epithelial, flavor="seurat_v3",subset=True, n_top_genes=3000)
sc.tl.pca(adata_cancer_epithelial, svd_solver='arpack')
#sce.pp.harmony_integrate(adata_myeloid, "orig.ident")
sc.pp.neighbors(adata_cancer_epithelial, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata_cancer_epithelial)
sc.tl.rank_genes_groups(adata_cancer_epithelial, 'genevector_updated', method='t-test', use_raw=False)
sc.pl.rank_genes_groups_dotplot(adata_cancer_epithelial,use_raw=False,min_logfoldchange=0.25,title="Myeloids")
sc.pl.umap(adata_cancer_epithelial,color=["genevector_updated"],use_raw=False,title="Myeloids")
```

```python
sc.pl.violin(adata_Full,["CLIC3","JCHAIN","GZMB","FCER1G","CCR7","MZB1","JCHAIN","CD79A","IRF4"],groupby="genevector_updated",use_raw=False,rotation=90) 
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
#look only at 4 subgroups of Monocytes+Macrophages+PDCs+DCs
```

```python
adata_myeloid=adata_Full[(adata_Full.obs["genevector_updated"]=="Dendritic_cells")|(adata_Full.obs["genevector_updated"]=="Macrophages")|(adata_Full.obs["genevector_updated"]=="Monocytes")|(adata_Full.obs["genevector_updated"]=="Plasmacytoid_DCs")]
```

```python
 adata_myeloid.X[0,:].todense().sum()
```

```python
adata_myeloid
```

```python
adata_myeloid.obs.index=adata_myeloid.obs["index"]
```

```python
adata_myeloid.obs["genevector_updated"]
```

```python
fig, axes = plt.subplots(1, 1, figsize=(5,3), gridspec_kw={'wspace':0.2})
celltype_composition(adata_myeloid, axes,"orig.ident","genevector_updated", 
                     list(set(adata_myeloid.obs["orig.ident"])), legend=False)
```

```python
sc.tl.pca(adata_myeloid, svd_solver='arpack')
sce.pp.harmony_integrate(adata_myeloid, "orig.ident")
sc.pp.neighbors(adata_myeloid, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_myeloid)
sc.pl.umap(adata_myeloid,color=["genevector_updated"],use_raw=False)
```

```python
set(adata_myeloid.obs["genevector_updated"])
```

```python
order=['Dendritic_cells', 'Macrophages', 'Monocytes', 'Plasmacytoid_DCs']
```

```python
barplot(adata_myeloid,"genevector_updated","Cell Types",adata_myeloid.uns["genevector_updated_colors"],order=order)
```

```python
sc.pl.violin(adata_myeloid,["CSF1R","MAFB","APOE"],groupby="genevector_updated",use_raw=False,rotation=90) 
```

```python
#Recluster the myeloids
```

```python
#sc.pp.normalize_total(adata_myeloid, target_sum=1e4)
#sc.pp.log1p(adata_myeloid)

adata_myeloid.uns['log1p']["base"] = None
```

```python
#adata_myeloid.X[0,:].todense().sum()  different number just because you're calculating the highly variable genes one more time
```

```python
sc.pp.highly_variable_genes(adata_myeloid, flavor="seurat_v3",subset=True, n_top_genes=3000)
sc.tl.pca(adata_myeloid, svd_solver='arpack')
#sce.pp.harmony_integrate(adata_myeloid, "orig.ident")
sc.pp.neighbors(adata_myeloid, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata_myeloid)
sc.tl.leiden(adata_myeloid,resolution=0.2)
sc.tl.rank_genes_groups(adata_myeloid, 'leiden', method='t-test', use_raw=False)
sc.pl.rank_genes_groups_dotplot(adata_myeloid,use_raw=False,min_logfoldchange=0.25,title="Myeloids")
sc.pl.umap(adata_myeloid,color=["leiden","genevector_updated"],use_raw=False,title="Myeloids")
```

```python
sc.pl.violin(adata_myeloid,color=["SPP1","leiden"],use_raw=False,ncols=3)
```

```python
sc.pl.violin(adata_myeloid,["CD1A","APOE","TREM2","S100AB","VCAM"],groupby="genevector_updated",use_raw=False,rotation=90) 
```

```python
fig, axes = plt.subplots(1, 1, figsize=(5,3), gridspec_kw={'wspace':0.2})
celltype_composition(adata_myeloid, axes,"genevector_updated","leiden", 
                     list(set(adata_myeloid.obs["genevector_updated"])), legend=False)
```

```python
fig, axes = plt.subplots(1, 1, figsize=(5,3), gridspec_kw={'wspace':0.2})
celltype_composition(adata_myeloid, axes,"leiden","genevector_updated", 
                     list(set(adata_myeloid.obs["leiden"])), legend=False)
```

```python
adata_myeloid.obs["leiden"]
```

```python
adata_myeloid.var
```

```python
adata_myeloid
```

```python
adata_myeloid_complete=adata_Full[(adata_Full.obs["genevector_updated"]=="Dendritic_cells")|(adata_Full.obs["genevector_updated"]=="Macrophages")|(adata_Full.obs["genevector_updated"]=="Monocytes")|(adata_Full.obs["genevector_updated"]=="Plasmacytoid_DCs")]
```

```python
adata_myeloid_complete
```

```python
adata_myeloid_complete.obs=adata_myeloid.obs
adata_myeloid_complete.obsm=adata_myeloid.obsm
```

```python
adata_myeloid_complete
```

```python
sc.pl.violin(adata_myeloid_complete,["CD1A","APOE","TREM2","CD1C","LAMP3","S100A8"],groupby="genevector_updated",use_raw=False,rotation=90,ncols=3) 
```

```python
for i in range(0,23,3):
    sc.pl.violin(adata_myeloid_complete,[genes[i],genes[(i+1)],genes[(i+2)]],groupby="genevector_updated",use_raw=False,rotation=90) 
```

```python

```

```python
#Cellphonedb Myeloids
adata_myeloid_Cancer=adata_Full[(adata_Full.obs["genevector_updated"]=="Dendritic_cells")|(adata_Full.obs["genevector_updated"]=="Macrophages")|(adata_Full.obs["genevector_updated"]=="Monocytes")|(adata_Full.obs["genevector_updated"]=="Plasmacytoid_DCs")|(adata_Full.obs["genevector_updated"]=="Cancer")]
adata_myeloid_Cancer.obs["celltype"]=adata_myeloid_Cancer.obs["genevector_updated"]
adata_myeloid_Cancer.write("/work/shah/users/Mobina/Project_Squamous/Squamousdata/adata_myeloid_Cancer_cellphonedb.h5ad")
```

```python
set(adata_myeloid_Cancer.obs["genevector_updated"])
```

```python
#bsub -J "jobname2" -R "rusage[mem=8]" -R "select[type==CentOS7]" -W 02:00 -n 16 -o output -e error sh cellphonedb.sh
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
