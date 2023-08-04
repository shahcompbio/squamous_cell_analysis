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
adata = sc.read("../Project_Squamous/Squamousdata/Final_data.h5ad")
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
adata.obs.groupby('orig.ident').plot(kind='scatter', x='n_genes_by_counts', y='pct_counts_mt',subplots=True)
```

```python
adata.obs
```

```python
adata.obs.groupby('orig.ident').plot(kind='scatter', x='total_counts', y='n_genes_by_counts',subplots=True,s=1,cmap= "tab10",c="genevector")
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
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2500,batch_key="orig.ident")
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
"CD3E" in adata.var.index.tolist()
```

```python
#adata.X = adata.X.todense()
adata_original=adata.copy()
```

```python

```

```python
adata
```

```python
dataset = GeneVectorDataset(adata)
```

```python
#import pickle
#mi_scores = pickle.load(open("mi_scores.p","rb"))
```

```python
#dataset.mi_scores=mi_scores
```

```python
# batch_size=1000000000000
#cmps =GeneVector(dataset,
#                   output_file="genes_1200removed.vec",
#                   c=10,  
#                   emb_dimension=100)

```

```python
#import pickle
#pickle.dump(dict(dataset.mi_scores),open("mi_scores.p","wb"))
#mi_scores = pickle.load(open("mi_scores.p","rb"))
#cmps.train(1000, threshold=1e-6)

#for _ in range(100):
  # cmps.train(60, threshold=1e-6) # run for 1000 iterations or loss delta below 1e-6.
```

<!-- #raw -->
cmps.plot()
<!-- #endraw -->

```python
embed = GeneEmbedding("genes_doublet_removed.vec", dataset, vector="average")
```

```python
embed
```

```python
"CD3E" in adata.var.index.tolist()
```

```python
gdata = embed.get_adata(resolution=40)
metagenes = embed.get_metagenes(gdata)
```

```python
gdata
```

```python
cembed = CellEmbedding(dataset, embed)

```

cembed

```python
adata = cembed.get_adata()
```

```python
adata
```

```python
#adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
```

```python
#adata.write("../Project_Squamous/Squamousdata/squmous_data_genevector_embed.h5ad")
```

```python
#adata = sc.read("../Project_Squamous/Squamousdata/squmous_data_genevector_embed.h5ad")
```

```python
"KLRD1" in adata.var.index
```

```python
adata_original=adata.copy()
```

```python
embed.compute_similarities("LYZ")[:20]
```

```python
#sc.pl.pca(adata,color=["FCER1G","LYZ","AIF1","APOE","G0S2","TKT","KIT","CLV"],use_raw=False,ncols=2)
```

```python
"CLU" in adata.var.index.tolist()
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
adata = phenotype_probability(cembed,adata,markers)
```

```python
adata.obs
```

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

```python
sc.pl.umap(adata,color="orig.ident")
```

```python
sc.pl.umap(adata,color=["genevector"],use_raw=False)
```

```python
ct = []
cells = []
for x in set(adata.obs["genevector"]):
    cells.append(len(adata[adata.obs["genevector"]==x].obs.index))
    ct.append(x)
df = pandas.DataFrame.from_dict({"Cells":cells,"Cell Type":ct})
import seaborn as sns
sns.barplot(data=df,x="Cell Type",y="Cells")
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
fig, axes = plt.subplots(1, 1, figsize=(5,3), gridspec_kw={'wspace':0.2})
celltype_composition(adata, axes, "orig.ident","genevector", 
                     list(set(adata.obs["orig.ident"])), legend=False)
```

```python
sc.tl.rank_genes_groups(adata, 'genevector', method='t-test', use_raw=False)
```

```python
sc.pl.rank_genes_groups_dotplot(adata,use_raw=False,swap_axes=True,min_logfoldchange=1.)
```

```python
vmarkers = dict()
vgenes = []
for ph,genes in markers.items():
    vgenes += list(set(adata.var.index).intersection(set(genes)))
print(genes)
```

```python
vgenes
```

```python
sc.pl.dotplot(adata,np.unique(vgenes),groupby="genevector", use_raw=False)
```

```python
adata_Full = adata_QC.copy()
sc.pp.normalize_total(adata_Full, target_sum=1e4)
sc.pp.log1p(adata_Full)
sc.pp.calculate_qc_metrics(adata_Full,inplace=True)
```

```python
adata_Full.shape
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
adata_Full.shape
```

```python
adata_Full.obs=adata.obs
```

```python
adata_Full.obsm=adata.obsm
```

```python
adata_Full
```

```python
adata_Full.obs
```

```python
sc.pl.umap(adata_Full,color=["genevector","NCAM1","CD3D","KLRG1"],use_raw=False,ncols=1,vmax=1)
```

```python
sc.pl.umap(adata_Full,color=["CSF1R","LILRB4","TYROBP","MAFB","APOE"],use_raw=False,ncols=2)
```

```python
sc.pl.violin(adata_Full,"PTPRC",groupby="genevector",use_raw=False)
```

```python
markers.keys()
```

```python
ph='T Cell'
Data = adata_Full[adata_Full.obs["genevector"]==ph]
#sc.pp.highly_variable_genes(Data, flavor="seurat_v3",subset=True, n_top_genes=2500)
sc.tl.pca(Data, svd_solver='arpack')
```

```python
sc.pp.neighbors(Data, n_neighbors=10, n_pcs=40)
sc.tl.umap(Data)
sc.tl.leiden(Data,resolution=0.05)
sc.pl.umap(Data,color=["leiden","orig.ident"],use_raw=False,title=ph)
```

```python
sc.pl.umap(adata_Full,color=["EPCAM","PTPRC","genevector"],use_raw=False)
```

```python
sc.tl.rank_genes_groups(Data, 'leiden', method='t-test', use_raw=False)
sc.pl.rank_genes_groups_dotplot(Data,use_raw=False,min_logfoldchange=1.,title=ph)
```

```python
Data
```

```python
fig, axes = plt.subplots(1, 1, figsize=(5,3), gridspec_kw={'wspace':0.2})
celltype_composition(Data, axes, "orig.ident","leiden", 
                     list(set(Data.obs["orig.ident"])), legend=False)
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
#!pip install harmonypy
```

```python
Data_celltype=[]
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
  Data_celltype.append(Data)
# print(ph)
```

```python
Data_celltype_2=[]
for ph,genes in markers.items():
  Data = adata_Full[adata_Full.obs["genevector"]==ph]
  #sc.pp.highly_variable_genes(Data, flavor="seurat_v3",subset=True, n_top_genes=2500,batch_key="orig.ident")
  sc.tl.pca(Data, svd_solver='arpack')
  sc.pp.neighbors(Data, n_neighbors=10, n_pcs=40)
  sc.tl.umap(Data)
  sc.tl.leiden(Data,resolution=0.05)
  sc.tl.rank_genes_groups(Data, 'leiden', method='t-test', use_raw=False)
  sc.pl.rank_genes_groups_dotplot(Data,use_raw=False,min_logfoldchange=1.,title=ph)
  sc.pl.umap(Data,color=["leiden","orig.ident"],use_raw=False,title=ph)
  fig, axes = plt.subplots(1, 1, figsize=(5,3), gridspec_kw={'wspace':0.2})
  celltype_composition(Data, axes, "orig.ident","leiden", 
                     list(set(Data.obs["orig.ident"])), legend=False)
  #sc.pl.umap(Data,color=["orig.ident"],use_raw=False,title=ph)  
  Data_celltype_2.append(Data)
# print(ph)
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
Data_celltype.append(Data)
```

```python
sc.pl.umap(adata,color=["genevector"],use_raw=False)
```

```python
Data_celltype[0].obs
```

```python
sns.violinplot(x="orig.ident",y="pct_counts_mt",data=adata_Full.obs)
plt.xticks(rotation=90)
plt.show()
```

```python
sns.violinplot(x="orig.ident",y="pct_counts_mt",data=Data_celltype[0].obs)
plt.xticks(rotation=90)
plt.show()

```

```python
for i in range(6):
    print(set(Data_celltype[i].obs["leiden"]))
```

```python
markers.keys()
```

```python
di_Tcell= {"0": "", "1": "", "2": "", "3":""}
di_Bcell= {"0": "" ,"1": "", "2": ""}
di_Plasma = {"0": "", "1": ""}
di_Myeloid= {"0":"","1": "", "2": "", "3":"","4":"","5":""}
di_Fibroblast= {"0":"","1": "", "2": "", "3":""}
```

```python
Data_celltype[0].obs.replace({"leiden": di_Tcell})
```

```python
Data_celltype_concat= pd.concat(Data_celltype, ignore_index=True)
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

```python

```

```python
#####Cellphonedb
adata_Full.obs["celltype"]=adata_Full.obs["genevector"]
```

```python
adata_Full.obs
```

```python
df_meta = pd.DataFrame(data={'Cell':list(adata.obs.index),
                             'cell_type':[ i for i in adata.obs['celltype']]
                            })
```

```python
df_meta.set_index('Cell', inplace=True)
```

```python
df_meta
```

```python
df_meta.to_csv('../Project_Squamous/Squamousdata/squamous_metadata.tsv', sep = '\t')
```

```python
adata_Full.__dict__['_raw'].__dict__['_var'] = adata_Full.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
```

```python
adata.write("../Project_Squamous/Squamousdata/squamous_cellphonedb.h5ad")
```

```python
adata.var
```

```python

```
