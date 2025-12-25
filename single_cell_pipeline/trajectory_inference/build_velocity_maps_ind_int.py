import scanpy as sc
import pandas as pd
from scipy.io import mmread
import anndata as ad
import scvelo as scv 

import matplotlib.pyplot as plt



# Load counts
dr15 = mmread("dr15.mtx").T.tocsr()  # transpose if needed
dr18 = mmread("dr18.mtx").T.tocsr()  # transpose if needed
dr23 = mmread("dr23.mtx").T.tocsr()  # transpose if needed

# Load genes
genes15 = pd.read_csv("dr15_genes.csv")
genes18 = pd.read_csv("dr18_genes.csv")
genes23 = pd.read_csv("dr23_genes.csv")

adata15 = ad.AnnData(X=dr15)
adata15.var_names = genes15.iloc[:,0]
adata18 = ad.AnnData(X=dr18)
adata18.var_names = genes18.iloc[:,0]
adata23 = ad.AnnData(X=dr23)
adata23.var_names = genes23.iloc[:,0]



# Load metadata
meta15 = pd.read_csv("dr15_meta.csv", index_col=0)
adata15.obs = meta15
meta18 = pd.read_csv("dr18_meta.csv", index_col=0)
adata18.obs = meta18
meta23 = pd.read_csv("dr23_meta.csv", index_col=0)
adata23.obs = meta23

# Load PCA / UMAP if available
pca15 = pd.read_csv("dr15_pca.csv", index_col=0)
umap15 = pd.read_csv("dr15_umap.csv", index_col=0)
adata15.obsm["X_pca"] = pca15.values
adata15.obsm["X_umap"] = umap15.values
pca18 = pd.read_csv("dr18_pca.csv", index_col=0)
umap18 = pd.read_csv("dr18_umap.csv", index_col=0)
adata18.obsm["X_pca"] = pca18.values
adata18.obsm["X_umap"] = umap18.values
pca23 = pd.read_csv("dr23_pca.csv", index_col=0)
umap23 = pd.read_csv("dr23_umap.csv", index_col=0)
adata23.obsm["X_pca"] = pca23.values
adata23.obsm["X_umap"] = umap23.values

# Plot using default settings


ldata1 = scv.read('run_count_JR_SCR_15.loom', cache=False)
ldata2 = scv.read('run_count_JR_SCR_18.loom', cache=False)
ldata3 = scv.read('run_count_JR_SCR_23.loom', cache=False)


#merging h5ad data with loom files

# make variable names unique

ldata1.var_names_make_unique()
ldata2.var_names_make_unique()
ldata3.var_names_make_unique()

#adata 15, 18 and 23 merged
#ldata = ldata1.concatenate([ldata2, ldata3])


adata15 = scv.utils.merge(adata15, ldata1)
adata18 = scv.utils.merge(adata18, ldata2)
adata23 = scv.utils.merge(adata23, ldata3)

adata15.var_names_make_unique()
adata18.var_names_make_unique()
adata23.var_names_make_unique()

scv.pp.filter_and_normalize(adata15)
scv.pp.moments(adata15)
scv.pp.filter_and_normalize(adata18)
scv.pp.moments(adata18)
scv.pp.filter_and_normalize(adata23)
scv.pp.moments(adata23)



metadata_column = "cluster"  

color_dict = {
    "Diencephalon": "#ffd7a8a1",
    "Forebrain": "#97a900ff",
    "Lens placode": "#2fb600ff",
    "Neural retina": "#569cc9ff",
    "Neural tube": "#00a6ffff",
    "Optic Stalk": "#00c0b7ff",
    "Periderm": "#c8ffc3ff",
    "Proliferative region": "#db8e00ff",
    "RPE": "#ad16c699",
    "Telencephalon": "#fba6a7ff",
    "Trigeminal placode": "#ff63b6ff",
    "Neural crest": "#6fa690ff",
    "primary neuron": "#b19641ff",
    "Midbrain": "#b19641ff",
    "Hindbrain": "#b19641ff"
}




color_dict_15 = {
    "Periderm": "#c8ffc3ff",
    "Proliferative region": "#db8e00ff",
    "Neural retina": "#569cc9ff",
    "Diencephalon": "#ffd7a8a1",
    "Telencephalon": "#fba6a7ff",
    "RPE": "#ad16c699",
    "Optic Stalk": "#00c0b7ff",
    "Neural tube": "#00a6ffff",
    "Lens placode": "#2fb600ff",
    "Forebrain": "#97a900ff",
    "Trigeminal placode": "#ff63b6ff",
    
}    



# Transform cluster column into categorical (same ordering as in the dictionary)
cluster_order = list(color_dict_15.keys())

adata15.obs[metadata_column] = pd.Categorical(
    adata15.obs[metadata_column],
    categories=cluster_order,
    ordered=True
)
# Create colo palette (same ordering)
palette = [color_dict_15[cat] for cat in cluster_order if cat in adata15.obs[metadata_column].cat.categories]

# Plot 15h
sc.pl.umap(
    adata15,
    color=metadata_column,
    palette=palette,
    legend_loc='on data',
    size=20
)

color_dict_18 = {
    "Neural retina": "#569cc9ff",
    "RPE": "#ad16c699",
    "Trigeminal placode": "#ff63b6ff",
    "Lens placode": "#2fb600ff",
    "Optic Stalk": "#00c0b7ff",
    "Telencephalon": "#fba6a7ff",
    "Proliferative region": "#db8e00ff",
    "Primary neuron": "#b19641ff",
    "Diencephalon": "#ffd7a8a1",
    "Hindbrain": "#b19641ff",
    "Periderm": "#c8ffc3ff",
    "Midbrain": "#b19641ff",
    

}


cluster_order = list(color_dict_18.keys())

adata18.obs[metadata_column] = pd.Categorical(
    adata18.obs[metadata_column],
    categories=cluster_order,
    ordered=True
)

palette = [color_dict_18[cat] for cat in cluster_order if cat in adata18.obs[metadata_column].cat.categories]

# Plot 18h
sc.pl.umap(
    adata18,
    color=metadata_column,
    palette=palette,
    legend_loc='on data',
    size=20 
)

color_dict_23 = {
    "Diencephalon": "#ffd7a8a1",
    "Neural retina": "#569cc9ff",
    "Proliferative region": "#db8e00ff",
    "RPE": "#ad16c699",
    "Neural crest": "#6fa690ff",
    "Optic Stalk": "#00c0b7ff",
    "Lens placode": "#2fb600ff",
    "Telencephalon": "#fba6a7ff"
  
}


cluster_order = list(color_dict_23.keys())
adata23.obs[metadata_column] = pd.Categorical(
    adata23.obs[metadata_column],
    categories=cluster_order,
    ordered=True
)

palette = [color_dict_23[cat] for cat in cluster_order if cat in adata23.obs[metadata_column].cat.categories]

# Plot 23h
sc.pl.umap(
    adata23,
    color=metadata_column,
    palette=palette,
    legend_loc='on data',
    size=20  
)


## RNA velocity dynamical mode

scv.tl.recover_dynamics(adata15, n_jobs=8)
scv.tl.recover_dynamics(adata18, n_jobs=8)
scv.tl.recover_dynamics(adata23, n_jobs=8)
#scv.tl.velocity(adata, mode="dynamical", var_names=adata.var_names)

scv.tl.velocity(adata15, mode="dynamical")
scv.tl.velocity(adata18, mode="dynamical")
scv.tl.velocity(adata23, mode="dynamical")

scv.tl.velocity_graph(adata15)
scv.tl.velocity_graph(adata18)
scv.tl.velocity_graph(adata23)


scv.pl.velocity_embedding_stream(
    adata15, basis="umap", legend_fontsize=12, title="", smooth=0.8, min_mass=4, color='cluster', show=False)
plt.savefig("velocity_umap_15.svg", format="svg", bbox_inches="tight")
plt.close()

scv.pl.velocity_embedding_stream(
    adata18, basis="umap", legend_fontsize=12, title="", smooth=0.8, min_mass=4, color='cluster', show=False)
plt.savefig("velocity_umap_18.svg", format="svg", bbox_inches="tight")
plt.close()

scv.pl.velocity_embedding_stream(
    adata23, basis="umap", legend_fontsize=12, title="", smooth=0.8, min_mass=4, color='cluster', show=False)
plt.savefig("velocity_umap_23.svg", format="svg", bbox_inches="tight")
plt.close()
