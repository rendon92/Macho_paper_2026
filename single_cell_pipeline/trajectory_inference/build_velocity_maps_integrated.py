import scanpy as sc
import pandas as pd
from scipy.io import mmread
import anndata as ad
import scvelo as scv
import matplotlib.pyplot as plt

# Load integrated counts
X = mmread("/data/scMultiome/paper/files/github_code_paper/single_cell_pipeline/trajectory_inference/integrated.mtx").T.tocsr()  # transpose if needed

# Load genes
genes = pd.read_csv("/data/scMultiome/paper/files/github_code_paper/single_cell_pipeline/trajectory_inference/integrated_genes.csv")
adata = ad.AnnData(X=X)
adata.var_names = genes.iloc[:,0]

# Load metadata
meta = pd.read_csv("/data/scMultiome/paper/files/github_code_paper/single_cell_pipeline/trajectory_inference/meta.csv", index_col=0)
adata.obs = meta

# Load PCA / UMAP if available
pca = pd.read_csv("/data/scMultiome/paper/files/github_code_paper/single_cell_pipeline/trajectory_inference/pca.csv", index_col=0)
umap = pd.read_csv("/data/scMultiome/paper/files/github_code_paper/single_cell_pipeline/trajectory_inference/umap.csv", index_col=0)
adata.obsm["X_pca"] = pca.values
adata.obsm["X_umap"] = umap.values

import scanpy as sc

# Plot using default settings


ldata1 = scv.read('/data/scRNAseq/experiment/yard/run_cellranger_count/run_count_JR_SCR_15/velocyto/run_count_JR_SCR_15.loom', cache=False)
ldata2 = scv.read('/data/scRNAseq/experiment/yard/run_cellranger_count/run_count_JR_SCR_18/velocyto/run_count_JR_SCR_18.loom', cache=False)
ldata3 = scv.read('/data/scRNAseq/experiment/yard/run_cellranger_count/run_count_JR_SCR_23/velocyto/run_count_JR_SCR_23.loom', cache=False)


#merging h5ad data with loom files

# make variable names unique

ldata1.var_names_make_unique()
ldata2.var_names_make_unique()
ldata3.var_names_make_unique()

#adata 15, 18 and 23 merged
ldata = ldata1.concatenate([ldata2, ldata3])

adata = scv.utils.merge(adata, ldata)

adata.var_names_make_unique()

scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)


# Color dictionary
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

# Convert dictionary to list (same ordering)
metadata_column = "cluster"  

cluster_order = list(color_dict.keys())

adata.obs[metadata_column] = pd.Categorical(
	adata.obs[metadata_column],
	categories=cluster_order,
	ordered=True
)


palette = [color_dict[cat] for cat in cluster_order if cat in adata.obs[metadata_column].cat.categories]

# Plot


sc.pl.umap(
	adata,
	color=metadata_column,
	palette=palette,
	legend_loc='on data',
	size=20  
)

#### velocity (dinamical mode)
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)


scv.pl.velocity_embedding_stream(
    adata, basis="umap", legend_fontsize=12, title="", smooth=0.8, min_mass=4, color='cluster', show=False)
plt.savefig("velocity_umap_integrated.svg", format="svg", bbox_inches="tight")
plt.close()


scv.pl.velocity_embedding_stream(
    adata, basis="umap", legend_fontsize=12, title="", smooth=0.8, min_mass=4, color='cluster')

