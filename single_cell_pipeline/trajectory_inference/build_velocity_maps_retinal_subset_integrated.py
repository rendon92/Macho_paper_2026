import scanpy as sc
import pandas as pd
from scipy.io import mmread
import anndata as ad
import scvelo as scv
import matplotlib.pyplot as plt



# Load integrated counts
X = mmread("/data/scMultiome/paper/files/github_code_paper/single_cell_pipeline/trajectory_inference/retina_integrated.mtx").T.tocsr()  # transpose if needed

# Load genes
genes = pd.read_csv("/data/scMultiome/paper/files/github_code_paper/single_cell_pipeline/trajectory_inference/retina_integrated_genes.csv")
adata = ad.AnnData(X=X)
adata.var_names = genes.iloc[:,0]

# Load metadata
meta = pd.read_csv("/data/scMultiome/paper/files/github_code_paper/single_cell_pipeline/trajectory_inference/retina_integrated_meta.csv", index_col=0)
adata.obs = meta

# Load PCA / UMAP if available
pca = pd.read_csv("/data/scMultiome/paper/files/github_code_paper/single_cell_pipeline/trajectory_inference/retina_integrated_pca.csv", index_col=0)
umap = pd.read_csv("/data/scMultiome/paper/files/github_code_paper/single_cell_pipeline/trajectory_inference/retina_integrated_umap.csv", index_col=0)
adata.obsm["X_pca"] = pca.values
adata.obsm["X_umap"] = umap.values

import scanpy as sc

# Plot using default settings

ldata1 = scv.read('/data/scRNAseq/experiment/yard/run_cellranger_count/run_count_JR_SCR_15/velocyto/run_count_JR_SCR_15.loom', cache=False)
ldata2 = scv.read('/data/scRNAseq/experiment/yard/run_cellranger_count/run_count_JR_SCR_18/velocyto/run_count_JR_SCR_18.loom', cache=False)
ldata3 = scv.read('/data/scRNAseq/experiment/yard/run_cellranger_count/run_count_JR_SCR_23/velocyto/run_count_JR_SCR_23.loom', cache=False)

ldata1.var_names_make_unique()
ldata2.var_names_make_unique()
ldata3.var_names_make_unique()

#adata 15, 18 and 23 merged
ldata = ldata1.concatenate([ldata2, ldata3])

adata = scv.utils.merge(adata, ldata)

adata.var_names_make_unique()

scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)





# Color dictionary same as R code
color_dict = {
 
    "Neural retina": "#db8e00ff",
    "Proliferative region": "#569cc9ff",
}

# Transform dictionary to list (same ordering)
metadata_column = "cluster"  # change to correct name
categories = adata.obs[metadata_column].unique()  # if categorical
palette = [color_dict[cat] for cat in categories]

# Plot
sc.pl.umap(adata, color=metadata_column, palette=palette, legend_loc='on data')


scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)


scv.pl.velocity_embedding_stream(
	adata, basis="umap", legend_fontsize=12, title="", smooth=0.8, min_mass=4, color='cluster', show=False)
plt.savefig("velocity_umap_retina.svg", format="svg", bbox_inches="tight")
plt.close()


scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80)


from collections import Counter


num_genes = 1000

# Select top genes according to several parameteres
top_genes_likelihood = adata.var['fit_likelihood'].sort_values(ascending=False).index[:num_genes]
top_genes_alpha = adata.var['fit_alpha'].sort_values(ascending=False).index[:num_genes]
top_genes_beta = adata.var['fit_beta'].sort_values(ascending=False).index[:num_genes]
top_genes_gamma = adata.var['fit_gamma'].sort_values(ascending=False).index[:num_genes]
top_genes_variance = adata.var['fit_variance'].sort_values(ascending=False).index[:num_genes]


all_top_genes = list(top_genes_likelihood) + list(top_genes_alpha) + \
                list(top_genes_beta) + list(top_genes_gamma) + list(top_genes_variance)

gene_counts = Counter(all_top_genes)
driver_genes = [gene for gene, count in gene_counts.items() if count >= 3]


with open("/data/scMultiome/paper/figure_all_data/figure1/retinal_markers.tsv", "r") as file:
    gene_list = [line.strip() for line in file.readlines()]


common_genes = set(gene_list).intersection(driver_genes)

cg = scv.pl.heatmap(
    adata,
    var_names=common_genes,
    sortby='latent_time',
    col_color='cluster',
    n_convolve=100,
    row_cluster=False,
    font_scale=0.4,
    show=False
)

scv.pl.heatmap(
    adata,
    var_names=common_genes,
    sortby='latent_time',
    col_color='cluster',
    n_convolve=100,
    row_cluster=False,
    font_scale=0.4,
    show=False
)
plt.savefig("velocity_genes_retina_heatmap.png", dpi=300, bbox_inches="tight")
plt.close()



ordered_common_genes = cg.data2d.index.tolist()


import pandas as pd
pd.DataFrame(ordered_common_genes, columns=['Gene']).to_csv('/data/scMultiome/paper/figure_all_data/figure1/driver_DEGs_ordered.csv', index=False)

