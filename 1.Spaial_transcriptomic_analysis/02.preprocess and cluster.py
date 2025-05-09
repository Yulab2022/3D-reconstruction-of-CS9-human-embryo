import squidpy as sq
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt


adata=sc.read_h5ad("cs9_embryo.h5ad")
adata.obsm['spatial'] = np.array(adata.obs[['spatial_X','spatial_Y']])
adata.var["pct_counts_mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, percent_top=None, qc_vars=["pct_counts_mt"], log1p=False, inplace=True)

adata.raw=adata
adata

##filter
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
adata = adata[adata.obs["percent.mito"] < 5].copy()

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata= adata[:,adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sq.gr.spatial_neighbors(adata)

##spatially constrained clustering (SCC)
conn = adata.obsp['connectivities']
conn.data[conn.data > 0] = 1
adj = conn + adata.obsp['spatial_connectivities']
adj.data[adj.data  > 0] = 1


sc.tl.umap(adata)
sc.tl.leiden(adata, adjacency=adj, resolution = 2, key_added = "clusters_1.2_adj" )

meta_data=adata.obs
meta_data.to_csv("meta_data_cs9_embryo.csv",index=True)

###umap for FigS1C
plt.rcParams['figure.figsize']=[30.0,10.0]
sc.pl.umap(adata, color=["clusters_1.2_adj"],add_outline=False,size=30)
plt.savefig("sq_umap_cs9_embryo.pdf")


###spatial plot for FigS1B
plt.rcParams['figure.figsize']=[8.0,64.0]
sc.pl.spatial(adata,spot_size=50, color='clusters_1.2_adj')
plt.savefig("sq_spatial_cs9_embryo.pdf")

adata.write("cs9_embryo_annotation.h5ad")