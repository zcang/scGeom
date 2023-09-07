# scGeom
Geometric and topological analysis of cell states in single-cell data

## Installation
cd to the code directory and run `pip install .`

## Basic usage

_Import packages_
```
import scgeom as sg
import scanpy as sc
from scipy.spatial import distance_matrix
```
_Curvature and topology of cell network_ \
Let `adata` be an anndata object of single-cell data.
```
# Preprocessing
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.pp.pca(adata)
D = distance_matrix(adata.obsm['X_pca'], adata.obsm['X_pca'])
# Computing curvature
sg.tl.neighbor_graph(adata, D=D, knn_k=30, weighted=False, graph_name='curvature_graph')
sg.tl.graph_curvature(adata, curvature_method='orc', orc_alpha=0.5, graph_name='curvature_graph')
# Computing topology with local weighted persistent homology
sc.pp.neighbors(adata)
sg.tl.cell_network_topology(adata, network_name='connectivities', embedding_name='X_pca', nb_method='knn',
    nb_knn=30, max_filtration=20, max_dim=1, method='lwph', uns_name='diagrams_cell_network_lwph)
# or relative local persistent homology
sg.tl.cell_network_topology(adata, network_name = 'connectivities', embedding_name = 'X_pca', rlph_method = 'knn',
    nb_knn=100, rl_knn=20, max_filtration = 20, max_dim=1, method='rlph', rlph_base='nb_knn', uns_name='diagrams_cell_network_rlph')
```

_Topology of gene networks_ \
Let `adata` be an anndata object of original data of a human dataset.
```
# Preprocessing
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sg.pp.get_gene_network(adata, species='human', database='scent_17')
# Computing vertex-based clique complex persistent homology
sg.tl.gene_network_topology(adata, max_filtration=None, filtration='vbcl', vbcl_edge_weight='min', 
    uns_name='diagrams_vbcl',n_cores=20)
# Computing edge-weighted Rips complex persistent homology with cell-specific networks
sc.pp.highly_variable_genes(adata)
adata = adata[:,adata.var.highly_variable]
sg.tl.cell_specific_gene_network_topology(adata, compute_csn=True, csn_weighted=True, uns_name='diagrams_ewvr')
# Generate features from the computed diagrams using, for example, persistence entropy
sg.tl.generate_topology_feature(adata, method='pe', diagram_name='diagrams_ewvr', dims=[0,1], feature_name='X_topo_ewvr_pe')
```

# Documentation and examples
Detailed api documentation can be found [here](docs/_build/markdown/api/index.md).
Several examples can be found in the `example` folder.

# Reference


