# scgeom.tl.graph_curvature

### scgeom.tl.graph_curvature(adata, graph_name='scgeom-graph', curvature_method='orc', orc_alpha=0.5, orc_method='OTD', node_curvature_name='scgeom-node_curvature', edge_curvature_name='scgeom-edge_curvature')

Compute graph curvature for each cell. We use the package GraphRicciCurvature [Ni2019] to compute the curvatures.

* **Parameters:**
  * **adata** (`AnnData`) – The `anndata` object of gene expression data.
  * **graph_name** (`str`) – The graph stored in `adata.obsp[graph_name]` will be used.
  * **curvature_method** (`str`) – ‘orc’ for Ollivier-Ricci curvature and ‘frc’ for Forman-Ricci curvature.
  * **orc_alpha** (`float`) – The alpha paramter for computing Ollivier-Ricci curvature. The weight of mass to be put on the center node when defining a mass distribution of a node’s neighborhood.
  * **orc_method** (`str`) – The method for computing distance between mass distributions.
    ‘OTD’ for optimal transport distance (earth mover’s distance),
    ‘Sinkhorn’ for entropy regularized approximation of optimal transport distance,
    ‘ATD’ for average transportation distance.
  * **node_curvature_name** (`str`) – The computed node curvature will be stored in `adata.obs[node_curvature_name]`.
  * **edge_curvature_name** (`str`) – The computed edge curvature will be stored in `adata.obsp[edge_curvature_name]`.

### References

[Ni2019] Ni, C. C., Lin, Y. Y., Luo, F., & Gao, J. (2019). Community detection on networks with Ricci flow. Scientific reports, 9(1), 9984.
