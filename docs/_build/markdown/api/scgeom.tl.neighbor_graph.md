# scgeom.tl.neighbor_graph

### scgeom.tl.neighbor_graph(adata, X=None, D=None, graph_name='scgeom-graph', graph_method='knn', knn_k=10, weighted=True, D_type='distance')

Construct a cell graph for curvature computation.

* **Parameters:**
  * **adata** (`AnnData`) – The `anndata` object of gene expression data.
  * **X** (`Optional`[`array`]) – A low dimensional embedding of data of shape `n_cell x n_dim` where Euclidean distance will be used to quantify similarity.
    If `X` is given, `D` will be ignored.
  * **D** (`Optional`[`array`]) – A `n_cell x n_cell` distance or similarity matrix. If `X` is given, `D` will be ignored.
  * **graph_name** (`str`) – The name to store the cell graph in `adata.obsp[graph_name]`.
  * **graph_method** (`str`) – Currently, only ‘knn’ (k-nearest neighbor graph) is implemented.
  * **knn_k** (`int`) – The number of neighbors when buidling the knn graph.
  * **weighted** (`bool`) – Whether to include weights in the sparse matrix representing the graph.
  * **D_type** (`str`) – ‘distance’ is `D` is a distance/dissimilarity matrix or ‘similarity’ is `D` represents similarity among points.
* **Returns:**
  **A** – A sparse matrix representing the cell graph stored in `adata.obsp[graph_name]`.
* **Return type:**
  sp.coo_matrix
