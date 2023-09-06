# scgeom.tl.cell_network_topology

### scgeom.tl.cell_network_topology(adata, network_name='connectivities', embedding_name='X_pca', method='lwph', rlph_method='self', rlph_base='all', metric='euclidean', nb_method='knn', nb_distance=20, rl_distance=10, nb_knn=20, rl_knn=10, filtration_value='euclidean', max_filtration=inf, max_dim=1, n_cores=20, uns_name='cell_network_diagrams', return_diagrams=False)

Compute topology for each cell on the cell network.

* **Parameters:**
  * **adata** (`AnnData`) – The `adata` object of gene expression data.
  * **network_name** (`str`) – The network of cells, e.g. a knn network to be used in `adata.obsp[network_name]`.
  * **embedding_name** (`str`) – The low-dimensional embedding of cells to be used to compute filtration values, stored in `adata.obsm[embedding_name]`.
  * **method** (`str`) – The persistent homology method to use.
    ‘lwph’ for local weighted persistent homology which computes Vietoris-Rips complex based PH on a local neighborhood of each cell.
    ‘rlph’ for relative local persistent homology which computes the relative persistent homology of the cell, i.e., the global topology relative to a cell or its local neighborhood.
  * **rlph_method** (`str`) – When `method` is set to ‘rlph’, choose ‘self’ to compute persistent homology relative only to the cell,
    ‘knn’ for relative to a neighborhood defined by k-nearest neighbors, or ‘distance’ for relative to a neighborhood defined by distance cutoff.
  * **rlph_base** (`str`) – What points to include when computing the global topology. ‘all’ for including all points,
    ‘nb_knn’ or ‘nb_distance’ for a relatively large neighborhood defined by k-nearest neighbors or distance cutoff, respectively.
  * **metric** (`str`) – The metric used to define base point set for computing global topology when `rlph_base` is set to ‘nb_knn’ or ‘nb_distance’.
  * **nb_method** (`str`) – The method to define neighborhood of cells when `method` is set to ‘lwph’. ‘knn’ for k-nearest neighbors or ‘distance’ for distance cutoff.
  * **nb_distance** (`float`) – If `nb_method` is set to ‘distance’ in ‘lwph’ or `rlph_base` is set to ‘nb_distance’ in ‘rlph’, the distance cutoff to use.
  * **rl_distance** (`float`) – If `method` is set to ‘rlph’ and `rlph_method` is set to ‘distance’, the distance cutoff to use to define the local cell neighborhood.
  * **nb_knn** (`int`) – If `nb_method` is set to ‘knn’ in ‘lwph’ or `rlph_base` is set to ‘nb_knn’ in ‘rlph’, the k value to use.
  * **rl_knn** (`int`) – If `method` is set to ‘rlph’ and `rlph_method` is set to ‘knn’, the k value to use to define the local cell neighborhood.
  * **filtration_value** (`str`) – When ‘lwph’ is used, the filtration value to use. Currently, only Euclidean distance in the low-dimensional embedding is used.
  * **max_filtration** (`float`) – The maximum filtration value. If `None`, it is set to the maximum node weight.
  * **max_dim** (`int`) – The maximum dimension for computing persistent homology.
  * **n_cores** (`int`) – Number of cores to use for computing persistent homology for each cell.
  * **return_diagrams** (`bool`) – Whether to return the computed persistence diagrams.
  * **uns_name** (`str`) – The persistence diagrams are stored in `adata.uns[uns_name]`.
* **Returns:**
  **diagrams** – A dictionary of diagrams to be stored in `adata.uns[uns_name]`.
  For ‘lwph’, the H1 diagram for the cell named ‘A’ can be accessed through `diagrams['A']['dgm_h1]`
  and the corresponding max_filtration value is in `diagrams['A']['max_filtration']`.

  For ‘rlph’, the H1 relative PH results for base point set relative to local point set is in `diagrams['A']['dgm_h1']`.
  Additionally, the H1 regular PH results for the base point set is in `diagrams['A']['dgm_h1_base']`.
* **Return type:**
  dict
