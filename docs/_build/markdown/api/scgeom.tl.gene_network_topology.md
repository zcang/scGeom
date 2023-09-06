# scgeom.tl.gene_network_topology

### scgeom.tl.gene_network_topology(adata, network_name='gene_network', filtration='vbcl', max_filtration=None, max_dim=1, n_cores=20, vbcl_edge_weight='min', return_diagrams=False, uns_name='diagrams')

Compute gene network topology for node weighted gene networks.

* **Parameters:**
  * **adata** (`AnnData`) – The `anndata` object of gene expression data.
  * **network_name** (`str`) – The network structure to be used should be available in `adata.uns[network_name]`. Use `pp.get_gene_network` to obtain a knowledge-based gene network to be shared by each cell.
  * **filtration** (`str`) – The filtration method for computing persistent homology. ‘vbcl’ for vertex-based clique complex filtration. Currently, only ‘vbcl’ is implemented.
  * **max_filtration** (`Optional`[`float`]) – The maximum filtration value. If `None`, it is set to the maximum node weight.
  * **max_dim** (`int`) – The maximum dimension for computing persistent homology.
  * **n_cores** (`int`) – Number of cores to use for computing persistent homology for each cell.
  * **vbcl_edge_weight** (`str`) – How filtration value for edges are assigned.
    ‘min’ will assign the smaller node weight to the edge.
    ‘mean’ will use the average of two node weights for each edge.
  * **return_diagrams** (`bool`) – Whether to return the computed persistence diagrams.
  * **uns_name** (`str`) – The persistence diagrams are stored in `adata.uns[uns_name]`.
* **Returns:**
  **diagrams** – A dictionary of diagrams to be stored in `adata.uns[uns_name]`.
  For example, the H1 diagram for the cell named ‘A’ can be accessed through `diagrams['A']['dgm_h1]`
  and the corresponding max_filtration value is in `diagrams['A']['max_filtration']`.
* **Return type:**
  dict
