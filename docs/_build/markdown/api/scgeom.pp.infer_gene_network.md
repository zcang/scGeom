# scgeom.pp.infer_gene_network

### scgeom.pp.infer_gene_network(adata, method='spearman', threshold=0.2, network_name='gene_network', return_network=False)

Simple computation of gene networks with Spearman’s correlation coefficient.

* **Parameters:**
  * **adata** (`AnnData`) – The data matrix of shape `n_obs` × `n_var`.
  * **method** (`str`) – Method to use, currently, only Spearman’s r is implemented.
  * **threshold** (`float`) – The threshold of absolute value of correlation coefficients to keep an edge.
  * **network_name** (`str`) – The name of the network to be stored in `adata.uns[network_name]`.
  * **return_network** (`bool`) – Whether to return the network.
* **Returns:**
  **gene_network** – A pandas DataFrame of the LR pairs with the three columns representing the ligand, receptor, and the signaling pathway name, respectively.
* **Return type:**
  dict
