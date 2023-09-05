# scgeom.pp.cell_specific_gene_network

### scgeom.pp.cell_specific_gene_network(adata, method='csn', csn_alpha=0.01, csn_boxsize=0.1, csn_weighted=True, network_name='cell_specific_gene_network', return_network=False)

Compute a cell-specific gene network for each cell.

* **Parameters:**
  * **adata** (`AnnData`) – The AnnData object of gene expression data.
  * **method** (`str`) – The method to used. ‘csn’ for cell-specific network [Dai2019] or ‘loccsn’ for local cell-specific network [Wang2021].
  * **csn_alpha** (`float`) – When method is ‘csn’, the alpha value (significance for keeping an edge).
  * **csn_boxsize** (`float`) – When method is ‘csn’, the size of box to define neiborhood of gene expression.
  * **csn_weighted** (`bool`) – When method is ‘csn’, whether use weighted output or binary output.
  * **network_name** (`str`) – The name of the dictionary of networks to be stored in `adata.uns[network_name]`. Use the `obs_names` of cells to access the network for each cell.
  * **return_network** (`bool`) – Whether to return a copy of the dictionary of networks.
* **Returns:**
  **gene_networks** – A dictionary of cell-specific networks with cell names in `adata.obs_names` as dictionary keyssss.
* **Return type:**
  dict

### References

[Dai2019] Dai, H., Li, L., Zeng, T., & Chen, L. (2019). Cell-specific network constructed by single-cell RNA sequencing data. Nucleic acids research, 47(11), e62-e62.

[Wang2021] Wang, X., Choi, D., & Roeder, K. (2021). Constructing local cell-specific networks from single-cell data. Proceedings of the National Academy of Sciences, 118(51), e2113178118.
