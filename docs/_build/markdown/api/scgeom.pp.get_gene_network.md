# scgeom.pp.get_gene_network

### scgeom.pp.get_gene_network(adata, database='scent_17', species='human')

Get the gene network from knowledge database from SCENT [Banerji2013].

* **Parameters:**
  * **adata** (`AnnData`) – The `anndata` object of gene expression data.
  * **database** (`str`) – Choose between ‘scent_17’ and ‘scent_13’ for the SCENT network updated in 2017 or 2013 respectively.
  * **species** (`str`) – ‘human’ or ‘mouse’
* **Returns:**
  **gene_network** – A gene network stored in `adata.uns[gene_network]`. The adjacency matrix of the network is in `adata.uns[gene_network]['network']`.
  The genes are in `adata.uns[gene_network]['genes']`.
* **Return type:**
  dict

### References

[Banerji2013] Banerji, C. R., Miranda-Saavedra, D., Severini, S., Widschwendter, M., Enver, T., Zhou, J. X., & Teschendorff, A. E. (2013). Cellular network entropy as the energy potential in Waddington’s differentiation landscape. Scientific reports, 3(1), 3039.
