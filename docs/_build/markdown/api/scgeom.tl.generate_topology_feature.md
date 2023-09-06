# scgeom.tl.generate_topology_feature

### scgeom.tl.generate_topology_feature(adata, method='tp', inf_value='replace_max', pi_pixel_sizes=[0.5, 0.5], dims=[0, 1], pi_return_images=False, diagram_name=None, feature_name=None, pi_birth_ranges=None, pi_pers_ranges=None, bc_n_intervals=[30, 30], bc_v_mins=[None, None], bc_v_maxs=[None, None])

Construct structured features from the unstructured output of persistent homology (persistence diagrams).

* **Parameters:**
  * **adata** (`AnnData`) – The `anndata` object of the gene expression data.
  * **method** (`str`) – The method to use for computing features.
    ‘tp’ for total persistence (summation of death - birth values),
    ‘std’ for standard deviation of persistences,
    ‘pe’ for persistence entropy,
    ‘bc’ for Betti curves,
    ‘pi’ for persistence images.
  * **inf_value** (`str`) – if `inf_value` is set to ‘replace_max’, infinities in persistence diagrams will be replaced by the max filtration values before computing the features.
  * **pi_pixel_sizes** (`list`) – The pixel sizes to use for each dimension for persistence images. Should have the same length with `dims`.
  * **dims** (`list`) – Compute features for which dimensions of persistence diagrams.
  * **pi_return_images** (`bool`) – Whether to return the persistence image results.
  * **diagram_name** (`Optional`[`str`]) – The diagram to use for cell ‘A’ should be available at `adata.uns[diagram_name]['A']`.
  * **feature_name** (`Optional`[`str`]) – The computed features will be stored in `adata.obsm[feature_name]`.
  * **pi_birth_ranges** (`Optional`[`list`]) – The birth ranges to use in persistence images for each dimension. Should have the same length with `dims`. If `None`, it will be determined by the PI algorithm.
  * **pi_pers_ranges** (`Optional`[`list`]) – The persistence ranges to use in persistence images for each dimension. Should have the same length with `dims`. If `None`, it will be determined by the PI algorithm.
  * **bc_n_intervals** (`list`) – The number of equal-length intervals to compute Betti curves. Should have the same length with `dims`.
  * **bc_v_mins** (`list`) – The left ends to compute Betti curves. Should have the same length with `dims`.
  * **bc_v_maxs** (`list`) – The right ends to compute Betti curves. Should have the same length with `dims`.
* **Returns:**
  **feature** – A `n_cell x n_feature` matrix in `adata.obsm[feature_name]`.
* **Return type:**
  np.ndarray
