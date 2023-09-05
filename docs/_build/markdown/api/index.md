<!-- commot documentation master file, created by
sphinx-quickstart on Sat Feb 20 12:08:49 2021.
You can adapt this file completely to your liking, but it should at least
contain the root `toctree` directive. -->

# API

## Preprocessing: pp

| [`pp.get_gene_network`](scgeom.pp.get_gene_network.md#scgeom.pp.get_gene_network)(adata[, database, species])                 | Get the gene network from knowledge database.                                |
|-------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------|
| [`pp.infer_gene_network`](scgeom.pp.infer_gene_network.md#scgeom.pp.infer_gene_network)(adata[, method, ...])                 | Simple computation of gene networks with Spearman's correlation coefficient. |
| [`pp.cell_specific_gene_network`](scgeom.pp.cell_specific_gene_network.md#scgeom.pp.cell_specific_gene_network)(adata[, ...]) | Compute a cell-specific gene network for each cell.                          |
