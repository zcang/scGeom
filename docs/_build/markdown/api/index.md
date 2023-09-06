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

## Tools: tl

| [`tl.neighbor_graph`](scgeom.tl.neighbor_graph.md#scgeom.tl.neighbor_graph)(adata[, X, D, graph_name, ...])                | Construct a cell graph for curvature computation.                                                         |
|----------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------|
| [`tl.graph_curvature`](scgeom.tl.graph_curvature.md#scgeom.tl.graph_curvature)(adata[, graph_name, ...])                   | Compute graph curvature for each cell.                                                                    |
| [`tl.gene_network_topology`](scgeom.tl.gene_network_topology.md#scgeom.tl.gene_network_topology)(adata[, ...])             | Compute gene network topology for node weighted gene networks.                                            |
| [`tl.cell_network_topology`](scgeom.tl.cell_network_topology.md#scgeom.tl.cell_network_topology)(adata[, ...])             | Compute topology for each cell on the cell network.                                                       |
| [`tl.generate_topology_feature`](scgeom.tl.generate_topology_feature.md#scgeom.tl.generate_topology_feature)(adata[, ...]) | Construct structured features from the unstructured output of persistent homology (persistence diagrams). |
