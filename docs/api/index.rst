.. commot documentation master file, created by
   sphinx-quickstart on Sat Feb 20 12:08:49 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. module:: scgeom
.. automodule:: scgeom
   :noindex:

API
==================================

Preprocessing: pp
-------------------

.. module:: scgeom.pp
.. currentmodule:: scgeom

.. autosummary::
   :toctree: .

   pp.get_gene_network
   pp.infer_gene_network
   pp.cell_specific_gene_network


Tools: tl
---------

.. module:: scgeom.tl
.. currentmodule:: scgeom

.. autosummary::
   :toctree: .
  
   tl.neighbor_graph
   tl.graph_curvature
   tl.gene_network_topology
   tl.cell_specific_gene_network
   tl.cell_network_topology
   tl.generate_topology_feature
