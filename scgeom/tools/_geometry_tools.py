import numpy as np
import networkx as nx
import anndata
from scipy import sparse
from sklearn.neighbors import kneighbors_graph

from GraphRicciCurvature.OllivierRicci import OllivierRicci
from GraphRicciCurvature.FormanRicci import FormanRicci

def neighbor_graph(
    adata: anndata.AnnData,
    X: np.array = None,
    D: np.array = None,
    graph_name: str = 'scgeom-graph',
    graph_method: str = 'knn',
    knn_k: int = 10,
    weighted: bool = True,
    D_type: str = 'distance'
):
    """
    Construct a cell graph for curvature computation.

    Parameters
    ----------
    adata
        The ``anndata`` object of gene expression data.
    X
        A low dimensional embedding of data of shape ``n_cell x n_dim`` where Euclidean distance will be used to quantify similarity. 
        If ``X`` is given, ``D`` will be ignored.
    D
        A ``n_cell x n_cell`` distance or similarity matrix. If ``X`` is given, ``D`` will be ignored.
    graph_name
        The name to store the cell graph in ``adata.obsp[graph_name]``.
    graph_method
        Currently, only 'knn' (k-nearest neighbor graph) is implemented.
    knn_k
        The number of neighbors when buidling the knn graph.
    weighted
        Whether to include weights in the sparse matrix representing the graph.
    D_type
        'distance' is ``D`` is a distance/dissimilarity matrix or 'similarity' is ``D`` represents similarity among points.
    
    Returns
    -------
    A: sp.coo_matrix
        A sparse matrix representing the cell graph stored in ``adata.obsp[graph_name]``.

    """
    if graph_method == 'knn':
        if not X is None:
            A = kneighbors_graph(X, knn_k, mode='connectivity',
                metric='minkowski', p=2, include_self=False)
        elif not D is None:
            data = []; i = []; j = []
            for ii in range(D.shape[0]):
                if D_type == 'distance':
                    tmp_idx = np.argsort(D[ii,:])
                elif D_type == 'similarity':
                    tmp_idx = np.argsort(-D[ii,:])
                for jj in range(1,knn_k+1):
                    if weighted:
                        data.append(D[ii,tmp_idx[jj]])
                    else:
                        data.append(1.0)
                    i.append(ii)
                    j.append(tmp_idx[jj])
            A = sparse.coo_matrix((data,(i,j)), shape=(D.shape[0], D.shape[0]))
            A = A.tocsr()

    adata.obsp[graph_name] = A

def graph_curvature(
    adata: anndata.AnnData,
    graph_name: str = 'scgeom-graph',
    curvature_method: str = 'orc',
    orc_alpha: float = 0.5,
    orc_method: str = 'OTD',
    node_curvature_name: str = 'scgeom-node_curvature',
    edge_curvature_name: str = 'scgeom-edge_curvature'
):
    """
    Compute graph curvature for each cell. We use the package GraphRicciCurvature [Ni2019]_ to compute the curvatures.

    Parameters
    ----------
    adata
        The ``anndata`` object of gene expression data.
    graph_name
        The graph stored in ``adata.obsp[graph_name]`` will be used.
    curvature_method
        'orc' for Ollivier-Ricci curvature and 'frc' for Forman-Ricci curvature.
    orc_alpha
        The alpha paramter for computing Ollivier-Ricci curvature. The weight of mass to be put on the center node when defining a mass distribution of a node's neighborhood.
    orc_method
        The method for computing distance between mass distributions. 
        'OTD' for optimal transport distance (earth mover's distance), 
        'Sinkhorn' for entropy regularized approximation of optimal transport distance,
        'ATD' for average transportation distance.
    node_curvature_name
        The computed node curvature will be stored in ``adata.obs[node_curvature_name]``.
    edge_curvature_name
        The computed edge curvature will be stored in ``adata.obsp[edge_curvature_name]``.
    
    References
    ----------

    [Ni2019] Ni, C. C., Lin, Y. Y., Luo, F., & Gao, J. (2019). Community detection on networks with Ricci flow. Scientific reports, 9(1), 9984.

    """
    A = adata.obsp[graph_name]
    G = nx.from_scipy_sparse_array(A)
    if curvature_method == 'orc':
        orc = OllivierRicci(G, alpha=orc_alpha, method=orc_method, verbose='ERROR')
        orc.compute_ricci_curvature()
        edges,weights = zip(*nx.get_edge_attributes(orc.G,'ricciCurvature').items())
        ndweights = project_to_node(edges, weights, A.shape[0])
    elif curvature_method == 'frc':
        frc = FormanRicci(G, verbose='ERROR')
        frc.compute_ricci_curvature()
        edges,weights = zip(*nx.get_edge_attributes(frc.G,'formanCurvature').items())
        ndweights = project_to_node(edges, weights, A.shape[0])
    
    data = []; i = []; j = []
    for idx in range(len(edges)):
        data.append(weights[idx])
        i.append(edges[idx][0])
        j.append(edges[idx][1])
    C = sparse.coo_matrix((data, (i,j)), shape=(A.shape[0], A.shape[0]))

    adata.obs[node_curvature_name] = ndweights
    adata.obsp[edge_curvature_name] = C.tocsr()




def project_to_node(edges, weights, n):
    nedges = np.zeros([n])
    node_weights = np.zeros([n])
    for i in range(len(edges)):
        node_weights[edges[i][0]] += weights[i]
        node_weights[edges[i][1]] += weights[i]
        nedges[edges[i][0]] += 1
        nedges[edges[i][1]] += 1
    # node_weights = node_weights / nedges
    return node_weights    