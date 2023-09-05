import numpy as np
import networkx as nx
from scipy import sparse
from sklearn.neighbors import kneighbors_graph

from GraphRicciCurvature.OllivierRicci import OllivierRicci
from GraphRicciCurvature.FormanRicci import FormanRicci

def neighbor_graph(
    adata,
    X = None,
    D = None,
    graph_name = 'scgeom-graph',
    graph_method = 'knn',
    knn_k = 10,
    weighted = True,
    D_type = 'distance'
):
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
    adata,
    graph_name = 'scgeom-graph',
    curvature_method = 'orc',
    orc_alpha = 0.5,
    orc_method = 'OTD',
    node_curvature_name = 'scgeom-node_curvature',
    edge_curvature_name = 'scgeom-edge_curvature'
):
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