import networkx as nx
import numpy as np
from sklearn.neighbors import kneighbors_graph

def knn_graph_nx(D,k):
    """Construct a k-nearest-neighbor graph as networkx object.

    :param D: a distance matrix for constructing the knn graph
    :type D: class:`numpy.ndarray`
    :param k: number of nearest neighbors
    :type k: int
    :return: a knn graph object and a list of edges
    :rtype: class:`networkx.Graph`, list of tuples
    """
    G_nx = nx.Graph()
    G_nx.add_nodes_from([i for i in range(D.shape[0])])
    edge_list = []
    for i in range(D.shape[0]):
        tmp_ids = np.argsort(D[i,:])
        for j in range(1,k+1):
            G_nx.add_edge(i,tmp_ids[j])
            edge_list.append((i,tmp_ids[j]))
    return G_nx, edge_list

