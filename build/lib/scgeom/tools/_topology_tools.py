from platform import node
from re import I
import numpy as np
import anndata
import gudhi
import networkx as nx
from scipy import sparse
from scipy.sparse import csr_matrix
from scipy.stats import norm
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist, squareform
from ripser import ripser
import dionysus as dys
import matplotlib.pyplot as plt
from multiprocessing import Process, Manager

from .._persistent_homology import ph_process_vbcl
from .._persistent_homology import ph_process_lwph
from .._persistent_homology import ph_process_rlph
from .._persistent_homology import ph_process_ewvr
from .._persistent_homology import ph_feature_persimages
from .._persistent_homology import ph_feature_total_persistence
from .._persistent_homology import ph_feature_persistence_std
from .._persistent_homology import ph_feature_persistence_entropy
from .._persistent_homology import ph_feature_betti_curve

def gene_network_topology(
    adata: anndata.AnnData,
    network_name: str = 'gene_network',
    filtration: str = 'vbcl',
    max_filtration: float = None,
    max_dim: int = 1,
    n_cores: int = 20,
    vbcl_edge_weight: str = 'min',
    return_diagrams: bool = False,
    uns_name: str = 'diagrams',
):  
    """
    Compute gene network topology for node weighted gene networks.

    Parameters
    ----------
    adata
        The ``anndata`` object of gene expression data.
    network_name
        The network structure to be used should be available in ``adata.uns[network_name]``. Use ``pp.get_gene_network`` to obtain a knowledge-based gene network to be shared by each cell.
    filtration
        The filtration method for computing persistent homology. 'vbcl' for vertex-based clique complex filtration. Currently, only 'vbcl' is implemented.
    max_filtration
        The maximum filtration value. If ``None``, it is set to the maximum node weight.
    max_dim
        The maximum dimension for computing persistent homology.
    n_cores
        Number of cores to use for computing persistent homology for each cell.
    vbcl_edge_weight
        How filtration value for edges are assigned. 
        'min' will assign the smaller node weight to the edge. 
        'mean' will use the average of two node weights for each edge.
    return_diagrams
        Whether to return the computed persistence diagrams.
    uns_name
        The persistence diagrams are stored in ``adata.uns[uns_name]``.
    
    Returns
    -------
    diagrams: dict
        A dictionary of diagrams to be stored in ``adata.uns[uns_name]``.
        For example, the H1 diagram for the cell named 'A' can be accessed through ``diagrams['A']['dgm_h1]`` 
        and the corresponding max_filtration value is in ``diagrams['A']['max_filtration']``.

    """

    manager = Manager()

    genes = adata.uns[network_name]['genes']
    A = adata.uns[network_name]['network']
    A = sparse.coo_matrix(A, dtype=float)
    X_gnet = np.array( adata[:,genes].X )
    ncell = adata.shape[0]
    # ngene = X_gnet.shape[1]
    persimages = []
    diagrams = manager.dict()
    D_maxs = []
    batch_size = n_cores
    cell_names = list(adata.obs_names)
    if filtration == 'vbcl':
        for i in range(0, ncell, batch_size):
            print(i)
            processes = [Process(target=ph_process_vbcl, \
                args=(A,X_gnet,j,diagrams,max_filtration,max_dim,vbcl_edge_weight,cell_names)) \
                for j in range(i, min(i+batch_size,ncell))]
            for process in processes:
                process.start()
            for process in processes:
                process.join()
        print('Done', flush=True)
    # elif filtration == 'wvr':
    #     for i in range(0, ncell, batch_size):
    #         print(i)
    #         processes = [Process(target=ph_process_wvr, args=(A,X_gnet,j,diagrams,max_filtration,max_dim)) for j in range(i, min(i+batch_size,ncell))]
    #         for process in processes:
    #             process.start()
    #         for process in processes:
    #             process.join()
    #     print('Done', flush=True)
    
    adata.uns[uns_name] = dict( diagrams )

    if return_diagrams:
        return diagrams
    else:
        return 

def cell_specific_gene_network_topology(
    adata: anndata.AnnData,
    network_name: str = 'cell_specific_gene_network',
    filtration: str = 'ewvr',
    max_filtration: float = None,
    max_dim: int = 1,
    n_cores: int = 20,
    return_diagrams: bool = False,
    uns_name: str = 'diagrams',
    compute_csn: bool = False,
    csn_alpha: float = 0.01,
    csn_boxsize: float = 0.1,
    csn_weighted: float = True,
    print_every_n_cell: str = 100,
):
    """
    Compute gene network topology for edge-weighted networks such as the cell-specific network.

    Parameters
    ----------
    adata
        The ``adata`` object of gene expression data.
    network_name
        The cell-specific network for the cell with cell name 'A' should be available in ``adata.uns[network_name]['A']``.
    filtration
        The filtration method for computing persistent homology. 'ewvr' for edge-weighted Vietoris-Rips complex filtration. Currently, only 'ewvr' is implemented.
    max_filtration
        The maximum filtration value. If ``None``, it is set to the maximum node weight.
    max_dim
        The maximum dimension for computing persistent homology.
    n_cores
        Number of cores to use for computing persistent homology for each cell.
    return_diagrams
        Whether to return the computed persistence diagrams.
    uns_name
        The persistence diagrams are stored in ``adata.uns[uns_name]``.
    compute_csn
        Whether to compute cell-specific network here. If ``True``,  precomputed CSN in ``adata.uns[network_name]`` is ignored.
    csn_alpha
        The alpha parameter in CSN computation, i.e., the p-value threshold for including a connection.
    csn_boxsize
        The boxsize for defining neighborhoods of gene expression in CSN method.
    csn_weighted
        Whether to output the weighted cell-specific networks.

    Returns
    -------
    diagrams: dict
        A dictionary of diagrams to be stored in ``adata.uns[uns_name]``.
        For example, the H1 diagram for the cell named 'A' can be accessed through ``diagrams['A']['dgm_h1]`` 
        and the corresponding max_filtration value is in ``diagrams['A']['max_filtration']``.    

    """
    
    manager = Manager()

    ncell = adata.shape[0]
    diagrams = manager.dict()
    batch_size = n_cores
    cell_names = list(adata.obs_names)
    if filtration == 'ewvr':
        if compute_csn:
            data = np.array( adata.X ).T
            weighted = csn_weighted
            boxsize = csn_boxsize
            alpha = csn_alpha
            n1, n2 = data.shape
            c = np.arange(n2)
            upper = np.zeros((n1, n2))
            lower = np.zeros((n1, n2))
            for i in range(n1):
                s1 = np.sort(data[i, :])
                s2 = np.argsort(data[i,:])
                n3 = int(n2 - np.sum(np.sign(s1)))
                h = round_half_up(boxsize / 2 * np.sum(np.sign(s1)))
                k = 0
                while k < n2:
                    s = 0
                    while k + s + 1 < n2 and s1[k + s + 1] == s1[k]:
                        s += 1
                    if s >= h:
                        upper[i, s2[k:k + s + 1]] = data[i, s2[k]]
                        lower[i, s2[k:k + s + 1]] = data[i, s2[k]]
                    else:
                        upper[i, s2[k:k + s + 1]] = data[i, s2[min(n2-1, k + s + h)]]
                        lower[i, s2[k:k + s + 1]] = data[i, s2[max(n3 * int(n3 > h), k - h)]]
                    k = k + s + 1
            p = -norm.ppf(alpha, 0, 1)
        tmp_cnt = 0
        for i in range(0, ncell, batch_size):
            if tmp_cnt * print_every_n_cell <= i:
                print("cell ", i)
                tmp_cnt += 1
            if not compute_csn:
                processes = [Process(target=ph_process_ewvr, \
                    args=(adata.uns[network_name][cell_names[j]],j,diagrams,max_filtration,max_dim,cell_names)) \
                    for j in range(i, min(i+batch_size,ncell))]
            elif compute_csn:
                csn = {}
                for k in range(i, min(i+batch_size, ncell)):
                    B = np.zeros((n1, n2))
                    for j in range(n2):
                        B[:, j] = (data[:, j] <= upper[:, k]) & (data[:, j] >= lower[:, k])
                    a = np.sum(B, axis=1)
                    d = (B.dot(B.T) * n2 - a.reshape(-1,1).dot(a.reshape(1,-1))) / np.sqrt((a.reshape(-1,1).dot(a.reshape(1,-1))) * ((n2 - a).reshape(-1,1)).dot((n2 - a).reshape(1,-1)) / (n2 - 1) + np.finfo(float).eps)
                    np.fill_diagonal(d, 0)
                    if weighted:
                        tmp_mat = csr_matrix( d * (d > 0) )
                        csn[k] = tmp_mat
                    else:
                        tmp_mat = csr_matrix((d > p).astype(float))
                        csn[k] = tmp_mat

                processes = [Process(target=ph_process_ewvr, \
                    args=(csn[j],j,diagrams,max_filtration,max_dim,cell_names)) \
                    for j in range(i, min(i+batch_size,ncell))]
                
            for process in processes:
                process.start()
            for process in processes:
                process.join()
        print('Done', flush=True)

    adata.uns[uns_name] = dict( diagrams )
    
    if return_diagrams:
        return diagrams
    else:
        return



def cell_network_topology(
    adata: anndata.AnnData,
    network_name: str = 'connectivities',
    embedding_name: str = 'X_pca',
    method: str = 'lwph',
    rlph_method: str = 'self',
    rlph_base: str = 'all',
    metric: str = 'euclidean',
    nb_method: str = 'knn',
    nb_distance: float = 20,
    rl_distance: float = 10,
    nb_knn: int = 20,
    rl_knn: int = 10,
    filtration_value: str = 'euclidean',
    max_filtration: float = np.inf,
    max_dim: int = 1,
    n_cores: int = 20,
    uns_name: str = 'cell_network_diagrams',
    return_diagrams: bool = False,
):
    """
    Compute topology for each cell on the cell network.

    Parameters
    ----------
    adata
        The ``adata`` object of gene expression data.
    network_name
        The network of cells, e.g. a knn network to be used in ``adata.obsp[network_name]``.
    embedding_name
        The low-dimensional embedding of cells to be used to compute filtration values, stored in ``adata.obsm[embedding_name]``.
    method
        The persistent homology method to use.
        'lwph' for local weighted persistent homology which computes Vietoris-Rips complex based PH on a local neighborhood of each cell.
        'rlph' for relative local persistent homology which computes the relative persistent homology of the cell, i.e., the global topology relative to a cell or its local neighborhood.
    rlph_method
        When ``method`` is set to 'rlph', choose 'self' to compute persistent homology relative only to the cell, 
        'knn' for relative to a neighborhood defined by k-nearest neighbors, or 'distance' for relative to a neighborhood defined by distance cutoff.
    rlph_base
        What points to include when computing the global topology. 'all' for including all points, 
        'nb_knn' or 'nb_distance' for a relatively large neighborhood defined by k-nearest neighbors or distance cutoff, respectively.
    metric
        The metric used to define base point set for computing global topology when ``rlph_base`` is set to 'nb_knn' or 'nb_distance'.
    nb_method
        The method to define neighborhood of cells when ``method`` is set to 'lwph'. 'knn' for k-nearest neighbors or 'distance' for distance cutoff.
    nb_distance
        If ``nb_method`` is set to 'distance' in 'lwph' or ``rlph_base`` is set to 'nb_distance' in 'rlph', the distance cutoff to use.
    rl_distance
        If ``method`` is set to 'rlph' and ``rlph_method`` is set to 'distance', the distance cutoff to use to define the local cell neighborhood.
    nb_knn
        If ``nb_method`` is set to 'knn' in 'lwph' or ``rlph_base`` is set to 'nb_knn' in 'rlph', the k value to use.
    rl_knn
        If ``method`` is set to 'rlph' and ``rlph_method`` is set to 'knn', the k value to use to define the local cell neighborhood.
    filtration_value
        When 'lwph' is used, the filtration value to use. Currently, only Euclidean distance in the low-dimensional embedding is used.
    max_filtration
        The maximum filtration value. If ``None``, it is set to the maximum node weight.
    max_dim
        The maximum dimension for computing persistent homology.
    n_cores
        Number of cores to use for computing persistent homology for each cell.
    return_diagrams
        Whether to return the computed persistence diagrams.
    uns_name
        The persistence diagrams are stored in ``adata.uns[uns_name]``.
    
    Returns
    -------
    diagrams: dict
        A dictionary of diagrams to be stored in ``adata.uns[uns_name]``.
        For 'lwph', the H1 diagram for the cell named 'A' can be accessed through ``diagrams['A']['dgm_h1]`` 
        and the corresponding max_filtration value is in ``diagrams['A']['max_filtration']``. 
        
        For 'rlph', the H1 relative PH results for base point set relative to local point set is in ``diagrams['A']['dgm_h1']``.
        Additionally, the H1 regular PH results for the base point set is in ``diagrams['A']['dgm_h1_base']``.

    """

    A = adata.obsp[network_name]
    X = adata.obsm[embedding_name]
    D = distance_matrix(X,X)

    manager = Manager()

    ncell = A.shape[0]
    diagrams = manager.dict()
    batch_size = n_cores
    cell_names = list(adata.obs_names)

    if method == 'lwph':
        for i in range(0, ncell, batch_size):
            print(i)
            processes = [Process(target=ph_process_lwph, \
                args=(A,D,j,diagrams,nb_method,nb_distance,nb_knn, \
                filtration_value,max_dim,max_filtration,cell_names)) \
                for j in range(i, min(i+batch_size,ncell))]
            for process in processes:
                process.start()
            for process in processes:
                process.join()
        print('Done', flush=True)
    
    elif method == 'rlph':
        dists = pdist(X, metric)
        f = dys.fill_rips(dists, max_dim+1, max_filtration)
        for i in range(0, ncell, batch_size):
            print(i)
            processes = [Process(target=ph_process_rlph, \
                args=(X, f, j, diagrams, metric, rlph_method, \
                rlph_base, nb_knn, rl_knn, nb_distance, rl_distance, \
                max_dim, max_filtration, cell_names)) \
                for j in range(i, min(i+batch_size,ncell))]
            for process in processes:
                process.start()
            for process in processes:
                process.join()
        print('Done', flush=True)

    adata.uns[uns_name] = dict( diagrams )

    if return_diagrams:
        return diagrams
    else:
        return 


def preprocess_persistence_diagrams(
    adata,
    diagram_name = None,
    inf_value = None,
    dims = [0,1],
    method = 'prune',
    prune_min_per = [0.1,0.1],
):
    if method == 'prune':
        ph_prune_diagrams(adata,
            diagram_name = diagram_name,
            dims = dims,
            inf_value = inf_value,
            min_per = prune_min_per)
    return

def generate_topology_feature(
    adata: anndata.AnnData,
    method: str = 'tp',
    inf_value: str = 'replace_max',
    pi_pixel_sizes: list = [0.5, 0.5],
    dims: list = [0,1],
    pi_return_images: bool = False,
    diagram_name: str = None,
    feature_name: str = None,
    pi_birth_ranges: list = None,
    pi_pers_ranges: list = None,
    bc_n_intervals: list = [30,30],
    bc_v_mins: list = [None, None],
    bc_v_maxs: list = [None, None],
):
    """
    Construct structured features from the unstructured output of persistent homology (persistence diagrams).

    Parameters
    ----------
    adata
        The ``anndata`` object of the gene expression data.
    method
        The method to use for computing features.
        'tp' for total persistence (summation of death - birth values),
        'std' for standard deviation of persistences,
        'pe' for persistence entropy,
        'bc' for Betti curves,
        'pi' for persistence images.
    inf_value
        if ``inf_value`` is set to 'replace_max', infinities in persistence diagrams will be replaced by the max filtration values before computing the features.
    pi_pixel_sizes
        The pixel sizes to use for each dimension for persistence images. Should have the same length with ``dims``.
    dims
        Compute features for which dimensions of persistence diagrams.
    pi_return_images
        Whether to return the persistence image results.
    diagram_name
        The diagram to use for cell 'A' should be available at ``adata.uns[diagram_name]['A']``.
    feature_name
        The computed features will be stored in ``adata.obsm[feature_name]``.
    pi_birth_ranges
        The birth ranges to use in persistence images for each dimension. Should have the same length with ``dims``. If ``None``, it will be determined by the PI algorithm.
    pi_pers_ranges
        The persistence ranges to use in persistence images for each dimension. Should have the same length with ``dims``. If ``None``, it will be determined by the PI algorithm.
    bc_n_intervals
        The number of equal-length intervals to compute Betti curves. Should have the same length with ``dims``.
    bc_v_mins
        The left ends to compute Betti curves. Should have the same length with ``dims``.
    bc_v_maxs
        The right ends to compute Betti curves. Should have the same length with ``dims``.

    Returns
    -------
    feature: np.ndarray
        A ``n_cell x n_feature`` matrix in ``adata.obsm[feature_name]``.
    """
    ncell = adata.shape[0]
    if method == 'pi':
        images, persimage_info = ph_feature_persimages(adata,
            inf_value = inf_value,
            pixel_sizes = pi_pixel_sizes,
            dims = dims,
            diagram_name = diagram_name,
            birth_ranges = pi_birth_ranges,
            pers_ranges = pi_pers_ranges)
        feature = np.empty([ncell,0])
        for idim in range(len(dims)):
            image = images[idim]
            feature = np.concatenate((feature, image.reshape(ncell,-1)), axis=1)
        adata.uns[feature_name] = persimage_info
        adata.obsm[feature_name] = feature
        if pi_return_images and method == 'pi':
            return images, persimage_info
        else:
            return
    elif method == 'pe':
        X_pe = ph_feature_persistence_entropy(adata,
            diagram_name = diagram_name,
            dims = dims,
            inf_value = inf_value)
        if X_pe.shape[1] == 1:
            adata.obs[feature_name] = X_pe.reshape(-1)
        else:
            adata.obsm[feature_name] = X_pe
    elif method == 'tp':
        X_tp = ph_feature_total_persistence(adata,
            diagram_name = diagram_name,
            dims = dims,
            inf_value = inf_value)
        if X_tp.shape[1] == 1:
            adata.obs[feature_name] = X_tp.reshape(-1)
        else:
            adata.obsm[feature_name] = X_tp
    elif method == 'std':
        X_std = ph_feature_persistence_std(adata,
            diagram_name = diagram_name,
            dims = dims,
            inf_value = inf_value)
        if X_std.shape[1] == 1:
            adata.obs[feature_name] = X_std.reshape(-1)
        else:
            adata.obsm[feature_name] = X_std
    elif method == 'bc':
        X_bc_tmp, xs_list = ph_feature_betti_curve(adata,
            diagram_name = diagram_name,
            dims = dims,
            inf_value = inf_value,
            n_intervals = bc_n_intervals,
            v_mins = bc_v_mins,
            v_maxs = bc_v_maxs)
        X_bc = np.concatenate(tuple([X_bc_tmp[i] for i in range(len(dims))]), axis=1)
        adata.obsm[feature_name] = X_bc
        adata.uns[feature_name] = {'xs': xs_list}
    return
    

































# =============================================================================

# def local_persistent_homology(
#     adata,
    
# )

def gene_network_persistent_homology(
    adata,
    network_name = 'gene_network',
    filtration = 'rips_1',
):  
    """
    rips_1: 
    """
    genes = adata.uns[network_name]['genes']
    A_network = adata.uns[network_name]['network']
    X_gnet = np.array( adata[:,genes].X )
    for icell in range(adata.shape[0]):
        nd_weight = X_gnet[icell,:]
        # if filtration == 'rips_1':
            
    return

def persistent_homology(
    adata,
    graph_name = None,
    node_filtration = None,
    dtm_p = 2,
    ph_filtration = 'dtm',
    ph_thresh = np.inf,
    ph_plot = False,
    ph_name = 'scgeom-ph'
):
    G = nx.from_scipy_sparse_matrix( adata.obsp[graph_name] )
    f = np.array( adata.obs[node_filtration] ).reshape(-1)
    diagram, D_rips = value_based_persistent_homology(G, f,
        dtm_p=dtm_p, thresh=ph_thresh, plot=ph_plot, filtration=ph_filtration)
    adata.uns[ph_name] = {'diagram': diagram, 'D': D_rips}

def persistent_homology_density(
    adata,
    D,
    node_filtration = None,
    ph_k = 15,
    ph_p = 2,
    ph_thresh = np.inf,
    ph_plot = False,
    ph_name = 'scgeom-ph'
):
    f = np.array(adata.obs[node_filtration]).reshape(-1)
    diagram, D_rips = persistent_homology_dtm(D, f, p=ph_p, thresh=ph_thresh,
        k=ph_k, plot=ph_plot)
    adata.uns[ph_name] = {'diagram': diagram, 'D': D_rips}

def persistent_homology_tomato(
    adata,
    G,
    node_filtration = None,
    ph_name = 'scgeom-ph',
    ph_plot = False
):
    f = np.array(adata.obs[node_filtration]).reshape(-1)
    n = len(f)

    sxt = gudhi.SimplexTree()
    
    for ind in range(n):
        sxt.insert([ind], filtration=f[ind])
        nei = list(G.neighbors(ind))
        for idx in nei:
            sxt.insert([ind, idx], filtration=np.max([f[ind], f[idx]]))
    
    sxt.initialize_filtration()
    sxt.persistence()
    dig, res = sxt.persistence(), []
    bars = []
    diagram = []
    for ele in dig:
        if ele[0] == 0:
            birth = -ele[1][0]
            death = -ele[1][1]
            bar = (0, (birth, death))
            bars.append([birth, death])
            diagram.append(bar)
            res.append(bar)
    bars = np.asarray(bars)
    if ph_plot:
        for i in range(bars.shape[0]):
            plt.plot([bars[i,0]], [max(0,bars[i,1])], "bo")
        plt.plot([0,np.max(bars[:,0])], [0,np.max(bars[:,0])], 'k')
        plt.show()

    adata.uns[ph_name] = {'diagram': diagram}

def clustering_topology(
    adata,
    tau,
    graph_name = None,
    node_filtration = None, # should be density
):
    G = nx.from_scipy_sparse_matrix(adata.obsp[graph_name])
    f = np.array( adata.obs[node_filtration] ).reshape(-1)
    n = len(G)
    G_srt = nx.Graph()
    G_srt.add_nodes_from([i for i in range(n)])
    idmap_srt2org = np.argsort(-f)
    idmap_org2srt = {}
    for i in range(len(idmap_srt2org)):
        idmap_org2srt[idmap_srt2org[i]] = i
    org_edges = list(G.edges())
    for org_edge in org_edges:
        G_srt.add_edge(idmap_org2srt[org_edge[0]], idmap_org2srt[org_edge[1]])
    f_srt = -np.sort(-f)

    U = UnionFind()
    r = np.empty([n], int)
    g = np.empty([n], int)
    e = -1
    for i in range(n):
        nb_ind = np.array(list(G_srt.neighbors(i)), int)
        nb_low_ind = nb_ind[np.where(nb_ind < i)[0]]
        # print(i,nb_ind,nb_low_ind)
        if nb_low_ind.shape[0] == 0:
            U.insert_objects([i])
            r[i] = i
        else:
            g[i] = np.min(nb_low_ind)
            ei = U.find(g[i])
            U.union(ei, i)
            for j in nb_low_ind:
                e = U.find(j)
                if e!=ei and min(f_srt[r[e]], f_srt[r[ei]]) < f_srt[i] + tau:
                    U.union(e,ei)
                    ei = U.find(e)
                    r[ei] = ei
    label_srt = np.empty([n], int)
    for i in range(n):
        label_srt[i] = U.find(i)
    label = np.empty_like(label_srt)
    for i in range(n):
        label[i] = idmap_srt2org[label_srt[idmap_org2srt[i]]]

    return label

def clustering_persistent_homology(
    adata,
    graph_name = None,
    node_filtration = None,
    ph_name = None,
    persistence = None
):
    G = nx.from_scipy_sparse_matrix(adata.obsp[graph_name])
    n = len(G)
    D = adata.uns[ph_name]['D']
    f = np.array( adata.obs[node_filtration] ).reshape(-1)
    G_srt = nx.Graph()
    G_srt.add_nodes_from(i for i in range(n))
    idmap_srt2org = np.argsort(f)
    idmap_org2srt = {}
    for i in range(len(idmap_srt2org)):
        idmap_org2srt[idmap_srt2org[i]] = i
    org_edges = list(G.edges())
    for org_edge in org_edges:
        i = idmap_org2srt[org_edge[0]]
        j = idmap_org2srt[org_edge[1]]
        G_srt.add_edge(i, j, weight=D[org_edge[0], org_edge[1]])
    f_srt = -np.sort(-f)

    U = UnionFind()
    r = np.empty([n], int)
    g = np.empty([n], int)
    e = -1
    for i in range(n):
        nb_ind = np.array(list(G_srt.neighbors(i)), int)
        nb_low_ind = nb_ind[np.where(nb_ind < i)[0]]
        # print(i,nb_ind,nb_low_ind)
        if nb_low_ind.shape[0] == 0:
            U.insert_objects([i])
            r[i] = i
        else:
            g[i] = np.min(nb_low_ind)
            ei = U.find(g[i])
            U.union(ei, i)
            for j in nb_low_ind:
                e = U.find(j) 
                value = max(f_srt[r[e]], f_srt[r[ei]])
                if e!=ei and  value < f_srt[i] + persistence:
                    U.union(e,ei)
                    ei = U.find(e)
                    r[ei] = ei
    label_srt = np.empty([n], int)
    for i in range(n):
        label_srt[i] = U.find(i)
    label = np.empty_like(label_srt)
    for i in range(n):
        label[i] = idmap_srt2org[label_srt[idmap_org2srt[i]]]

    return label

def confidence_band(
    adata,
    graph_name = None,
    node_filtration = None,
    dtm_p = 2,
    ph_filtration = 'dtm',
    ph_thresh = np.inf,
    ph_plot = False,
    ph_name = 'scgeom-ph',
    repeat = 200,
    alpha=0.05,
    wasserstein_p = np.inf,
    random_state=0
):
    np.random.seed(random_state)
    G = nx.from_scipy_sparse_matrix( adata.obsp[graph_name] )
    f = np.array( adata.obs[node_filtration] ).reshape(-1)
    dgm_base, _ = value_based_persistent_homology(G, f,
        dtm_p=dtm_p, thresh=ph_thresh, plot=ph_plot, filtration=ph_filtration)
    dgm_base_0 = []
    for bar in dgm_base:
        if bar[0] == 0:
            dgm_base_0.append([bar[1][0], bar[1][1]])
    n = len(f)
    ind = np.arange(n)
    A = adata.obsp[graph_name]
    diff_0 = []
    for i in range(repeat):
        shuffled_ind = np.unique(np.random.choice(ind, n))
        A_sub = sparse.csr_matrix( A.toarray()[shuffled_ind,:][:,shuffled_ind] )
        G_sub = nx.from_scipy_sparse_matrix( A_sub )
        dgm_apx,_ = value_based_persistent_homology(G_sub, f[shuffled_ind],
            dtm_p=dtm_p, thresh=ph_thresh, plot=ph_plot, filtration=ph_filtration)
        dgm_apx_0 = []
        for bar in dgm_apx:
            if bar[0] == 0:
                dgm_apx_0.append([bar[1][0], bar[1][1]])
        d_0 = gudhi.bottleneck_distance(dgm_base_0, dgm_apx_0)
        diff_0.append(d_0)
    diff_0 = np.array(diff_0, float)
    return np.quantile(diff_0, 1-alpha)

def decide_ncluster(
    diagram,
    min_persistence = 10
):
    n = 0
    for i in range(len(diagram)):
        tmp = diagram[i]
        dim = tmp[0]
        birth = tmp[1][0]
        death = tmp[1][1]
        if np.abs(death-birth) >= min_persistence:
            n += 1
    return n

class UnionFind:
    # http://code.activestate.com/recipes/215912-union-find-data-structure/
    # Initialization
    def __init__(self):

        self.num_weights = {}
        self.parent_pointers = {}
        self.num_to_objects = {}
        self.objects_to_num = {}
        
    # Insert objects among the already existing ones
    def insert_objects(self, objects):

        for object in objects: self.find(object)
            
    # Find a given object / build it if non-existing
    def find(self, object):

        if not object in self.objects_to_num:

            obj_num = len(self.objects_to_num)
            self.num_weights[obj_num] = 1
            self.objects_to_num[object] = obj_num
            self.num_to_objects[obj_num] = object
            self.parent_pointers[obj_num] = obj_num

            return object
        
        stk = [self.objects_to_num[object]]
        par = self.parent_pointers[stk[-1]]

        while par != stk[-1]:
            stk.append(par)
            par = self.parent_pointers[par]

        for i in stk: self.parent_pointers[i] = par
            
        return self.num_to_objects[par]
    
    # Link two different objects in a same distinct set
    def union(self, object1, object2):

        o1p = self.find(object1)
        o2p = self.find(object2)
        
        if o1p != o2p:
        	
            on1 = self.objects_to_num[o1p]
            on2 = self.objects_to_num[o2p]
            w1 = self.num_weights[on1]
            w2 = self.num_weights[on2]

            if w1 < w2: o1p, o2p, on1, on2, w1, w2 = o2p, o1p, on2, on1, w2, w1

            self.num_weights[on1] = w1+w2
            del self.num_weights[on2]
            self.parent_pointers[on2] = on1

def round_half_up(x):
    if x == 0.5:
        return 1
    else:
        return int(round(x))