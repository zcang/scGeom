from re import L
import numpy as np
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist, squareform
from scipy import sparse
from sklearn.neighbors import KDTree
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
import gudhi
import gudhi.representations
import matplotlib.pyplot as plt
from ripser import ripser
import networkx as nx
# import PersistenceImages.persistence_images as pimg
import dionysus as dys

from multiprocessing import Process, Manager

def ph_process_vbcl(A,
    X,
    icell,
    diagrams,
    max_filtration,
    max_dim, 
    edge_weight,
    cell_names,
):
    """
    Vertex based Clique-complex persistent homology
    """
    nd_weight = X[icell,:]
    D = A.copy()
    nd_weight_max = np.max(nd_weight)
    nd_weight = nd_weight_max - nd_weight
    for i in range(len(A.data)):
        if edge_weight == 'min':
            D.data[i] = max(nd_weight[A.row[i]], nd_weight[A.col[i]])
        elif edge_weight == 'mean':
            D.data[i] = 0.5 * (nd_weight[A.row[i]] + nd_weight[A.col[i]])

    st = gudhi.SimplexTree()
    for i in range(len(nd_weight)):
        st.insert([i], filtration=nd_weight[i])
    for i in range(len(D.data)):
        st.insert([D.row[i],D.col[i]], filtration=D.data[i])
    st.persistence()
    dgm_dim0 = st.persistence_intervals_in_dimension(0)
    
    D = sparse.csr_matrix(D)
    if max_filtration is None:
        max_filtration = nd_weight_max - 1e-10

    ph_result = ripser(D, distance_matrix=True, maxdim=max_dim, thresh=max_filtration)
    dgms = ph_result['dgms']

    dgms[0] = dgm_dim0
    if len(dgms) == 2:
        diagrams[cell_names[icell]] = {'max_filtration': max_filtration, 
            'dgm_h0': dgms[0], 'dgm_h1': dgms[1]}

def ph_process_lwph(
    A, # adjacency graph of cells
    D, # distance matrix of the cells from the embedding
    icell,
    diagrams,
    nb_method, # 'distance or k nearest neighbors
    nb_distance,
    nb_knn,
    filtration_value, # 'euclidean' or 'graph'
    max_dim,
    max_filtration,
    cell_names,
):
    """
    Local weighted persistent homology
    """
    if nb_method == 'distance':
        idx = np.where(D[icell,:] <= nb_distance)[0]
    elif nb_method == 'knn':
        idx = np.argsort(D[icell,:])[:nb_knn]
    if filtration_value == 'euclidean':
        D_sub = D[idx,:][:,idx]
    
    ph_result = ripser(D_sub, distance_matrix=True, maxdim=max_dim, thresh=max_filtration)
    dgms = ph_result['dgms']

    if len(dgms) == 1:
        diagrams[cell_names[icell]] = {'max_filtration': max_filtration, 
            'dgm_h0': dgms[0]}
    if len(dgms) == 2:
        diagrams[cell_names[icell]] = {'max_filtration': max_filtration, 
            'dgm_h0': dgms[0], 'dgm_h1': dgms[1]}
    if len(dgms) == 3:
        diagrams[cell_names[icell]] = {'max_filtration': max_filtration, 
            'dgm_h0': dgms[0], 'dgm_h1': dgms[1], 'dgm_h2': dgms[2]}

def ph_process_rlph(
    X,
    f,
    icell,
    diagrams,
    metric,
    relative_method,
    rlph_base,
    nb_knn,
    rl_knn,
    nb_distance,
    rl_distance,
    max_dim,
    max_filtration,
    cell_names,
):
    """
    Relative local persistent homology
    """
    if rlph_base == 'nb_distance':
        tmp_dists = distance_matrix(X, X[icell,:].reshape(1,X.shape[1])).flatten()
        idx = list( np.where(tmp_dists <= nb_distance)[0] )
    elif rlph_base == 'nb_knn':
        tmp_dists = distance_matrix(X, X[icell,:].reshape(1,X.shape[1])).flatten()
        idx = list( np.argsort(tmp_dists)[:nb_knn+1] )
    if rlph_base in ['nb_distance', 'nb_knn']:
        idx.remove(icell)
        idx = [icell] + idx
        idx = np.array(idx, int)
        XX = X[idx,:].copy()
        dists = pdist(XX, metric)
        f = dys.fill_rips(dists, max_dim+1, max_filtration)
    elif rlph_base == 'all':
        XX = X.copy()
    
    if rlph_base == 'all':
        tmp_icell = icell
    elif rlph_base in ['nb_distance','nb_knn']:
        tmp_icell = 0

    if relative_method == 'self':
        s = dys.Simplex([tmp_icell])
        f1 = dys.Filtration([s])
    elif relative_method == 'knn':
        tmp_dists = distance_matrix(XX, XX[tmp_icell,:].reshape(1,X.shape[1])).flatten()
        neighbors = set(list(np.argsort(tmp_dists)[:rl_knn+1]))
        f1 = dys.Filtration([s for s in f if set(list(s)).issubset(neighbors)])
    elif relative_method == 'distance':
        tmp_dists = distance_matrix(XX, XX[tmp_icell,:].reshape(1,X.shape[1])).flatten()
        neighbors = set(list(np.where(tmp_dists <= rl_distance)[0]))
        f1 = dys.Filtration([s for s in f if set(list(s)).issubset(neighbors)])

    m = dys.homology_persistence(f, relative = f1)
    tmp_dgm = dys.init_diagrams(m, f)
    dgms_rel = []
    for i in range(max_dim+1):
        dgms_rel.append(np.array( [[pt.birth, pt.death] for pt in tmp_dgm[i]] ))

    if rlph_base in ['nb_distance', 'nb_knn']:
        m = dys.homology_persistence(f)
        tmp_dgm = dys.init_diagrams(m,f)
        dgms_base = []
        for i in range(max_dim+1):
            dgms_base.append(np.array( [[pt.birth, pt.death] for pt in tmp_dgm[i]] ))

    if rlph_base == 'all':     
        if len(dgms_rel) == 1:
            diagrams[cell_names[icell]] = {'max_filtration': max_filtration, 
                'dgm_h0': dgms_rel[0]}
        if len(dgms_rel) == 2:
            diagrams[cell_names[icell]] = {'max_filtration': max_filtration, 
                'dgm_h0': dgms_rel[0], 'dgm_h1': dgms_rel[1]}
        if len(dgms_rel) == 3:
            diagrams[cell_names[icell]] = {'max_filtration': max_filtration, 
                'dgm_h0': dgms_rel[0], 'dgm_h1': dgms_rel[1], 'dgm_h2': dgms_rel[2]}
    elif rlph_base in ['nb_distance', 'nb_knn']:
        if len(dgms_rel) == 1:
            diagrams[cell_names[icell]] = {'max_filtration': max_filtration, 
                'dgm_h0': dgms_rel[0], 
                'dgm_h0_base': dgms_base[0]}
        if len(dgms_rel) == 2:
            diagrams[cell_names[icell]] = {'max_filtration': max_filtration, 
                'dgm_h0': dgms_rel[0], 'dgm_h1': dgms_rel[1],
                'dgm_h0_base': dgms_base[0], 'dgm_h1_base': dgms_base[1]}
        if len(dgms_rel) == 3:
            diagrams[cell_names[icell]] = {'max_filtration': max_filtration, 
                'dgm_h0': dgms_rel[0], 'dgm_h1': dgms_rel[1], 'dgm_h2': dgms_rel[2],
                'dgm_h0_base': dgms_base[0], 'dgm_h1_base': dgms_base[1], 'dgm_h2_base': dgms_base[2]}
    
def ph_process_ewvr(
    A, # edge weight larger weight means higher correlation
    icell,
    diagrams,
    max_filtration,
    max_dim,
    cell_names,
):
    """
    edge weighted VR filtration
    """
    edge_weight_max = A.max()
    A_filtration = A.copy()
    A_filtration.data = edge_weight_max - A_filtration.data
    if max_filtration is None:
        max_filtration = edge_weight_max - 1e-10
    ph_result = ripser(A_filtration, distance_matrix=True, maxdim=max_dim, thresh=max_filtration)
    dgms = ph_result['dgms']

    if len(dgms) == 2:
        diagrams[cell_names[icell]] = {'max_filtration': max_filtration,
            'dgm_h0': dgms[0], 'dgm_h1': dgms[1]}

def ph_feature_persimages(adata,
    inf_value = 'replace_max',
    pixel_sizes = [0.5,0.5],
    dims = [0,1],
    diagram_name = 'diagrams',
    birth_ranges = None,
    pers_ranges = None,
):
    ncell = adata.shape[0]
    diagrams = adata.uns[diagram_name]
    cell_names = list(adata.obs_names)
    images = []
    if birth_ranges is None:
        birth_ranges = [None for i in range(len(dims))]
    if pers_ranges is None:
        pers_ranges = [None for i in range(len(dims))]
    used_birth_ranges = []
    used_pers_ranges = []
    for idim in range(len(dims)):
        dim = dims[idim]
        dgm_list = []
        for i in range(ncell):
            cell_name = cell_names[i]
            max_filtration = diagrams[cell_name]['max_filtration']
            dgm = diagrams[cell_name]['dgm_h%d' % dim].copy()
            if inf_value == 'replace_max':
                dgm[np.where(dgm==np.inf)] = max_filtration
            dgm_list.append(dgm)
        
        pers_imager = pimg.PersistenceImager(birth_range=birth_ranges[idim], 
            pers_range=pers_ranges[idim], pixel_size=pixel_sizes[idim])
        pers_imager.weight = pimg.weighting_fxns.persistence
        pers_imager.weight_params = {'n': 2}
        if dim == 0:
            for i in range(ncell):
                dgm_list[i][:,0] = np.random.rand(dgm_list[i].shape[0])*pixel_sizes[idim] * 1E-5
        pers_imager.fit(dgm_list, skew=True)        
        if not birth_ranges[idim] is None:
            pers_imager.birth_range = birth_ranges[idim]
        if not pers_ranges[idim] is None:
            pers_imager.pers_range = pers_ranges[idim]
        used_birth_ranges.append(pers_imager.birth_range)
        used_pers_ranges.append(pers_imager.pers_range)
        for i in range(ncell):
            img = pers_imager.transform(dgm_list[i])
            if i == 0:
                image_arr = np.empty([ncell]+list(img.shape), float)
            image_arr[i,:,:] = img[:,:]
        images.append(image_arr)
    
    persimage_info = {"dims":dims, 
        "birth_ranges":used_birth_ranges, 
        "pers_ranges":used_pers_ranges, 
        "pixel_sizes":pixel_sizes}
    return images, persimage_info

def ph_feature_total_persistence(
    adata,
    diagram_name = None,
    dims = [0,1],
    inf_value = 'replace_max',
):
    diagrams = adata.uns[diagram_name]
    cell_names = list( adata.obs_names )
    ncell = len(cell_names)
    total_persistence = np.zeros([ncell,len(dims)], float)
    for idim in range(len(dims)):
        dim = dims[idim]
        for icell in range(ncell):
            cell_name = cell_names[icell]
            max_filtration = diagrams[cell_name]['max_filtration']
            dgm = diagrams[cell_name]['dgm_h%d' % dim].copy()
            if inf_value == 'replace_max':
                dgm[np.where(dgm==np.inf)] = max_filtration
            total_persistence[icell,idim] = np.sum(np.abs(dgm[:,1]-dgm[:,0]))
    return total_persistence

def ph_feature_persistence_std(
    adata,
    diagram_name = None,
    dims = [0,1],
    inf_value = 'replace_max',
):
    diagrams = adata.uns[diagram_name]
    cell_names = list( adata.obs_names )
    ncell = len(cell_names)
    persistence_std = np.zeros([ncell,len(dims)], float)
    for idim in range(len(dims)):
        dim = dims[idim]
        for icell in range(ncell):
            cell_name = cell_names[icell]
            max_filtration = diagrams[cell_name]['max_filtration']
            dgm = diagrams[cell_name]['dgm_h%d' % dim].copy()
            if inf_value == 'replace_max':
                dgm[np.where(dgm==np.inf)] = max_filtration
            persistence_std[icell,idim] = np.sum(np.abs(dgm[:,1]-dgm[:,0]))
    return persistence_std
        
def ph_feature_persistence_entropy(
    adata,
    diagram_name = None,
    dims = [0,1],
    inf_value = 'replace_max',
):
    diagrams = adata.uns[diagram_name]
    cell_names = list( adata.obs_names )
    ncell = len(cell_names)
    X_pe = np.empty([ncell,len(dims)], float)
    for idim in range(len(dims)):
        dim = dims[idim]
        dgm_list = []
        for i in range(ncell):
            cell_name = cell_names[i]
            max_filtration = diagrams[cell_name]['max_filtration']
            dgm = diagrams[cell_name]['dgm_h%d' % dim].copy()
            if inf_value == 'replace_max':
                dgm[np.where(dgm==np.inf)] = max_filtration+1e-5
            dgm_list.append(dgm)
        PE = gudhi.representations.Entropy()
        pe = PE.fit_transform(dgm_list)
        X_pe[:,idim] = pe[:,0]
    return X_pe

def ph_feature_betti_curve(
    adata,
    diagram_name = None,
    dims = [0,1],
    inf_value = 'replace_max',
    n_intervals = [20,20],
    v_mins = [None, None],
    v_maxs = [None, None],
):
    diagrams = adata.uns[diagram_name]
    cell_names = list( adata.obs_names )
    ncell = len(cell_names)
    X_bc = [np.zeros([ncell,n_intervals[i]]) for i in range(len(dims)) ]
    xs_list = []
    for idim in range(len(dims)):
        dim = dims[idim]
        dgm_list = []
        for i in range(ncell):
            cell_name = cell_names[i]
            max_filtration = diagrams[cell_name]['max_filtration']
            dgm = diagrams[cell_name]['dgm_h%d' % dim].copy()
            if inf_value == 'replace_max':
                dgm[np.where(dgm==np.inf)] = max_filtration
            dgm_list.append(dgm)
        if v_mins[idim] is None:
            tmp_vals = [np.min(dgm_list[i][:,0]) for i in range(ncell)]
            v_min = np.min(tmp_vals)
        else:
            v_min = v_mins[idim]
        if v_maxs[idim] is None:
            tmp_vals = [np.max(dgm_list[i][:,1]) for i in range(ncell)]
            v_max = np.max(tmp_vals)
        else:
            v_max = v_maxs[idim]
        xs = np.linspace(v_min, v_max, n_intervals[idim])
        xs_list.append(xs)
        bc = gudhi.representations.BettiCurve(predefined_grid=xs)
        bettis = bc.fit_transform(dgm_list)
        X_bc[idim][:,:] = bettis[:,:]
    return X_bc, xs_list

def ph_prune_diagrams(
    adata,
    diagram_name = None,
    dims = None,
    min_per = None,
    inf_value = None,
):
    ncell = adata.shape[0]
    cell_names = list( adata.obs_names )
    for idim in range(len(dims)):
        dim = dims[idim]
        for icell in range(ncell):
            cell_name = cell_names[icell]
            tmp_dgm = adata.uns[diagram_name][cell_name]['dgm_h%d' % dim].copy()
            per = tmp_dgm[:,1] - tmp_dgm[:,0]
            keep_idx = np.where(per >= min_per[idim])[0]
            adata.uns[diagram_name][cell_name]['dgm_h%d' % dim] = tmp_dgm[keep_idx,:]
            del tmp_dgm
    return