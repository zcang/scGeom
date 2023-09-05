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

def ph_process_wvr(
    A,
    X,
    icell,
    diagrams,
    max_filtration,
    max_dim,
):
    return

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








def value_based_persistent_homology(
    G,
    f,
    dtm_p = 2,
    thresh = np.inf,
    plot = True,
    filtration = 'dtm'
):
    n = f.shape[0]
    st = gudhi.SimplexTree()
    for i in range(n):
        if filtration == 'dtm':
            value = f[i]
        elif filtration == 'tomato':
            value = f[i]
        if value <= thresh:
            st.insert([i], filtration=value)
    I, J, A = [], [], []
    for (u, v, d) in G.edges(data=True):
        if filtration == 'dtm':
            value = filtration_value(dtm_p, f[u], f[v], d['weight'])
        elif filtration == 'tomato':
            value = np.mean([f[u], f[v]])
        if value  < thresh:
            st.insert([u,v], filtration=value)
            I.extend([u,v]); J.extend([v,u]); A.extend([value, value])
    I = np.array(I, int)
    J = np.array(J, int)
    A = np.array(A, float)
    D_rips = sparse.coo_matrix((A, (I, J)), shape=(n, n)).tocsr()         
    st.expansion(1)
    diagram = st.persistence()
    if plot:
        gudhi.plot_persistence_diagram(diagram)
        plt.show()
    return diagram, D_rips


def persistent_homology_clustering(
    G,
    f,
    persistence
):
    st = gudhi.SimplexTree()
    G_srt = nx.Graph()
    G_srt.add_nodes_from([i for i in range(len(G))])
    
    


def persistent_homology_dtm(
    D,
    f,
    p=2,
    thresh=np.inf,
    plot=True,
    k = None,
):
    n = f.shape[0]
    st = gudhi.SimplexTree()
    I,J,A = [],[],[]
    for i in range(n):
        value = f[i]
        if value < thresh:
            st.insert([i], filtration= value)
    if k is None:
        for i in range(n-1):
            for j in range(i+1, n):
                value = filtration_value(p, f[i], f[j], D[i][j])
                if value < thresh:
                    st.insert([i,j], filtration = value)
                    I.extend([i,j])
                    J.extend([j,i])
                    A.extend([value,value])
    else:
        for i in range(n):
            values = np.empty([n], float)
            for j in range(n):
                value = filtration_value(p, f[i], f[j], D[i][j])
                values[j] = value
            idx_srt = np.argsort(values)
            for j in range(1, k+1):
                if values[idx_srt[j]] < thresh:
                    st.insert([i, idx_srt[j]], filtration=values[idx_srt[j]])
                    I.extend([i,idx_srt[j]])
                    J.extend([idx_srt[j],i])
                    A.extend([values[idx_srt[j]], values[idx_srt[j]]])
        I = np.array(I, int)
        J = np.array(J, int)
        A = np.array(A, float)

    D_rips = sparse.coo_matrix((A, (I, J)), shape=(n, n)).tocsr()
    # ripser_diagram = ripser(D_rips, distance_matrix=True)['dgms']
    st.expansion(1)
    diagram = st.persistence()
    # for i in range(ripser_diagram[1].shape[0]):
    #     birth = ripser_diagram[1][i,0]
    #     death = ripser_diagram[1][i,1]
    #     diagram.append((1,(birth, death)))
    if plot:
        gudhi.plot_persistence_diagram(diagram)
        plt.show()
    return diagram, D_rips

def distance_to_measure_persistent_homology(D, f, m=0.05, p=2, thresh=np.inf, k=None, use_kdtree=False, plot=True):
    n = f.shape[0]
    # n = X.shape[0]
    # f = distance_to_measure(X, m=m, p=p)
    if not use_kdtree:
        # D = distance_matrix(X, X)
        st = gudhi.SimplexTree()
        I,J,A = [],[],[]
        for i in range(n):
            value = f[i]
            if value < thresh:
                st.insert([i], filtration = value)
        
        if k is None:
            for i in range(n-1):
                for j in range(i+1, n):
                    value = Filtration_value(p, f[i], f[j], D[i][j])
                    if value <= thresh:
                        st.insert([i,j], filtration=value)
                        I.extend([i,j])
                        J.extend([j,i])
                        A.extend([value, value])
        else:
            for i in range(n):
                values = np.empty([n], float)
                for j in range(n):
                    value = Filtration_value(p, f[i], f[j], D[i][j])
                    values[j] = value
                ind_srt = np.argsort(values)
                for j in range(1,k+1):
                    if values[ind_srt[j]] <= thresh:
                        st.insert([i,ind_srt[j]], filtration=values[ind_srt[j]])
                        I.extend([i,ind_srt[j]])
                        J.extend([ind_srt[j],i])
                        A.extend([values[ind_srt[j]], values[ind_srt[j]]])
        I = np.array(I, int)
        J = np.array(J, int)
        A = np.array(A, float)
    
    D_rips = sparse.coo_matrix((A, (I, J)), shape=(n, n)).tocsr()
    ripser_diagram = ripser(D_rips, distance_matrix=True)['dgms']
    st.expansion(1)
    diagram = st.persistence()
    for i in range(ripser_diagram[1].shape[0]):
        birth = ripser_diagram[1][i,0]
        death = ripser_diagram[1][i,1]
        diagram.append((1,(birth, death)))
    if plot:
        gudhi.plot_persistence_diagram(diagram)
        plt.show()
    return diagram, D_rips


def filtration_value(p, fx, fy, d, n = 10):
    '''
    Compute the filtrations values of the edge [x,y] in the weighted Rips filtration
    If p is not 1, 2 or 'np.inf, an implicit equation is solved
    The equation to solve is G(I) = d, where G(I) = (I**p-fx**p)**(1/p)+(I**p-fy**p)**(1/p)
    We use a dichotomic method
    
    Input:
    p: parameter of the weighted Rips filtration, in [1, +inf) or np.inf
    fx: filtration value of the point x
    fy: filtration value of the point y
    d: distance between the points x and y
    n: number of iterations of the dichotomic method
        
    Output: 
    val : filtration value of the edge [x,y], i.e. solution of G(I) = d    
    
    Example:
    Filtration_value(2.4, 2, 3, 5, 10)
    '''
    if p==np.inf:
        value = max([fx,fy,d/2])
    else:
        fmax = max([fx,fy])
        if d < (abs(fx**p-fy**p))**(1/p):
            value = fmax
        elif p==1:
            value = (fx+fy+d)/2
        elif p==2:
            if d == 0.0:
                value = 0.0
            else:
                value = np.sqrt( ( (fx+fy)**2 +d**2 )*( (fx-fy)**2 +d**2 ) )/(2*d)            
        else:
            Imin = fmax; Imax = (d**p+fmax**p)**(1/p)
            for i in range(n):
                I = (Imin+Imax)/2
                g = (I**p-fx**p)**(1/p)+(I**p-fy**p)**(1/p)
                if g<d:
                    Imin=I
                else:
                    Imax=I
            value = I
    return value

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

class tomato(object):
    def __init__(self, G=None, f=None):
        self.G = G
        self.f = f
        self.n = nx.number_of_nodes(G)
        self.g = np.empty([self.n],float)
        self.r = np.empty([self.n],float)

    
    def persistent_homology(self, plot=True):

        self.sxt = gudhi.SimplexTree()

        for ind in range(self.n):
            self.sxt.insert([ind], filtration=-self.f[ind])
            nei = list(self.G.neighbors(ind))
            for idx in nei:
                self.sxt.insert([ind, idx], filtration=np.max([-self.f[ind], -self.f[idx]]))
                
        self.sxt.initialize_filtration()
        self.sxt.persistence()
        dig, res = self.sxt.persistence(), []
        bars = []
        for ele in dig:
            if ele[0] == 0:
                birth = -ele[1][0]
                death = -ele[1][1]
                bar = (0, (birth, death))
                bars.append([birth, death])
                res.append(bar)
        # gudhi.plot_persistence_barcode(res)
        # gudhi.plot_persistence_diagram(res)
        bars = np.asarray(bars)
        self.bars = bars
        if plot:
            for i in range(bars.shape[0]):
                plt.plot([bars[i,0]], [max(0,bars[i,1])], "bo")
            plt.plot([0,np.max(bars[:,0])], [0,np.max(bars[:,0])], 'k')
            plt.show()

    def determine_tau(self, ncl):
        persis = np.sort(np.abs(self.bars[:,0]-self.bars[:,1]))
        # print(persis)
        if ncl >= len(persis):
            tau = 0.0
        else:
            tau = 0.5 * (persis[-ncl-1] + persis[-ncl])
        return tau

    def clustering(self, tau):

        G_srt = nx.Graph()
        G_srt.add_nodes_from([i for i in range(len(self.G))])
        idmap_srt2org = np.argsort(-self.f)
        idmap_org2srt = {}
        for i in range(len(idmap_srt2org)):
            idmap_org2srt[idmap_srt2org[i]] = i
        org_edges = list(self.G.edges())
        for org_edge in org_edges:
            G_srt.add_edge(idmap_org2srt[org_edge[0]], idmap_org2srt[org_edge[1]])
        f_srt = -np.sort(-self.f)

        U = UnionFind()
        r = np.empty([self.n], int)
        g = np.empty([self.n], int)
        e = -1
        for i in range(self.n):
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
        label_srt = np.empty([self.n], int)
        for i in range(self.n):
            label_srt[i] = U.find(i)
        label = np.empty_like(label_srt)
        for i in range(self.n):
            label[i] = idmap_srt2org[label_srt[idmap_org2srt[i]]]

        return label
