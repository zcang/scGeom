import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy import sparse
import pkgutil
import io
import gc
import sys
import locCSN

def compute_pairwise_scc(X1, X2):
    X1 = X1.argsort(axis=1).argsort(axis=1)
    X2 = X2.argsort(axis=1).argsort(axis=1)
    X1 = (X1-X1.mean(axis=1, keepdims=True))/X1.std(axis=1, keepdims=True)
    X2 = (X2-X2.mean(axis=1, keepdims=True))/X2.std(axis=1, keepdims=True)
    sccmat = np.empty([X1.shape[0], X2.shape[0]], float)
    for i in range(X1.shape[0]):
        for j in range(X2.shape[0]):
            c = np.dot( X1[i,:], X2[j,:]) / float(X1.shape[1])
            sccmat[i,j] = c
    return sccmat

import numpy as np
from scipy.stats import norm
from scipy.sparse import csr_matrix
np.set_printoptions(formatter={'float': lambda x: "{0:0.4f}".format(x)})
def round_half_up(x):
    if x == 0.5:
        return 1
    else:
        return int(round(x))

def csnet(data, c=None, alpha=0.01, boxsize=0.1, weighted=False):
    """
    Construction of cell-specific networks.

    Args:
    data: Gene expression matrix, rows = genes, columns = cells
    c: Construct the CSNs for all cells, set c = None (Default);
        Construct the CSN for cell k, set c = k
    alpha: Significant level (e.g., 0.001, 0.01, 0.05 ...)
           larger alpha leads to more edges, Default = 0.01
    boxsize: Size of neighborhood, Default = 0.1
    weighted: 1  edge is weighted
               0  edge is not weighted (Default)

    Returns:
    csn: Cell-specific network, the kth CSN is in csn[k]
         rows = genes, columns = genes
    """
    if weighted is None:
        weighted = False
    if boxsize is None:
        boxsize = 0.1
    if alpha is None:
        alpha = 0.01

    n1, n2 = data.shape
    if c is None:
        c = np.arange(n2)

    # Define the neighborhood of each plot
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

    # Construction of cell-specific network
    csn = {}
    p = -norm.ppf(alpha, 0, 1)
    for k in c:
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
        if k == 0:
            d_0 = d
            

        if k % 100 == 0:
            print(f"Cell {k} is completed")

    return csn

def infer_gene_network(
    adata: anndata.AnnData,
    method: str = 'spearman',
    threshold: float= 0.2,
    network_name: str = 'gene_network',
    return_network: bool = False,
):
    """
    Simple computation of gene networks with Spearman's correlation coefficient.

    Parameters
    ----------
    adata
        The data matrix of shape ``n_obs`` × ``n_var``.
    method
        Method to use, currently, only Spearman's r is implemented.
    threshold
        The threshold of absolute value of correlation coefficients to keep an edge.
    network_name
        The name of the network to be stored in ``adata.uns[network_name]``.
    return_network
        Whether to return the network.

    Returns
    -------
    gene_network : dict
        A pandas DataFrame of the LR pairs with the three columns representing the ligand, receptor, and the signaling pathway name, respectively.
    """
    
    X = np.array( adata.X )
    ngene = X.shape[1]
    if method == 'spearman':
        sccmat = compute_pairwise_scc(X.T, X.T)
        network = np.zeros([ngene, ngene], int)
        network[np.where(np.abs(sccmat)>=threshold)] = 1
        gene_network = {'method': method,
            'genes': list(adata.var_names),
            'network': sparse.csr_matrix(network),
            'sccmat': sccmat}
    adata.uns[network_name] = gene_network
    if return_network:
        return gene_network
    else:
        return

def cell_specific_gene_network(
    adata: anndata.AnnData,
    method: str = 'csn',
    csn_alpha: float = 0.01,
    csn_boxsize: float = 0.1,
    csn_weighted: bool = True,
    network_name: str = 'cell_specific_gene_network',
    return_network: bool = False,
):
    """
    Compute a cell-specific gene network for each cell.

    Parameters
    ----------
    adata
        The AnnData object of gene expression data.
    method
        The method to used. 'csn' for cell-specific network [Dai2019]_ or 'loccsn' for local cell-specific network [Wang2021]_.
    csn_alpha
        When method is 'csn', the alpha value (significance for keeping an edge).
    csn_boxsize
        When method is 'csn', the size of box to define neiborhood of gene expression.
    csn_weighted
        When method is 'csn', whether use weighted output or binary output.
    network_name
        The name of the dictionary of networks to be stored in ``adata.uns[network_name]``. Use the ``obs_names`` of cells to access the network for each cell.
    return_network
        Whether to return a copy of the dictionary of networks.
    
    Returns
    -------
    gene_networks: dict
        A dictionary of cell-specific networks with cell names in ``adata.obs_names`` as dictionary keys.

    References
    ----------

    .. [Dai2019] Dai, H., Li, L., Zeng, T., & Chen, L. (2019). Cell-specific network constructed by single-cell RNA sequencing data. Nucleic acids research, 47(11), e62-e62.
    .. [Wang2021] Wang, X., Choi, D., & Roeder, K. (2021). Constructing local cell-specific networks from single-cell data. Proceedings of the National Academy of Sciences, 118(51), e2113178118.

    """
    data = np.array(adata.X).T
    if method == 'csn':
        csn = csnet(data, alpha=csn_alpha, boxsize=csn_boxsize, weighted=csn_weighted)
    elif method == 'loccsn':
        csn = locCSN.csn(data, dev = True)
    cell_names = list(adata.obs_names)
    gene_networks = {}
    for i in range(len(cell_names)):
        gene_networks[cell_names[i]] = csn[i]
    
    if return_network:
        return gene_networks
    else:
        adata.uns[network_name] = gene_networks
        return