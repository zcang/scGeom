import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import pkgutil
import io

def get_gene_network(
    adata: anndata.AnnData,
    database: str = 'scent_17',
    species: str = 'human',
):
    """
    Get the gene network from knowledge database from SCENT [Banerji2013].

    Parameters
    ----------
    adata
        The ``anndata`` object of gene expression data.
    database
        Choose between 'scent_17' and 'scent_13' for the SCENT network updated in 2017 or 2013 respectively.
    species
        'human' or 'mouse'
    
    Returns
    -------
    gene_network: dict
        A gene network stored in ``adata.uns[gene_network]``. The adjacency matrix of the network is in ``adata.uns[gene_network]['network']``.
        The genes are in ``adata.uns[gene_network]['genes']``.

    References
    ----------
    [Banerji2013] Banerji, C. R., Miranda-Saavedra, D., Severini, S., Widschwendter, M., Enver, T., Zhou, J. X., & Teschendorff, A. E. (2013). Cellular network entropy as the energy potential in Waddington's differentiation landscape. Scientific reports, 3(1), 3039.
    """
    if database == 'scent_17':
        data = pkgutil.get_data(__name__, "_data/NetworkDatabase/SCENT_network/net17jan16_filtered_"+species+".npz")
        A = sparse.load_npz(io.BytesIO(data))
        data = pkgutil.get_data(__name__, "_data/NetworkDatabase/SCENT_network/net17jan16_filtered_genes_"+species+".csv")
        genes = list( pd.read_csv(io.BytesIO(data), header=None).values.flatten() )
    elif database == 'scent_13':
        data = pkgutil.get_data(__name__, "_data/NetworkDatabase/SCENT_network/net13jun12_filtered_"+species+".npz")
        A = sparse.load_npz(io.BytesIO(data))
        data = pkgutil.get_data(__name__, "_data/NetworkDatabase/SCENT_network/net13jun12_filtered_genes_"+species+".csv")
        genes = list( pd.read_csv(io.BytesIO(data), header=None).values.flatten() )
    
    adata_genes = list( adata.var_names )
    common_genes = list( set(adata_genes).intersection(set(genes)) )
    common_index = [genes.index(gene) for gene in common_genes]
    A_common = A[common_index,:][:,common_index]
    gene_network = {'database': database,
        'genes': common_genes,
        'network': A_common}
    adata.uns['gene_network'] = gene_network

    return