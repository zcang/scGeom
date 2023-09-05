import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import pkgutil
import io

def get_gene_network(
    adata,
    database = 'scent_17',
    species = 'human',
):
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
    common_genes = set(adata_genes).intersection(set(genes))
    common_index = [genes.index(gene) for gene in common_genes]
    A_common = A[common_index,:][:,common_index]
    gene_network = {'database': database,
        'genes': common_genes,
        'network': A_common}
    adata.uns['gene_network'] = gene_network

    return