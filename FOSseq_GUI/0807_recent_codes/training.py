from pathlib import Path

import numpy as np
import anndata
import pandas as pd
from scipy import stats
from scipy.sparse import vstack

from sklearn.preprocessing import normalize
from sklearn.neighbors import NearestNeighbors

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache


def list_file_name(option1):
    files = {
        "WMB-10Xv2": ['WMB-10Xv2-CTXsp', 'WMB-10Xv2-HPF', 'WMB-10Xv2-HY', 'WMB-10Xv2-Isocortex-1', 
                    'WMB-10Xv2-Isocortex-2', 'WMB-10Xv2-Isocortex-3', 'WMB-10Xv2-Isocortex-4',
                    'WMB-10Xv2-MB', 'WMB-10Xv2-OLF', 'WMB-10Xv2-TH'],
        "WMB-10Xv3": ['WMB-10Xv3-CB', 'WMB-10Xv3-CTXsp', 'WMB-10Xv3-HPF', 'WMB-10Xv3-HY',
                    'WMB-10Xv3-Isocortex-1', 'WMB-10Xv3-Isocortex-2', 'WMB-10Xv3-MB', 'WMB-10Xv3-MY',
                    'WMB-10Xv3-OLF', 'WMB-10Xv3-P', 'WMB-10Xv3-PAL', 'WMB-10Xv3-STR', 'WMB-10Xv3-TH'],
        "WMB-10XMulti": ['WMB-10XMulti']
    }
    result = []
    for opt in option1:
        result.extend(files.get(opt, []))
    return result

# ===============================
# HELPER FUNCTIONS
# ===============================
# By Option 1, 2, 3
def get_feature_matrix_label(o2, o3=0):
    if o3 == 0:
        feature_data_file = o2 + '/log2'
    else:
        feature_data_file = o2 + '/raw'

    return feature_data_file


# ===============================
# MAIN JOB
# ===============================
def get_mer_gene_data(abc_cache):
    file = abc_cache.get_data_path(
        directory='MERFISH-C57BL6J-638850',
        file_name='C57BL6J-638850/log2'
    )
    adata = anndata.read_h5ad(file, backed='r')

    return adata


def get_10x_adata_list(abc_cache, o1, o2, o3):
    adata_list = []
    cell_number = 0
    for dir in o1:
        for label in o2:
            if label.startswith(dir):
                feature_data_file = get_feature_matrix_label(label, o3)
                file = abc_cache.get_data_path(directory = dir, file_name = feature_data_file)
                adata = anndata.read_h5ad(file, backed = 'r')
                adata_list.append(adata)

                cell_number += len(adata.obs)
    
    return adata_list, cell_number


def get_common_genes(mdata, adata):
    common_genes = mdata.var_names.intersection(adata.var_names)
    return common_genes


def get_ref_df(genes, adata_list):
    sparse_matrices = []
    obs_names = []
    adata_indices = []
    for index, adata in adata_list:
        adata = adata[:, genes]
        sparse_matrices.append(adata.X)
        obs_names.extend(adata.obs_names)
        adata_indices.extend([index] * adata.n_obs)
    
    combined_sparse_matrix = vstack(sparse_matrices)
    combined_df = pd.DataFrame(combined_sparse_matrix.toarray(), index=obs_names, columns=genes)
    index_df = pd.DataFrame(
        {'adata_index': adata_indices},
        index=combined_df.index
    )

    return combined_df, index_df
    # filtered_dfs = []
    # for adata in adata_list:
    #     adata = adata[:, genes]
    #     adata_df = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
    #     filtered_dfs.append(adata_df)

    # combined_df = pd.concat(filtered_dfs)
    # return combined_df