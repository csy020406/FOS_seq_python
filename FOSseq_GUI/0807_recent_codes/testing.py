from pathlib import Path
import json

import pandas as pd
import numpy as np
from scipy import stats

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache


def get_mer_metadata(abc_cache):
    cell = abc_cache.get_metadata_dataframe(
        directory='MERFISH-C57BL6J-638850',
        file_name='cell_metadata_with_cluster_annotation'
    )
    cell.set_index('cell_label', inplace=True)
    cell = cell.loc[:, ['brain_section_label', 'cluster_alias']]

    # CCF coordinates
    coords = abc_cache.get_metadata_dataframe(
        directory='MERFISH-C57BL6J-638850-CCF',
        file_name='ccf_coordinates',
        dtype={"cell_label": str}
    )
    coords.set_index('cell_label', inplace=True)
    cell = cell.join(coords, how='inner')

    # Parcellation annotation
    parcellation_annotation = abc_cache.get_metadata_dataframe(
        directory='Allen-CCF-2020',
        file_name='parcellation_to_parcellation_term_membership_acronym'
    )                                                            
    parcellation_annotation.set_index('parcellation_index', inplace=True)
    parcellation_annotation.columns = ['parcellation_%s'% x for x in  parcellation_annotation.columns]
    parcellation_annotation = parcellation_annotation.loc[:, ['parcellation_substructure']]
    cell = cell.join(parcellation_annotation, on='parcellation_index')

    return cell



def get_mer_dat(mer_metadata, mer_adata, genes, roi):
    pred = (mer_metadata['parcellation_substructure'].isin(roi))

    # Use copy() to avoid SettingWithCopyWarning error
    filtered_meta = mer_metadata[pred].copy()

    asubset = mer_adata[filtered_meta.index, genes].to_memory
    mer_dat_df = pd.DataFrame(asubset.X.toarray(), index=asubset.obs_names, columns=asubset.var_names)

    return filtered_meta, mer_dat_df

# NOTE: 임시 함수
def cal_weighted_average(indices, distances, adata_df):
    weighted_averages = []
    for i in range(len(indices)):
        cell_indices = indices[i]
        cell_distances = distances[i]

        # weights = cell_distances - np.max(cell_distances) - 0.001
        # weights = -weights / np.max(np.abs(weights))
        # weights = weights / np.sum(weights)
        inverse_distances = 1 / cell_distances
        sum_inverse_distances = np.sum(inverse_distances)
        weights = inverse_distances / sum_inverse_distances

        weighted_avg = np.dot(weights, adata_df.iloc[cell_indices])
        weighted_averages.append(weighted_avg)

    return weighted_averages

# TODO: 구현
# def cal_weighted_average(ref_idx_df, indices, distances, adata_list, start, end):
#     num_cells = indices.shape[0]
#     num_genes = end - start
#     num_neighbors = 15
#     batch_size = 200

#     weights = distances #TODO: 구현

#     unique_indices = np.unique(indices.flatten())
#     unique_cell_labels = ref_idx_df.index[unique_indices]
#     cell_to_adata_index = {label: ref_idx_df.loc[label, 'adata_index'] for label in unique_cell_labels}

#     weighted_avg = np.zeros((num_cells, num_genes))

#     for adata_index, adata in enumerate(adata_list):
#         ############
#         relevant_labels = [label for label, idx in cell_to_adata_index.items() if idx == adata_index]
#         if not relevant_labels:
#             continue
#         obs_names = adata.obs_names.isin(relevant_labels)

#         # cell_label_to_idx = {label: adata.obs_names.get_loc(label) for label in relevant_labels}
        
#         # Load data in batches
#         for batch_start in range(start, end, batch_size):
#             batch_end = min(batch_start + batch_size, end)
#             gene_subset = adata[obs_names, batch_start:batch_end].to_memory()
            
#             # Iterate over each target cell
#             for target_idx in range(num_cells):
#                 # Get the neighbor indices and weights
#                 neighbor_indices = indices[target_idx]
#                 neighbor_weights = weights[target_idx]
                
#                 # Accumulate weighted sums for target cell
#                 weighted_sum = np.zeros(batch_end - batch_start)
#                 for neighbor_idx, weight in zip(neighbor_indices, neighbor_weights):
#                     neighbor_cell_label = ref_idx_df.index[neighbor_idx]
#                     neighbor_adata_index = cell_to_adata_index[neighbor_cell_label]
                    
#                     if neighbor_adata_index == adata_index:
#                         neighbor_cell_index = cell_label_to_idx[neighbor_cell_label]
#                         weighted_sum += weight * gene_subset[neighbor_cell_index, :]
                
#                 # Store the weighted average for this target cell
#                 weighted_avg[target_idx, batch_start:batch_end] = weighted_sum
    
#     # Convert the result to a DataFrame
#     weighted_avg_df = pd.DataFrame(weighted_avg, index=ref_idx_df.index[:num_cells], columns=range(start, end))
    
#     return weighted_avg_df


def read_user_coords(json_path, res=25):
    with open(json_path, 'r') as file:
        json_data = json.load(file)
    
    points = []

    for item in json_data:
        name = item['name']
        triplets = item['triplets']
        count = item['count']

        for i in range(count):
            if res == 25:   # For QUINT workflow results (resolution: 25 um)
                z = triplets[i*3] * 25                  # 25 * qx = z
                x = 13200 - (25 * triplets[i*3 + 1])    # 13200 - (25 * qy) = x
                y = 8000 - (25 * triplets[i*3 + 2])     # 8000 - (25 * qz) = y 
            else:           # Sholud be modified if another method is used to obtain CCF results
                z = triplets[i*3] * 10
                x = 13200 - (25 * triplets[i*3 + 1])
                y = 8000 - (25 * triplets[i*3 + 2])
            points.append([x, y, z, name])
    
    coords = pd.DataFrame(points, columns=['x', 'y', 'z', 'name'])
    return coords


def cfos_sorting(mer_dat_meta, cfoscfos_coords, d):
    pass


def _is_cell_near(cx, cy, cz, xx, yy, zz, d):
    dis = ((cx - xx)**2 + (cy - yy)**2 + (cz - zz)**2)**0.5
    if dis <= d:
        return True
    else:
        return False
    

#################


# def cfos_grouping(json_path, dis, struct_names, pc=None, cancel_event=None):
#     if pc: pc((0, "Reflecting user options ..."))
#     _change_json_path(json_path)
#     _change_d(dis)
#     _change_struct_names(struct_names)

#     if pc: pc((1, "Retrieving all cell location data ..."))
#     if cell_matrix.empty:
#         _generate_all_cell_matrix()

#     if pc: pc((2, "Extracting cell metadata ..."))
#     _generate_filtered_cell_matrix()

#     if pc: pc((3, "Reading user's coordinate information ..."))
#     cfos_coords = _read_user_coords()

#     if pc: pc((4, "Grouping cells ..."))
#     global cfos_o
#     global cfos_x

#     # Using ccf coords
#     for i in filtered_matrix.itertuples():
#         # Adjust 10 um resolution to Allen CCF
#         cx = i.x * 1000
#         cy = i.y * 1000
#         cz = i.z * 1000

#         for j in cfos_coords.itertuples():
#             all_zero = 1
#             if _is_cell_near(cx, cy, cz, j.x, j.y, j.z):
#                 all_zero = 0
#                 break

#         if all_zero:
#             if i.cluster_alias in cfos_x['cluster_alias'].values:
#                 cfos_x.loc[cfos_x['cluster_alias'] == i.cluster_alias, 'n'] += 1
#             else:
#                 new = pd.DataFrame({'cluster_alias': [i.cluster_alias], 'n': [1]})
#                 cfos_x = pd.concat([cfos_x, new], ignore_index = True)
#         else:
#             if i.cluster_alias in cfos_o['cluster_alias'].values:
#                 cfos_o.loc[cfos_o['cluster_alias'] == i.cluster_alias, 'n'] += 1
#             else:
#                 new = pd.DataFrame({'cluster_alias': [i.cluster_alias], 'n': [1]})
#                 cfos_o = pd.concat([cfos_o, new], ignore_index = True)

#     cfos_o = cfos_o.sort_values(by='cluster_alias').reset_index(drop=True)
#     cfos_x = cfos_x.sort_values(by='cluster_alias').reset_index(drop=True)

#     cfos_o = cfos_o.set_index('cluster_alias')
#     cfos_x = cfos_x.set_index('cluster_alias')
#     return cfos_o, cfos_x






