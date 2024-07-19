import os
import pandas as pd
from pathlib import Path
import numpy as np
import json
import matplotlib.pyplot as plt
import SimpleITK as sitk
import pathlib
from ast import literal_eval
from scipy import stats
import matplotlib.pyplot as plt

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

address = ''
json_path = ''
struct_name = ''

d = 0.25  # for _is_cell_near

cell_matrix = pd.DataFrame()
filtered_matrix = pd.DataFrame()

def change_address(s):
    global address
    address = s


def change_json_path(s):
    global json_path
    json_path = s


def change_struct_name(s):
    global struct_name
    struct_name = s


def change_d(n):
    global d
    d = n


def generate_all_cell_matrix():
    download_base = Path(address)
    abc_cache = AbcProjectCache.from_s3_cache(download_base)

    cell = abc_cache.get_metadata_dataframe(directory='MERFISH-C57BL6J-638850', file_name='cell_metadata_with_cluster_annotation')
    cell.set_index('cell_label', inplace=True)

    # extract your needs from cell metadata
    cell_extract = cell.loc[:, ['brain_section_label',
                                'cluster_alias',
                                'average_correlation_score']]

    # reconstructed coordinates
    reconstructed_coords = abc_cache.get_metadata_dataframe(
        directory='MERFISH-C57BL6J-638850-CCF',
        file_name='reconstructed_coordinates',
        dtype={"cell_label": str}
    )
    reconstructed_coords.rename(columns={'x': 'x_reconstructed',
                                         'y': 'y_reconstructed',
                                         'z': 'z_reconstructed'},
                                inplace=True)
    reconstructed_coords.set_index('cell_label', inplace=True)
    cell_joined = cell_extract.join(reconstructed_coords, how='inner')

    # parcellation annotation
    parcellation_annotation = abc_cache.get_metadata_dataframe(directory='Allen-CCF-2020',
                                                               file_name='parcellation_to_parcellation_term_membership_acronym')
    parcellation_annotation.set_index('parcellation_index', inplace=True)
    parcellation_annotation.columns = ['parcellation_%s'% x for x in  parcellation_annotation.columns]
    parcellation_annotation = parcellation_annotation.loc[:, ['parcellation_division',
                                                              'parcellation_structure',
                                                              'parcellation_substructure']]
    cell_joined = cell_joined.join(parcellation_annotation, on='parcellation_index')

    global cell_matrix
    cell_matrix = cell_joined


def generate_filtered_cell_matrix():
    pred = (cell_matrix['parcellation_substructure'] == struct_name)

    # use copy() to avoid SettingWithCopyWarning error
    global filtered_matrix
    filtered_matrix = cell_matrix[pred].copy()

    print(len(filtered_matrix))


def _read_user_coords():
    with open(json_path, 'r') as file:
        json_data = json.load(file)
    
    points = []

    for item in json_data:
        name = item['name']
        triplets = item['triplets']
        count = item['count']

        for i in range(count):
            # Should be MODIFIED by CCF converts
            x = triplets[i*3]
            y = triplets[i*3 + 1]
            z = triplets[i*3 + 2]
            points.append([x, y, z, name])
    
    coords = pd.DataFrame(points, columns=['x', 'y', 'z', 'name'])
    return coords


def _is_cell_near(cx, cy, cz, xx, yy, zz):
    dis = ((cx - xx)**2 + (cy - yy)**2 + (cz - zz)**2)**0.5
    if dis <= d:
        return True
    else:
        return False


# ======== Main job =========
def _cfos_grouping():
    cfos_coords = _read_user_coords()

    cfos_o = pd.DataFrame(columns=['cluster_alias', 'number'], dtype=int)
    cfos_x = pd.DataFrame(columns=['cluster_alias', 'number'], dtype=int)

    # uses reconstructed coords
    for i in filtered_matrix.itertuples():
        cx = i.x_reconstructed
        cy = i.y_reconstructed
        cz = i.z_reconstructed

        for j in cfos_coords.itertuples():
            all_zero = 1
            if _is_cell_near(cx, cy, cz, j.x, j.y, j.z):
                all_zero = 0
                break

        if all_zero:
            if i.cluster_alias in cfos_x['cluster_alias'].values:
                cfos_x.loc[cfos_x['cluster_alias'] == i.cluster_alias, 'number'] += 1
            else:
                new = pd.DataFrame({'cluster_alias': [i.cluster_alias], 'number': [1]})
                cfos_x = pd.concat([cfos_x, new], ignore_index = True)
        else:
            if i.cluster_alias in cfos_o['cluster_alias'].values:
                cfos_o.loc[cfos_o['cluster_alias'] == i.cluster_alias, 'number'] += 1
            else:
                new = pd.DataFrame({'cluster_alias': [i.cluster_alias], 'number': [1]})
                cfos_o = pd.concat([cfos_o, new], ignore_index = True)

    cfos_o = cfos_o.sort_values(by='cluster_alias').reset_index(drop=True)
    cfos_x = cfos_x.sort_values(by='cluster_alias').reset_index(drop=True)

    print(cfos_o.head(5))
    print(cfos_x.head(5))

def t_test():
    print()

if __name__ == '__main__':
    print('main')
    change_address('C:/programming_data/abc_download_root')
    change_json_path('C:/programming_data/test_resized/Coordinates/Test_resized_3D_combined.json')

    struct = 'LH'
    change_struct_name(struct)

    # optional
    # dis = 0.25
    # change_d(dis)

    generate_all_cell_matrix()
    generate_filtered_cell_matrix()
    _cfos_grouping()
    t_test()