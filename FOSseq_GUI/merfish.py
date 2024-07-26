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

# NOTE:
# The following code aligns the CCF as follows (in um):
# x: increases backward (posterior direction) from 0 to 13200
# y: increases downward (inferior direction) from 0 to 8000
# z: increases to the right from 0 to 11400

# MERFISH CCF_coordinates (mx, my, mz) is transformed as follows:
# (x, y, z) = (mx * 1000, my * 1000, mz * 1000)
# (10 um resolution)

# QUINT workflow result CCF (qx, qy, qz) is transformed as follows:
# (x, y, z) = (13200 - (25 * qy), 8000 - (25 * qz), 25 * qx)
# (25 um resolution)


address = ''
json_path = ''
struct_names = []

d = 25      # for _is_cell_near
count = 10  # for calculate_weighted_stats

cell_matrix = pd.DataFrame()
filtered_matrix = pd.DataFrame()
cfos_o = pd.DataFrame(columns=['cluster_alias', 'n'], dtype=int)
cfos_x = pd.DataFrame(columns=['cluster_alias', 'n'], dtype=int)

filtered_len = 0


def change_address(s):
    global address
    address = s


def _change_json_path(s):
    global json_path
    json_path = s


def _change_struct_names(s):
    global struct_names
    struct_names = s


def _change_d(n):
    global d
    d = n


def change_count(n):
    global count
    count = n


def _generate_all_cell_matrix():
    download_base = Path(address)
    abc_cache = AbcProjectCache.from_s3_cache(download_base)

    cell = abc_cache.get_metadata_dataframe(directory='MERFISH-C57BL6J-638850', file_name='cell_metadata_with_cluster_annotation')
    cell.set_index('cell_label', inplace=True)

    # extract your needs from cell metadata
    cell_extract = cell.loc[:, ['brain_section_label',
                                'cluster_alias',
                                'average_correlation_score']]

    # ccf coordinates
    coords = abc_cache.get_metadata_dataframe(
        directory='MERFISH-C57BL6J-638850-CCF',
        file_name='ccf_coordinates',
        dtype={"cell_label": str}
    )
    coords.set_index('cell_label', inplace=True)
    cell_joined = cell_extract.join(coords, how='inner')

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


def _generate_filtered_cell_matrix():
    pred = (cell_matrix['parcellation_substructure'].isin(struct_names))

    # use copy() to avoid SettingWithCopyWarning error
    global filtered_matrix
    filtered_matrix = cell_matrix[pred].copy()
    

def get_filtered_cell_number():
    global filtered_len
    filtered_len = len(filtered_matrix)

    return filtered_len


def _read_user_coords(res=25):
    with open(json_path, 'r') as file:
        json_data = json.load(file)
    
    points = []

    for item in json_data:
        name = item['name']
        triplets = item['triplets']
        count = item['count']

        for i in range(count):
            if res == 25:   # for QUINT workflow results (resolution: 25 um)
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


def _is_cell_near(cx, cy, cz, xx, yy, zz):
    dis = ((cx - xx)**2 + (cy - yy)**2 + (cz - zz)**2)**0.5
    if dis <= d:
        return True
    else:
        return False


# ======== Main job =========
def cfos_grouping(json_path, dis, struct_names, pc=None, cancel_event=None):
    if pc: pc((0, "Reflecting user options ..."))
    _change_json_path(json_path)
    _change_d(dis)
    _change_struct_names(struct_names)

    if pc: pc((1, "Retrieving all cell location data ..."))
    if cell_matrix.empty:
        _generate_all_cell_matrix()

    if pc: pc((2, "Extracting cell metadata ..."))
    _generate_filtered_cell_matrix()

    if pc: pc((3, "Reading user's coordinate information ..."))
    cfos_coords = _read_user_coords()

    if pc: pc((4, "Grouping cells ..."))
    global cfos_o
    global cfos_x

    # using ccf coords
    for i in filtered_matrix.itertuples():
        # Adjust 10 um resolution to Allen CCF
        cx = i.x * 1000
        cy = i.y * 1000
        cz = i.z * 1000

        for j in cfos_coords.itertuples():
            all_zero = 1
            if _is_cell_near(cx, cy, cz, j.x, j.y, j.z):
                all_zero = 0
                break

        if all_zero:
            if i.cluster_alias in cfos_x['cluster_alias'].values:
                cfos_x.loc[cfos_x['cluster_alias'] == i.cluster_alias, 'n'] += 1
            else:
                new = pd.DataFrame({'cluster_alias': [i.cluster_alias], 'n': [1]})
                cfos_x = pd.concat([cfos_x, new], ignore_index = True)
        else:
            if i.cluster_alias in cfos_o['cluster_alias'].values:
                cfos_o.loc[cfos_o['cluster_alias'] == i.cluster_alias, 'n'] += 1
            else:
                new = pd.DataFrame({'cluster_alias': [i.cluster_alias], 'n': [1]})
                cfos_o = pd.concat([cfos_o, new], ignore_index = True)

    cfos_o = cfos_o.sort_values(by='cluster_alias').reset_index(drop=True)
    cfos_x = cfos_x.sort_values(by='cluster_alias').reset_index(drop=True)

    cfos_o = cfos_o.set_index('cluster_alias')
    cfos_x = cfos_x.set_index('cluster_alias')
    return cfos_o, cfos_x


# Used in _get_weighted_stats
def _calculate_weighted_stats(gene_df, cfos_df):
    if gene_df.empty:
        return np.nan, np.nan, np.nan
    gene_df = gene_df.T
    gene_df.columns = gene_df.columns.droplevel(0)

    # filtering by count
    gene_df = gene_df[gene_df['count'] >= count]
    cfos_df = cfos_df.loc[gene_df.index]

    total_weight = cfos_df['n'].sum()
    weighted_mean = (gene_df['mean'] * cfos_df['n']).sum() / total_weight

    variance_numer = (cfos_df['n'] * (gene_df['std'] ** 2 + (gene_df['mean'] - weighted_mean) ** 2)).sum()
    weighted_variance = variance_numer / total_weight
    weighted_std = np.sqrt(weighted_variance)

    return weighted_mean, weighted_std, total_weight


# Called by t_test
def _get_weighted_stats(multi_gene_df, cfos_df):
    try:
        cluster_order = cfos_df.index
        multi_gene_df = multi_gene_df.loc[cluster_order]
    except KeyError as e:
        missing_indices = set(cluster_order) - set(multi_gene_df.index)
        print(f"Warning: The following indices are not in multi_gene_df and were ignored: {missing_indices}")
        multi_gene_df = multi_gene_df.reindex(cluster_order)

    # transpose to prevent FUTUREWARNING
    result = multi_gene_df.T.groupby(level=0).apply(lambda x: pd.Series(_calculate_weighted_stats(x, cfos_df)))
    result.columns = ['mean', 'std', 'n']
    result.rename_axis('gene', inplace=True)
    return result


# rna_seq is rna_seq_collection result from rsd
# tester would send it
def t_test(rna_seq=None, cfos_o=None, cfos_x=None, pc=None):
    weighted_stats_o = _get_weighted_stats(rna_seq, cfos_o)
    weighted_stats_x = _get_weighted_stats(rna_seq, cfos_x)

    t_test_results = pd.DataFrame(index=weighted_stats_o.index, columns=['fold_change', 'p_value'])

    total_loops = len(weighted_stats_o.index)
    for i, gene in enumerate(weighted_stats_o.index, 1):
        mean_o, std_o, n_o = weighted_stats_o.loc[gene]
        mean_x, std_x, n_x = weighted_stats_x.loc[gene]

        t_stat, p_value = stats.ttest_ind_from_stats(mean_o, std_o, n_o, mean_x, std_x, n_x)

        t_test_results.loc[gene, 'fold_change'] = mean_o - mean_x
        t_test_results.loc[gene, 'p_value'] = p_value

        if pc: pc(i, total_loops)

    t_test_results.rename_axis('gene', inplace=True)
    return t_test_results


def draw_volcano_plot(results):
    plt.figure(figsize=(10,6))

    fold_changes = results['fold_change']
    p_values = results['p_value']

    # neg_log_p_values = -np.log10(p_values)
    neg_log_p_values = p_values

    plt.scatter(fold_changes, neg_log_p_values, alpha=0.75)

    for gene, fold_change, neg_log_p_value in zip(results.index, fold_changes, neg_log_p_values):
        if abs(fold_change) > 1 or neg_log_p_value > 2:
            plt.text(fold_change, neg_log_p_value, gene, fontsize=0)

    plt.xlabel('Fold Change')
    plt.ylabel('-log10(p-value)')
    plt.title('Volcano Plot')
    plt.axhline(y=-np.log10(0.05), color='r', linestyle='--')  # p-value < 0.05 기준선
    plt.axvline(x=1, color='r', linestyle='--')  # Fold change > 1 기준선
    plt.axvline(x=-1, color='r', linestyle='--')  # Fold change < -1 기준선
    plt.show(block=True)


if __name__ == '__main__':
    # print('main')
    # change_address('C:/programming_data/abc_download_root')
    # change_json_path('C:/programming_data/test_output_ABA/Coordinates/Test_resized_ABA_3D_combined.json')

    # struct = 'LH'
    # change_struct_names(struct)

    # # optional
    # # dis = 0.25
    # # change_d(dis)
    # print("generating_all_cell_matrix")
    # generate_all_cell_matrix()
    # print("generating_filtered_cell_matrix")
    # generate_filtered_cell_matrix()
    # print("cfos_grouping")
    # cfos_grouping()
    # print("t_test")
    # t_test()

    # weight cal test
    data = {
    ('Gene1', 'mean'): [8.923037, 8.656472, 8.849827, 8.4526],
    ('Gene1', 'std'): [np.nan, 0.305846, 0.390072, 0.35243],
    ('Gene1', 'count'): [1, 2, 15, 34],
    ('Gene2', 'mean'): [0.000000, 2.138494, 0.000000, 4.5424],
    ('Gene2', 'std'): [np.nan, 3.703981, 0.000000, 5.451],
    ('Gene2', 'count'): [15, 3, 12, 54]
    }
    index = [565, 569, 579, 600]
    multi_gene_df = pd.DataFrame(data, index=index)
    multi_gene_df.index.name = 'cluster_alias'

    cfos_data1 = {
        'cluster_alias': [565, 569, 579],
        'n': [10, 15, 20]
    }

    cfos_data2 = {
        'cluster_alias': [565, 569, 579, 600, 625],
        'n': [12, 20, 35, 45, 50]
    }
    
    cfos_df1 = pd.DataFrame(cfos_data1)
    cfos_df1 = cfos_df1.set_index('cluster_alias')

    cfos_df2 = pd.DataFrame(cfos_data2)
    cfos_df2 = cfos_df2.set_index('cluster_alias')

    print("gene df\n")
    print(multi_gene_df.head(5))
    print("cfos_df\n")
    print(cfos_df1.head(5))
    print(cfos_df2.head(5))
    # 함수 호출

    result = t_test(multi_gene_df, cfos_df1, cfos_df2)
    draw_volcano_plot(result)