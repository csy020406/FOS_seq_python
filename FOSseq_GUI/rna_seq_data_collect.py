import pandas as pd
import anndata
import os
from pathlib import Path

import abc_atlas_access
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

address = ''
feature_directory = 'WMB-10Xv3'         # WMB-10Xv2         for default
feature_matrix_label = 'WMB-10Xv3-TH'   # WMB-10Xv2-TH      for default
feature_data_file = 'WMB-10Xv3-TH/log2' # WMB-10Xv2-TH/log2 for default

gene_start = 0
gene_end = 1

cell_number = 0

def change_address(s):
    global address
    address = s


# by Option 1, 2, 3
def _change_feature_matrix_label(o1, o2, o3=0):
    global feature_directory
    global feature_matrix_label
    global feature_data_file

    feature_directory = o1
    feature_matrix_label = o2
    if o3 == 0:
        feature_data_file = feature_matrix_label + '/log2'
    else:
        feature_data_file = feature_matrix_label + '/raw'


# by Option 4
def _change_gene_range(s, e):
    global gene_start
    global gene_end

    gene_start = s
    gene_end = e


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


def _update_cell_number(n, is_multi=0):
    global cell_number
    if is_multi:
        cell_number += n
    else:
        cell_number = n


def get_cell_number():
    return cell_number


# ======== Main job =========
def _do_data_collect(pc=None):
    download_base = Path(address)
    abc_cache = AbcProjectCache.from_s3_cache(download_base)

    # Metadata
    cell = abc_cache.get_metadata_dataframe(
        directory = 'WMB-10X',
        file_name = 'cell_metadata',
        dtype = {'cell_label': str}
    )
    cell.set_index('cell_label', inplace = True)

    cluster_details = abc_cache.get_metadata_dataframe(
        directory = 'WMB-taxonomy',
        file_name = 'cluster_to_cluster_annotation_membership_pivoted',
        keep_default_na = False
    )
    cluster_details.set_index('cluster_alias', inplace = True)

    cell_extended = cell.join(cluster_details, on = 'cluster_alias')    # all cell metadata
    pred = (cell_extended['feature_matrix_label'] == feature_matrix_label)
    cell_filtered = cell_extended[pred] # filtered cell metadata

    # Should be controlled by Option 1, 2, 3
    file = abc_cache.get_data_path(directory = feature_directory, file_name = feature_data_file)
    adata = anndata.read_h5ad(file, backed = 'r')

    _update_cell_number(0)
    _update_cell_number(len(adata.obs))

    if pc: pc(-1, 1)   # Done cell metadata generating
    
    def _create_expression_dataframe(ad, gf, cf):
        gdata = ad[:, gf.index].to_df()
        gdata.columns = gf.gene_symbol
        joined = cf.join(gdata)
        return joined
    
    def _aggregate_by_metadata(df, gene_filtered, values, sort = False):
        grouped = df.groupby(values)[gene_filtered].agg(['mean', 'std', 'count']) # don't count NaN
        if sort:
            grouped = grouped.sort_values(by = gene_filtered[0], ascending = False)
        return grouped
    
    global gene_end
    gene_end = min(gene_end, adata.shape[1])    # Compare with the number of total genes
    chunk_size = 200    # Change to improve performance

    # from gene_start to gene_end
    final_agg = None
    for i in range(gene_start, gene_end, chunk_size):
        j = min(i + chunk_size, gene_end)
        # Generate gene_filtered
        gene_filtered = adata.var.iloc[i:j, :]     # gene range change
        asubset = adata[:, gene_filtered.index].to_memory()

        # Gene expression data per cell
        ntexp = _create_expression_dataframe(asubset, gene_filtered, cell_filtered)
        agg = _aggregate_by_metadata(ntexp, gene_filtered.gene_symbol, values=['cluster_alias'])    # (clusters*gene) mean, std, count matrix

        if final_agg is None:
            final_agg = agg
        else:
            final_agg = final_agg.merge(agg, on='cluster_alias', how='outer')
        
        if pc: pc(j, gene_end)
    
    return final_agg


# ======== Main job (multi.ver) =========
def _do_data_collect_multi(o1, o2, o3, pc=None):
    download_base = Path(address)
    abc_cache = AbcProjectCache.from_s3_cache(download_base)

    # Metadata
    cell = abc_cache.get_metadata_dataframe(
        directory = 'WMB-10X',
        file_name = 'cell_metadata',
        dtype = {'cell_label': str}
    )
    cell.set_index('cell_label', inplace = True)

    cluster_details = abc_cache.get_metadata_dataframe(
        directory = 'WMB-taxonomy',
        file_name = 'cluster_to_cluster_annotation_membership_pivoted',
        keep_default_na = False
    )
    cluster_details.set_index('cluster_alias', inplace = True)

    cell_extended = cell.join(cluster_details, on = 'cluster_alias')    # all cell metadata
    pred = pd.Series([False] * len(cell_extended), index=cell_extended.index)

    # Should be controlled by Option 1, 2, 3
    adata_list = []
    _update_cell_number(0)
    for dir in o1:
        for f in o2:
            if f.startswith(dir):
                _change_feature_matrix_label(dir, f, o3)
                file = abc_cache.get_data_path(directory = feature_directory, file_name = feature_data_file)
                adata = anndata.read_h5ad(file, backed = 'r')
                adata_list.append(adata)

                pred |= (cell_extended['feature_matrix_label'] == feature_matrix_label)
                _update_cell_number(len(adata.obs), 1)

    cell_filtered = cell_extended[pred] # filtered cell metadata
    if pc: pc(-1, 1)   # Done cell metadata generating

    def _create_expression_dataframe(ad, gf, cf):
        gdata = ad[:, gf.index].to_df()
        gdata.columns = gf.gene_symbol
        joined = cf.join(gdata)
        return joined
    
    def _aggregate_by_metadata(df, gene_filtered, values, sort = False):
        grouped = df.groupby(values)[gene_filtered].agg(['mean', 'std', 'count']) # don't count NaN
        if sort:
            grouped = grouped.sort_values(by = gene_filtered[0], ascending = False)
        return grouped
    
    global gene_end
    gene_end = min(gene_end, adata_list[0].shape[1])    # Compare with the number of total genes
    chunk_size = 100    # Change to improve performance

    # from gene_start to gene_end
    final_agg = None
    for i in range(gene_start, gene_end, chunk_size):
        j = min(i + chunk_size, gene_end)
        # Generate gene_filtered
        gene_filtered = adata_list[0].var.iloc[i:j, :]     # gene range change

        asubset_list = []
        for a in adata_list:
            asubset = a[:, gene_filtered.index].to_memory()
            asubset_list.append(asubset)

        combined = anndata.concat(asubset_list, axis=0, join='outer')

        # Gene expression data per cell
        ntexp = _create_expression_dataframe(combined, gene_filtered, cell_filtered)
        agg = _aggregate_by_metadata(ntexp, gene_filtered.gene_symbol, values=['cluster_alias'])    # (clusters*gene) mean, std, count matrix

        if final_agg is None:
            final_agg = agg
        else:
            final_agg = final_agg.merge(agg, on='cluster_alias', how='outer')

        if pc: pc(j, gene_end)

    return final_agg


# Excute _do_data_collect or _do_data_collect_multi
def data_collect(o1, o2, o3, o4s, o4e, pc=None):
    if pc: pc(-1, 0)
    _change_gene_range(o4s, o4e)

    if len(o2) == 1:
        _change_feature_matrix_label(o1[0], o2[0], o3)
        agg = _do_data_collect(pc=pc)
    else:
        agg = _do_data_collect_multi(o1, o2, o3, pc=pc)
    return agg


if __name__ == '__main__':
    # main job
    print('main')
    change_address('C:/programming_data/abc_download_root')

    option1 = ['WMB-10Xv2', 'WMB-10Xv3']
    option2 = ['WMB-10Xv2-TH', 'WMB-10Xv2-CTXsp', 'WMB-10Xv3-TH']
    option3 = 0
    option4_start = 0
    option4_end = 10

    agg = data_collect(option1, option2, option3, option4_start, option4_end)
    print(agg)