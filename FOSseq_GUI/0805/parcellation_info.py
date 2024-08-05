import csv
import json

def build_tree_from_csv(file_path):
    tree = {}
    nodes = {}
    
    with open(file_path, mode='r') as file:
        reader = csv.DictReader(file)
        
        for row in reader:
            term_label = row['parcellation_term_label']
            term_acronym = row['parcellation_term_acronym']
            term_set_order = int(row['term_set_order'])
            parent_label = row['parent_term_label']
            
            # if term_set_order in (0, 1): continue
            if term_label not in nodes:
                node = {
                    "label": term_label,
                    "acronym": term_acronym,
                    "children": {}
                }
                nodes[term_label] = node
            
                if term_set_order == 0:
                    tree[term_label] = node
                else:
                    parent_node = nodes[parent_label]
                    parent_node['children'][term_label] = node
                
    return tree

def extract_3_level_tree(tree):
    def extract(node, level):
        if level == 0:
            return {}
        
        sub_tree = {}
        for child_label, child_node in node['children'].items():
            sub_tree.update(extract(child_node, level - 1))
        
        return {node['acronym']: sub_tree}
    
    extracted_tree = {}
    for root_label, root_node in tree.items():
        extracted_tree.update(extract(root_node, 5))
        
    return extracted_tree

def save_tree_to_py(tree, file_path):
    with open(file_path, 'w') as file:
        file.write('region_tree = ')
        file.write(json.dumps(tree, indent=4))

# Example usage
file_path = 'C:/programming_data/abc_download_root/metadata/Allen-CCF-2020/20230630/parcellation_to_parcellation_term_membership.csv'
tree = build_tree_from_csv(file_path)
three_level_tree = extract_3_level_tree(tree)
output_path = 'C:/programming_data/GUI_project/FOSseq_GUI/region_tree.py'
save_tree_to_py(three_level_tree, output_path)
