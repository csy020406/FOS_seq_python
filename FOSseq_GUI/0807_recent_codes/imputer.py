import threading
import queue
from pathlib import Path

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
from sklearn.preprocessing import normalize
from sklearn.neighbors import NearestNeighbors

import anndata
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt


from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

from region_tree import region_tree
import training as tr
import testing as tt


class Imputer:
    def __init__(
            self, k=15, method='cosine',
            root=None, abc_frame=None, train_frame=None, test_frame=None, result_frame=None
    ):
        # ====================
        # Windows
        # ====================
        self.root = root
        self.abc_frame = abc_frame
        self.train_frame = train_frame
        self.test_frame = test_frame
        self.result_frame = result_frame

        # Result Gene Range Frame
        self.result_gene_frame = tk.Frame(self.result_frame, bd=0)
        self.result_gene_frame.pack(padx=10, pady=(0,10))

        # Result Text Widget
        self.result_text_widget = tk.Text(result_frame, wrap=tk.WORD)
        self.result_text_widget.pack(padx=10, pady=10, expand=True, fill=tk.BOTH)

        # Result Matplotlib Figure
        self.fig = Figure(figsize=(10, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)

        self.result_canvas = FigureCanvasTkAgg(self.fig, master=self.result_frame)
        self.result_canvas.draw()
        self.result_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # ==============================
        # Variables
        # ==============================
        self.abc_cache = None

        self.queue = queue.Queue()

        self.nn = NearestNeighbors(n_neighbors=k, metric=method)

        # Training Results
        self.mer_adata = None
        self.ten_adata_list = None
        self.ten_cell_num = 0
        self.gene_markers = None
        self.ref_dat_df = None
        self.ref_idx_df = None

        # Testing Results
        self.mer_metadata = None
        self.mer_dat_meta = None
        self.mer_dat_df = None
        self.knn_indices = None
        self.knn_distances = None
        self.knn_weights = None
        self.cfos_ox = None

        # Imputing Results
        self.weighted_average_df = None

        # ==============================
        # ABC Frame (abc_frame)
        # ==============================
        self.abc_path_label = tk.Label(self.abc_frame, text="")
        self.abc_path_label.pack(pady=10)
        ttk.Button(abc_frame, text="Select ABC Download Path(*)", command=self.browse_folder).pack(pady=10)
        
        # ==============================
        # Training Frame (train_frame)
        # ==============================
        tk.Label(self.train_frame, text="Model Training").pack()

        # OPTION 1: 10Xv2, 10Xv3, 10XMulti
        self.option1_vars = [tk.IntVar() for _ in range(3)]
        option1 = ["WMB-10Xv2", "WMB-10Xv3", "WMB-10XMulti"]
        tk.Label(self.train_frame, text="Select Directories (multiple choice):").pack()
        for i, option in enumerate(option1):
            tk.Checkbutton(self.train_frame, text=option, variable=self.option1_vars[i], command=self.show_datafiles).pack(anchor='w')

        # OPTION 2: file names
        tk.Label(self.train_frame, text="Select Files:").pack()
        self.option2_listbox = tk.Listbox(self.train_frame, selectmode=tk.MULTIPLE, exportselection=False)
        self.option2_listbox.pack(fill=tk.BOTH, expand=True)
        self.scrollbar = tk.Scrollbar(self.option2_listbox, orient=tk.VERTICAL)
        self.scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.option2_listbox.config(yscrollcommand=self.scrollbar.set)
        self.scrollbar.config(command=self.option2_listbox.yview)

        # OPTION 3: data type
        self.option3_var = tk.IntVar()
        option3 = [("log2", 0), ("raw", 1)]
        tk.Label(self.train_frame, text="Select Data Type:").pack()
        for option, val in option3:
            tk.Radiobutton(self.train_frame, text=option, variable=self.option3_var, value=val).pack(anchor='w')

        # EXCUTE TRAINING button
        self.exe_tr_button = ttk.Button(self.train_frame, text="Excute Training", command=self.start_training)
        self.exe_tr_button.pack(pady=10)
        self.exe_tr_button.config(state=tk.DISABLED)


        # ==============================
        # Testing Frame (test_frame)
        # ==============================
        # OPTION 4: USER region of interests
        tk.Label(self.test_frame, text="Select your Regions of Interests\n").pack()

        self.region_info = region_tree
        self.region_frame = tk.Frame(self.test_frame, bd=0)
        self.region_frame.pack(padx=10, pady=(0, 10))

        self.division_listbox = tk.Listbox(self.region_frame, selectmode=tk.MULTIPLE, exportselection=False)
        self.division_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        self.structure_listbox = tk.Listbox(self.region_frame, selectmode=tk.MULTIPLE, exportselection=False)
        self.structure_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        self.substructure_listbox = tk.Listbox(self.region_frame, selectmode=tk.MULTIPLE, exportselection=False)
        self.substructure_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.populate_divisions()

        # Bind selection events
        self.division_listbox.bind('<<ListboxSelect>>', self.update_structures)
        self.structure_listbox.bind('<<ListboxSelect>>', self.update_substructures)

        # OPTION 5: effective radius
        self.radius = tk.DoubleVar(value=25)
        tk.Label(self.test_frame, text="Enter Effective Radius (um)\n").pack()
        tk.Entry(self.test_frame, textvariable=self.radius).pack()

        # User JSON File Path
        self.file_json_path = tk.StringVar()
        ttk.Button(self.test_frame, text="Select JSON file Path (your cfos coords)", command=self.browse_file).pack(pady=10)
        self.path_label = tk.Label(self.test_frame, textvariable=self.file_json_path)
        self.path_label.pack(pady=10)

        # EXCUTE TESTING button
        self.exe_test_button = ttk.Button(self.test_frame, text="Excute Testing", command=self.start_testing)
        self.exe_test_button.pack(pady=10)

        # ==========================================================================================
        # Result Frame (result_frame, result_gene_frame, result_text_widget, result_canvas)
        # ==========================================================================================
        # OPTION 6: gene range
        self.range_start = tk.IntVar()
        self.range_end = tk.IntVar()
        tk.Label(self.result_gene_frame, text="Enter Range (0-32284)\n(-1, -1) for all genes, but not recommended:").pack()
        tk.Entry(self.result_gene_frame, textvariable=self.range_start).pack()
        tk.Entry(self.result_gene_frame, textvariable=self.range_end).pack()

        # EXCUTE IMPUTING button
        self.exe_impute_button = ttk.Button(self.result_gene_frame, text="Excute Imputing", command=self.start_imputing)
        self.exe_impute_button.pack(pady=10)
    
    # =============================
    # ABC Frame Functions
    # =============================
    def browse_folder(self):
        folder_selected = filedialog.askdirectory()
        if folder_selected:
            self.abc_path_label.config(text=folder_selected)
            
            download_base = Path(folder_selected)
            self.abc_cache = AbcProjectCache.from_s3_cache(download_base)

            self.exe_tr_button.config(state=tk.NORMAL)

    # =============================
    # Training Functions
    # =============================
    def show_datafiles(self):
        # Get option 1
        user_option1 = [opt for var, opt in zip(self.option1_vars, ["WMB-10Xv2", "WMB-10Xv3", "WMB-10XMulti"]) if var.get() == 1]
        if not user_option1:
            return
        datafiles = tr.list_file_name(user_option1)

        # Update listbox
        self.option2_listbox.delete(0, tk.END)
        for file_name in datafiles:
            self.option2_listbox.insert(tk.END, file_name)

    def start_training(self):
        # ==== Reflect USER OPTIONs ====
        # OPTION 1
        user_option1 = [opt for var, opt in zip(self.option1_vars, ["WMB-10Xv2", "WMB-10Xv3", "WMB-10XMulti"]) if var.get() == 1]
        if not user_option1:
            messagebox.showerror("Error", "Please select at least one option in the set.")
            return
        
        # OPTION 2
        user_option2 = [self.option2_listbox.get(i) for i in self.option2_listbox.curselection()]
        if not user_option2:
            messagebox.showerror("Error", "Please select at least one file in the set.")
            return
        
        # OPTION 3
        user_option3 = self.option3_var.get()   # 0 for "log2", 1 for "raw"

        # ==== All Options are ready ====
        self.user_o1 = user_option1
        self.user_o2 = user_option2
        self.user_o3 = user_option3

        self.thread = threading.Thread(target=self.training)
        self.thread.start()

    def training(self):
        try:
            if self.mer_adata is None:
                self.mer_adata = tr.get_mer_gene_data(self.abc_cache)
            self.ten_adata_list, self.ten_cell_num = tr.get_10x_adata_list(
                self.abc_cache,
                self.user_o1, self.user_o2, self.user_o3
            )
            self.gene_markers = tr.get_common_genes(self.mer_adata, self.ten_adata_list[0])
            self.ref_dat_df, self.ref_idx_df = tr.get_ref_df(self.gene_markers, self.ten_adata_list)

            norm_ref_dat = normalize(self.ref_dat_df, norm='l2')
            self.nn.fit(norm_ref_dat)
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.queue.put("Error")
        finally:
            self.queue.put("Task Done")
    
    # =============================
    # Testing Functions
    # =============================
    def browse_file(self):
        file_selected = filedialog.askopenfilename(
            filetypes=[("JSON files", "*.json")],
            title="Select a JSON file"
        )
        if file_selected:
            self.file_json_path.set(file_selected)

    def populate_divisions(self):
        self.division_listbox.delete(0, tk.END)

        for organ, categories in self.region_info.items():
            for category, divisions in categories.items():
                for division in divisions:
                    self.division_listbox.insert(tk.END, division)
                    
        self.structure_listbox.delete(0, tk.END)
        self.substructure_listbox.delete(0, tk.END)
    
    def update_structures(self, event):
        selected_divisions = [self.division_listbox.get(i) for i in self.division_listbox.curselection()]
        self.structure_listbox.delete(0, tk.END)
        self.substructure_listbox.delete(0, tk.END)

        for organ, categories in self.region_info.items():
            for category, divisions in categories.items():
                for division in divisions:
                    if division in selected_divisions:
                        for structure in self.region_info[organ][category][division].keys():
                            self.structure_listbox.insert(tk.END, structure)

    def update_substructures(self, event):
        selected_structures = [self.structure_listbox.get(i) for i in self.structure_listbox.curselection()]
        self.substructure_listbox.delete(0, tk.END)

        for organ, categories in self.region_info.items():
            for category, divisions in categories.items():
                for division, structures in divisions.items():
                    for structure, substructures in structures.items():
                        if structure in selected_structures:
                            for substructure in substructures:
                                self.substructure_listbox.insert(tk.END, substructure)


    def start_testing(self):
        # ==== Reflect USER OPTIONs ====
        # User JSON File Path
        json_path = self.file_json_path.get()
        if not json_path:
            messagebox.showerror("Error", "Please select a JSON file path.")
            return
        
        # OPTION 4
        user_roi = [self.substructure_listbox.get(i) for i in self.substructure_listbox.curselection()]
        if not user_roi:
            messagebox.showerror("Error", "Please select at least one substructure in the set.")
            return
        
        # OPTION 5
        radius = self.radius.get()
        if radius < 0:
            messagebox.showerror("Error", "Please enter a valid double.")
            return
        if radius > 250:
            messagebox.showwarning("Warning", "Your specified value may be too large.")
        
        # ==== All Options are ready ====
        self.json_path = json_path
        self.user_roi = user_roi
        self.user_radius = radius

        self.thread = threading.Thread(target=self.testing)
        self.thread.start()

    def testing(self):
        try:
            if self.mer_metadata is None:
                self.mer_metadata = tt.get_mer_metadata(self.abc_cache)
            self.mer_dat_meta, self.mer_dat_df = tt.get_mer_dat(self.mer_metadata, self.mer_adata, self.gene_markers, self.user_roi)
            
            norm_dat = normalize(self.mer_dat_df, norm='l2')
            self.knn_indices, self.knn_distances = self.nn.kneighbors(norm_dat)

            cfos_coords = tt.read_user_coords(self.json_path)
            self.cfos_ox = tt.cfos_sorting(self.user_radius)
        except Exception as e:
            messagebox.showerror("Error", str(e))
        finally:
            self.queue.put("Task Done")


    # =============================
    # Imputing Functions
    # =============================
    def start_imputing(self):
        # OPTION 6
        start = self.range_start.get()
        end = self.range_end.get()
        if (start == -1 and end == -1):
            start = 0
            end = 33000     # For mis-calculation
        elif not (0 <= start <= 32284 and 0 <= end <= 33000 and start <= end):
            messagebox.showerror("Error", "Please enter valid integers in the range 1-32285.")
            return
        
        self.gene_start = start
        self.gene_end = end

        self.thread = threading.Thread(target=self.imputing)
        self.thread.start()

    def imputing(self):
        try:
            # TODO: 변경!!!!!!!!!!
            weighted_average = tt.cal_weighted_average(self.knn_indices, self.knn_distances, self.ref_dat_df)
            
            self.weighted_average_df = pd.DataFrame(weighted_average, index=self.mer_dat_df.index, columns=self.ref_dat_df.columns)


        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.queue.put("Error")
        finally:
            self.queue.put("Task Done")


    def generate_seq_data(self, directory, file_name, roi):
        # self.rna_seq_data
        file = self.abc_cache.get_data_path(
            directory = directory,
            file_name = file_name
        )
        self.rna_seq_data = anndata.read_h5ad(file, backed='r')

        # self.mer_seq_data

        cell = abc_cache.get_metadata_dataframe(directory='MERFISH-C57BL6J-638850', file_name='cell_metadata_with_cluster_annotation')
        cell.set_index('cell_label', inplace=True)

        # Extract your needs from cell metadata
        cell_extract = cell.loc[:, ['brain_section_label',
                                    'cluster_alias']]           # TODO: 'average_correlation_score'?

        # CCF coordinates
        coords = abc_cache.get_metadata_dataframe(
            directory='MERFISH-C57BL6J-638850-CCF',
            file_name='ccf_coordinates',
            dtype={"cell_label": str}
        )
        coords.set_index('cell_label', inplace=True)
        cell_joined = cell_extract.join(coords, how='inner')

        # Parcellation annotation
        parcellation_annotation = abc_cache.get_metadata_dataframe(directory='Allen-CCF-2020',
                                                                file_name='parcellation_to_parcellation_term_membership_acronym')
        parcellation_annotation.set_index('parcellation_index', inplace=True)
        parcellation_annotation.columns = ['parcellation_%s'% x for x in  parcellation_annotation.columns]
        parcellation_annotation = parcellation_annotation.loc[:, ['parcellation_substructure']]         # TODO: need 'parcellation_division',' parcellation_structure'?
        cell_joined = cell_joined.join(parcellation_annotation, on='parcellation_index')

        global cell_matrix
        cell_matrix = cell_joined





