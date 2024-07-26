import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
import threading
import queue
import pandas as pd
import os

from tester import Tester

import merfish as mer
from region_tree import region_tree


class Sorter:
    def __init__(self, root=None, option_frame=None, cfos_frame=None, cfos_text_widget=None, tester=None):
        self.root = root
        self.option_frame = option_frame
        self.cfos_frame = cfos_frame
        self.cfos_text_widget = cfos_text_widget

        self.tester = tester

        self.label = tk.Label(self.option_frame, text="t-test")
        self.label.pack(padx=10, pady=10)

        self.cfos_button_frame = tk.Frame(self.cfos_frame, bd=0)
        self.cfos_button_frame.pack(padx=10, pady=(0, 10))

        self.load_button = ttk.Button(self.cfos_button_frame, text="Load", command=self.load_file)
        self.load_button.pack(side=tk.LEFT, padx=(0, 5))

        self.save_button = ttk.Button(self.cfos_button_frame, text="Save", command=self.save_file)
        self.save_button.pack(side=tk.RIGHT, padx=(5, 0))

        # USER ABC File Path
        self.file_abc_path = None

        # User JSON File Path
        self.file_json_path = tk.StringVar()
        tk.Button(self.option_frame, text="Select JSON file Path (your cfos coords)", command=self.browse_file).pack(pady=10)
        self.path_label = tk.Label(self.option_frame, textvariable=self.file_json_path)
        self.path_label.pack(pady=10)

        # OPTION 1: USER region of interests
        self.region_info = region_tree
        self.region_frame = tk.Frame(self.option_frame, bd=0)
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

        # OPTION 2: effective radius
        self.radius = tk.DoubleVar(value=25)
        tk.Label(self.option_frame, text="Enter Effective Radius (um)\n").pack()
        tk.Entry(self.option_frame, textvariable=self.radius).pack()

        # Excute button
        self.cfos_button = ttk.Button(self.option_frame, text="Excute cFOS grouping", command=self.start_cfos_grouping)
        self.cfos_button.pack(pady=10)

    def browse_file(self):
        file_selected = filedialog.askopenfilename(
            filetypes=[("JSON files", "*.json")],
            title="Select a JSON file"
        )
        if file_selected:
            self.file_json_path.set(file_selected)

    def set_abc_path(self, new_path):
        self.file_abc_path = new_path

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


    # ====== PROCESS MANAGE ======
    def show_progress_window(self):
        self.progress_window = tk.Toplevel(self.root)
        self.progress_window.title("Processing")
        self.progress_window.geometry("500x150")
        self.progress_window.resizable(False, False)

        self.progress_window.protocol("WM_DELETE_WINDOW", self.on_close_attempt)

        self.progress_window.columnconfigure(0, weight=1)
        self.progress_window.columnconfigure(1, weight=1)
        self.progress_window.rowconfigure(1,weight=1)
        
        self.progress_label = tk.Label(self.progress_window, text="cFOS grouping Excuted.\nPreparing...", 
                                       font=("Arial", 10),
                                       anchor="w",
                                       justify="left")
        self.progress_label.grid(row=0, column=0, columnspan=2, padx=10, pady=(10,0), sticky="w")

        self.progress_var = tk.DoubleVar()

        self.progress_bar = ttk.Progressbar(self.progress_window, variable=self.progress_var, maximum=100)
        self.progress_bar.grid(row=1, column=0, columnspan=2, padx=10, pady=(0, 10), sticky="ew")

        # self.cancel_button = ttk.Button(self.progress_window, text="Cancel", command=self.cancel_task)
        # self.cancel_button.grid(row=2, column=1, padx=10, pady=10, sticky="se")

        # Check periodically if the task is done
        self.root.after(100, self.check_task_done)

    def on_close_attempt(self):
        # Optionally, you can show a message or simply ignore the close attempt
        pass

    def check_task_done(self):
        if self.thread.is_alive():
            self.root.after(100, self.check_task_done)
        else:
            self.progress_bar.stop()
            self.progress_window.destroy()
            self.cfos_button.config(state=tk.NORMAL)

            self.show_cfos_results()

    # def cancel_task(self):
    #     self.queue.put("Cancel")

    #     self.progress_bar.stop()
    #     self.progress_window.destroy()
    #     self.cfos_button.config(state=tk.NORMAL)
                

    # 'current' represents the current task
    # 0: Reflecting user options
    # 1: Retrieving all cell location data from MERFISH database
    # 2: Extracting cell metadata for the selected region
    # 3: Reading user's coordinate infromation
    # 4: cfos grouping
    def update_progress(self, current):
        if current == 0:
            self.progress_label.config(text="cFOS grouping Excuted.\nReflecting user options ...")
            self._gradual_update_progress(self.progress_var.get(), 10)
        elif current == 1:
            self.progress_label.config(text="cFOS grouping Excuted.\nRetrieving all cell location data from MERFISH database ...")
            self._gradual_update_progress(self.progress_var.get(), 50, 800)
        elif current == 2:
            self.progress_label.config(text="cFOS grouping Excuted.\nExtracting cell metadata for the selected region ...")
            self._gradual_update_progress(self.progress_var.get(), 55)
        elif current == 3:
            self.progress_label.config(text="cFOS grouping Excuted.\nReading user's coordinate infromation ...")
            self._gradual_update_progress(self.progress_var.get(), 60)
        else:
            self.progress_label.config(text="cFOS grouping Excuted.\nGrouping cells ...")
            self._gradual_update_progress(self.progress_var.get(), 100, 500)
        
    def _gradual_update_progress(self, from_value, to_value, rate=150):
        step = 1    # Increment step for gradual update
        if from_value < to_value:
            new_value = min(from_value + step, to_value)
            self.progress_var.set(new_value)
            self.progress_window.update_idletasks()     # Update GUI
            # Schedule next update
            self.progress_window.after(rate, self._gradual_update_progress, new_value, to_value)
        else:
            self.progress_var.set(to_value)

    # ====== MAIN JOB ======
    # Examine USER options
    # Excute cfos_grouping
    def start_cfos_grouping(self):
        # ==== Reflect USER OPTIONs ====
        # User File Path
        path = self.file_abc_path
        if not path:
            messagebox.showerror("Error", "Please select a download path.")
            return
        
        # User JSON File Path
        json_path = self.file_json_path.get()
        if not json_path:
            messagebox.showerror("Error", "Please select a JSON file path.")
            return
        
        # OPTION 1
        user_roi = [self.substructure_listbox.get(i) for i in self.substructure_listbox.curselection()]
        if not user_roi:
            messagebox.showerror("Error", "Please select at least one substructure in the set.")
            return
        
        # OPTION 2
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

        self.cfos_button.config(state=tk.DISABLED)

        self.show_progress_window()

        self.queue = queue.Queue()
        
        self.thread = threading.Thread(target=self.cfos_grouping)
        self.thread.start()


    # Excute mer.cfos_grouping
    def cfos_grouping(self):
        try:
            # Reflect USER options
            mer.change_address(self.file_abc_path)
            self.cfos_o, self.cfos_x = mer.cfos_grouping(self.json_path, self.user_radius, self.user_roi, pc=self.update_progress)
        except Exception as e:
            messagebox.showerror("Error", str(e))
        finally:
            self.queue.put("Task Done")

    def show_cfos_results(self):
        self.cfos_text_widget.delete("1.0", "end")

        option_summary = ("========== USER OPTION ==========\n" +
                            "USER path:\t" + f"{self.json_path}\n" +
                            "USER selected region(s):\t" + f"{self.user_roi}\n" +
                            "USER radius:\t" + f"{self.user_radius}\n" +
                            "\n")
        self.cfos_text_widget.insert(tk.END, option_summary)

        cell_number = mer.get_filtered_cell_number()

        if hasattr(self, 'cfos_o') and hasattr(self, 'cfos_x'):
            result_summary = ("========= RESULT SUMMARY =========\n" +
                              "TOTAL cells:\t" + f"{cell_number}\n" +
                              "cFOS (+) cells:\t" + f"{self.cfos_o['n'].sum()}\n" +
                              "cFOS (-) cells:\t" + f"{self.cfos_x['n'].sum()}\n" +
                              "\n")
            self.summary = option_summary + result_summary
            self.cfos_text_widget.insert(tk.END, result_summary)
            self.cfos_text_widget.insert(tk.END, f"[cFOS (+) groups]\n{self.cfos_o}\n")
            self.cfos_text_widget.insert(tk.END, "\n")
            self.cfos_text_widget.insert(tk.END, f"[cFOS (-) groups]\n{self.cfos_x}\n")

            self.tester.update_cfos_data(self.cfos_o, self.cfos_x)
        else:
            self.cfos_text_widget.insert(tk.END, "cFOS grouping was cancelled or failed.")

    def update_display(self):
        pd.set_option('display.multi_sparse', True)

        self.cfos_text_widget.delete(1.0, tk.END)
        if self.summary is not None:
            self.cfos_text_widget.insert(tk.END, f"{self.summary}\n")
        else:
            self.cfos_text_widget.insert(tk.END, "There is no loaded summary data.\n")

        if self.cfos_o is not None and self.cfos_x is not None:
            self.cfos_text_widget.insert(tk.END, f"[cFOS (+) groups]\n{self.cfos_o}\n")
            self.cfos_text_widget.insert(tk.END, "\n")
            self.cfos_text_widget.insert(tk.END, f"[cFOS (-) groups]\n{self.cfos_x}\n")
        else:
            self.cfos_text_widget.insert(tk.END, "There is no loaded cfos data.\n")

    def load_file(self):
        file_path = filedialog.askopenfilename(
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        print(file_path)
        if file_path:
            try:
                if "_cfos_o.csv" in file_path:
                    cfos_o_path = file_path
                    cfos_x_path = file_path.replace("_cfos_o", "_cfos_x")
                elif "_cfos_x.csv" in file_path:
                    cfos_x_path = file_path
                    cfos_o_path = file_path.replace("_cfos_x", "_cfos_o")
                else:
                    messagebox.showerror("Error", "Please check filename format.")
                    return

                self.cfos_o = pd.read_csv(cfos_o_path, index_col=0)
                self.cfos_x = pd.read_csv(cfos_x_path, index_col=0)

                # Check if proper data format
                if self.cfos_o.columns.get_level_values(0) != ['n'] or self.cfos_x.columns.get_level_values(0) != ['n']:
                    messagebox.showerror("Error", "Please check data format.")
                    return
                elif not self.cfos_o.index.to_series().apply(lambda x: isinstance(x, int)).all():
                    messagebox.showerror("Error", "Please check data format.")
                    return
                elif not self.cfos_x.index.to_series().apply(lambda x: isinstance(x, int)).all():
                    messagebox.showerror("Error", "Please check data format.")
                    return

                # Check txt file
                txt_file_path = cfos_o_path.replace("_cfos_o.csv", ".txt")
                if os.path.exists(txt_file_path):
                    with open(txt_file_path, 'r') as txt_file:
                        self.summary = txt_file.read()
                else:
                    self.summary = None

                self.update_display()
                self.tester.update_cfos_data(self.cfos_o, self.cfos_x)
                
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file:\n{e}")

    def save_file(self):
        if self.cfos_o is None or self.cfos_x is None:
            messagebox.showwarning("Save", "Data is not ready yet!")
            return

        self.cfos_o.reset_index(inplace=True)
        self.cfos_x.reset_index(inplace=True)

        file_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if file_path:
            try:
                cfos_o_path = os.path.splitext(file_path)[0] + "_cfos_o" + ".csv"
                cfos_x_path = os.path.splitext(file_path)[0] + "_cfos_x" + ".csv"

                self.cfos_o.to_csv(cfos_o_path, index=False)
                self.cfos_x.to_csv(cfos_x_path, index=False)

                if self.summary:
                    txt_file_path = os.path.splitext(file_path)[0] + ".txt"
                    with open(txt_file_path, 'w') as txt_file:
                        txt_file.write(self.summary)
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save file:\n{e}")
   

if __name__ == "__main__":
    # root = tk.Tk()
    # root.title("FOS seq")
    # root.geometry("900x600")

    # # Frames
    # root.rowconfigure(0, weight=1)
    # root.rowconfigure(1, weight=1)
    # root.columnconfigure(0, weight=0)
    # root.columnconfigure(1, weight=1)
    # root.columnconfigure(2, weight=1)

    # collector_frame = tk.Frame(root, borderwidth=1, relief="solid", width=500)
    # collector_frame.grid(row=0, column=0, rowspan=2, padx=10, pady=10, sticky='news')
    # collector_frame.grid_propagate(False)

    # tester_frame = tk.Frame(root, borderwidth=1, relief="solid")
    # tester_frame.grid(row=0, column=1, padx=10, pady=10, sticky='news')

    # sorter_option_frame = tk.Frame(root, borderwidth=1, relief="solid")
    # sorter_option_frame.grid(row=1, column=1, padx=10, pady=10, sticky='news')

    # sorter_result_frame = tk.Frame(root, borderwidth=1, relief="solid")
    # sorter_result_frame.grid(row=0, column=2, rowspan=2, padx=10, pady=10, sticky='news')

    # # Text widget
    # collection_text_widget = tk.Text(tester_frame, wrap=tk.WORD)
    # collection_text_widget.pack(padx=10, pady=10, expand=True, fill=tk.BOTH)

    # Classes
    # sorter_app = Sorter(root, sorter_option_frame, sorter_result_frame)
    
    # root.mainloop()

    sorter_app = Sorter()
    sorter_app.load_file()