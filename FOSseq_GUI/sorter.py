import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
import threading
import queue
import pandas as pd

from tester import Tester

import merfish as mer
from region_tree import region_tree


class Sorter:
    def __init__(self, root, option_frame, cfos_frame, cfos_text_widget, tester):
        self.root = root
        self.option_frame = option_frame
        self.cfos_frame = cfos_frame
        self.cfos_text_widget = cfos_text_widget

        self.tester = tester

        self.label = tk.Label(self.option_frame, text="t-test")
        self.label.pack(padx=10, pady=10)

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
    def show_warning_window(self, file_name, file_size):
        self.dialog = tk.Toplevel(self, root)
        self.dialog.title("Warning")

        self.dialog.geometry("500x200")

        message = (
            f"The required sizes of each ABC atlas metadata file '{file_name}' is {file_size} GB.\n\n"
            "YOU CANNOT INTERRUPT THIS PROCESS ONCE STARTED, so\n"
            "it is reccommended to DOUBLE-CHECK the ABC download path. If you have not previously\n"
            "downloaded this file to the same path, the download will begin and may take a long time."
        )
        self.warning_label = tk.Label(self.dialog, text=message, padx=10, pady=10, wraplength=480)
        self.warning_label.pack()

        self.button_frame = tk.Frame(self.dialog)
        self.button_frame.pack(pady=10)

        self.cancel_button = tk.Button(self.button_frame, text="Cancel", command=self.on_cancel, bg="lightblue")
        self.cancel_button.pack(side=tk.LEFT, padx=5)
        self.continue_button = tk.Button(self.button_frame, text="Continue", command=self.on_continue)
        self.continue_button.pack(side=tk.LEFT, padx=5)

        self.result = None

        self.root.wait_window(self.dialog)
        
        return self.result
    
    def on_cancel(self):
        self.result = 1
        self.dialog.destroy()

    def on_continue(self):
        self.result = 0
        self.dialog.destroy()
    
    def show_progress_window(self):
        self.progress_window = tk.Toplevel(self.root)
        self.progress_window.title("Processing")
        self.progress_window.geometry("500x150")
        self.progress_window.resizable(False, False)

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

        self.cancel_button = ttk.Button(self.progress_window, text="Cancel", command=self.cancel_task)
        self.cancel_button.grid(row=2, column=1, padx=10, pady=10, sticky="se")

        # Check periodically if the task is done
        self.root.after(100, self.check_task_done)

    def check_task_done(self):
        if self.thread.is_alive():
            self.root.after(100, self.check_task_done)
        else:
            self.progress_bar.stop()
            self.progress_window.destroy()
            self.cfos_button.config(state=tk.NORMAL)

            self.show_cfos_results()

    def cancel_task(self):
        self.queue.put("Cancel")

        self.progress_bar.stop()
        self.progress_window.destroy()
        self.cfos_button.config(state=tk.NORMAL)
                

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
                              "cFOS (+) cells:\t" + f"{len(self.cfos_o)}\n" +
                              "cFOS (-) cells:\t" + f"{len(self.cfos_x)}\n" +
                              "\n")
            # self.summary = option_summary + "TOTAL cells:\t" + f"{cell_number}\n"
            self.cfos_text_widget.insert(tk.END, result_summary)
            self.cfos_text_widget.insert(tk.END, f"[cFOS (+) groups]\n{self.cfos_o['n'].sum}\n")
            self.cfos_text_widget.insert(tk.END, "\n")
            self.cfos_text_widget.insert(tk.END, f"[cFOS (-) groups]\n{self.cfos_x['n'].sum}\n")

            self.tester.update_cfos_data(self.cfos_o, self.cfos_x)
        else:
            self.cfos_text_widget.insert(tk.END, "cFOS grouping was cancelled or failed.")
   

if __name__ == "__main__":
    root = tk.Tk()
    root.title("FOS seq")
    root.geometry("900x600")

    # Frames
    root.rowconfigure(0, weight=1)
    root.rowconfigure(1, weight=1)
    root.columnconfigure(0, weight=0)
    root.columnconfigure(1, weight=1)
    root.columnconfigure(2, weight=1)

    collector_frame = tk.Frame(root, borderwidth=1, relief="solid", width=500)
    collector_frame.grid(row=0, column=0, rowspan=2, padx=10, pady=10, sticky='news')
    collector_frame.grid_propagate(False)

    tester_frame = tk.Frame(root, borderwidth=1, relief="solid")
    tester_frame.grid(row=0, column=1, padx=10, pady=10, sticky='news')

    sorter_option_frame = tk.Frame(root, borderwidth=1, relief="solid")
    sorter_option_frame.grid(row=1, column=1, padx=10, pady=10, sticky='news')

    sorter_result_frame = tk.Frame(root, borderwidth=1, relief="solid")
    sorter_result_frame.grid(row=0, column=2, rowspan=2, padx=10, pady=10, sticky='news')

    # Text widget
    collection_text_widget = tk.Text(tester_frame, wrap=tk.WORD)
    collection_text_widget.pack(padx=10, pady=10, expand=True, fill=tk.BOTH)

    # Classes
    sorter_app = Sorter(root, sorter_option_frame, sorter_result_frame)
    
    root.mainloop()