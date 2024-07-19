import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk
import threading
import queue

from editor import Editor

import rna_seq_data_collect as rsd
import merfish as mer

class Collector:
    def __init__(self, root, frame, result_frame, editor):
        self.root = root
        self.editor = editor
        
        self.frame = frame
        self.label = ttk.Label(self.frame, text="10X RNA seq Data Collection")
        self.label.pack(padx=10, pady=10)

        self.result_frame = result_frame
        self.result_label = ttk.Label(self.result_frame, text="Data Collection Results")
        self.result_label.pack(padx=10, pady=10)

        self.text_widget = tk.Text(self.result_frame, wrap=tk.WORD)
        self.text_widget.pack(padx=10, pady=10, expand=True, fill=tk.BOTH)

        self.update_button = ttk.Button(self.result_frame, text="Update", command=lambda: self.editor.update_agg(self.agg, self.summary))
        self.update_button.pack(pady=10)
        self.update_button.config(state=tk.DISABLED)

        # User File Path
        self.file_path = tk.StringVar()
        tk.Button(self.frame, text="Select Download Path", command=self.browse_folder).pack(pady=10)
        self.path_label = tk.Label(self.frame, textvariable=self.file_path)
        self.path_label.pack(pady=10)

        # OPTION 1: 10Xv2, 10Xv3, 10XMulti
        self.option1_vars = [tk.IntVar() for _ in range(3)]
        option1 = ["WMB-10Xv2", "WMB-10Xv3", "WMB-10XMulti"]
        tk.Label(self.frame, text="Select Directories (multiple choice):").pack()
        for i, option in enumerate(option1):
            tk.Checkbutton(self.frame, text=option, variable=self.option1_vars[i], command=self.show_datafiles).pack(anchor='w')

        # OPTION 2: file names
        tk.Label(self.frame, text="Select Files:").pack()
        self.option2_listbox = tk.Listbox(self.frame, selectmode=tk.MULTIPLE, exportselection=False)
        self.option2_listbox.pack(fill=tk.BOTH, expand=True)
        self.scrollbar = tk.Scrollbar(self.option2_listbox, orient=tk.VERTICAL)
        self.scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.option2_listbox.config(yscrollcommand=self.scrollbar.set)
        self.scrollbar.config(command=self.option2_listbox.yview)

        # OPTION 3: data type
        self.option3_var = tk.IntVar()
        option3 = [("log2", 0), ("raw", 1)]
        tk.Label(self.frame, text="Select Data Type:").pack()
        for option, val in option3:
            tk.Radiobutton(self.frame, text=option, variable=self.option3_var, value=val).pack(anchor='w')
        
        # OPTION 4: gene range
        self.range_start = tk.IntVar()
        self.range_end = tk.IntVar()
        tk.Label(self.frame, text="Enter Range (0-32284)\n(-1, -1) for all genes, but not recommended:").pack()
        tk.Entry(self.frame, textvariable=self.range_start).pack()
        tk.Entry(self.frame, textvariable=self.range_end).pack()

        # Excute button
        self.exe_button = ttk.Button(self.frame, text="Excute Data Collecting", command=self.start_data_collect)
        self.exe_button.pack(pady=10)

    def browse_folder(self):
        folder_selected = filedialog.askdirectory()
        if folder_selected:
            self.file_path.set(folder_selected)

    def show_datafiles(self):
        # Get option 1
        user_option1 = [opt for var, opt in zip(self.option1_vars, ["WMB-10Xv2", "WMB-10Xv3", "WMB-10XMulti"]) if var.get() == 1]
        if not user_option1:
            return
        
        datafiles = rsd.list_file_name(user_option1)

        # Update listbox
        self.option2_listbox.delete(0, tk.END)
        for file_name in datafiles:
            self.option2_listbox.insert(tk.END, file_name)

    def show_progress_window(self):
        self.progress_window = tk.Toplevel(self.root)
        self.progress_window.title("Processing")
        self.progress_window.geometry("500x150")
        self.progress_window.resizable(False, False)

        self.progress_window.columnconfigure(0, weight=1)
        self.progress_window.columnconfigure(1, weight=1)
        self.progress_window.rowconfigure(1,weight=1)
        
        self.progress_label = tk.Label(self.progress_window, text="10X RNA seq Data Collection Excuted.\nPreparing...", 
                                       font=("Arial", 10),
                                       anchor="w",
                                       justify="left")
        self.progress_label.grid(row=0, column=0, columnspan=2, padx=20, pady=(10,0), sticky="w")

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
            self.exe_button.config(state=tk.NORMAL)
            self.update_button.config(state=tk.NORMAL)

            self.show_results()

    def cancel_task(self):
        self.queue.put("Cancel")

        self.progress_bar.stop()
        self.progress_window.destroy()
        self.exe_button.config(state=tk.NORMAL)

    # Excute data_collect
    def start_data_collect(self):
        # ==== Reflect USER OPTIONs ====
        # User File Path
        path = self.file_path.get()
        if not path:
            messagebox.showerror("Error", "Please select a download path.")
            return
        
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
        
        # OPTION 4
        start = self.range_start.get()
        end = self.range_end.get()
        if (start == -1 and end == -1):
            start = 0
            end = 33000 # for mis-calculation
        if not (0 <= start <= 32284 and 0 <= end <= 33000 and start <= end):
            messagebox.showerror("Error", "Please enter valid integers in the range 1-32285.")
            return

        # ==== All Options are ready ====
        self.path = path
        self.user_option1 = user_option1
        self.user_option2 = user_option2
        self.user_option3 = user_option3
        self.start = start
        self.end = end

        self.exe_button.config(state=tk.DISABLED)
        self.update_button.config(state=tk.DISABLED)
        self.show_progress_window()

        self.queue = queue.Queue()
        
        self.thread = threading.Thread(target=self.data_collect)
        self.thread.start()


    # Excute rsd.data_collect
    def data_collect(self):
        # current represents current task or current gene
        # total represents whethere current task done or total genes
        def update_progress(current, total):
            if current == -1:
                if total == 0:
                    self.progress_label.config(text="10X RNA seq Data Collection Excuted.\nGenerating Cell Metadata ...")
                else:
                    self.progress_var.set(10)
                    self.progress_window.update_idletasks()     # Update GUI
            else:
                self.progress_label.config(text="10X RNA seq Data Collection Excuted.\nAggregating Gene Expression data ...")
                progress = 10 + (current / total) * 90
                self.progress_var.set(progress)
                self.progress_window.update_idletasks()     # Update GUI

        try:
            rsd.change_address(self.path)
            mer.change_address(self.path)

            # Reflect USER options
            self.agg = rsd.data_collect(self.user_option1, self.user_option2, self.user_option3, self.start, self.end, pc=update_progress)

        except Exception as e:
            messagebox.showerror("Error", str(e))
        finally:
            self.queue.put("Task Done")

    def show_results(self):
        self.text_widget.delete("1.0", "end")

        option_summary = ("========== USER OPTION ==========\n" +
                            "USER path:\t" + f"{self.path}\n" +
                            "USER directories:\t" + f"{self.user_option1}\n" +
                            "USER files:\t" + f"{self.user_option2}\n" +
                            "USER data type option:\t" + f"{'raw' if self.user_option3 else 'log2'}\n" +
                            "Gene range:\t" + f"({self.start}-{self.end})\n" +
                            "\n")
        self.text_widget.insert(tk.END, option_summary)

        cell_number = rsd.get_cell_number()

        if hasattr(self, 'agg'):
            result_summary = ("========= RESULT SUMMARY =========\n" +
                                "TOTAL cells:\t" + f"{cell_number}\n")
            self.summary = option_summary + "TOTAL cells:\t" + f"{cell_number}\n"
            self.text_widget.insert(tk.END, result_summary)
            self.text_widget.insert(tk.END, f"{self.agg}\n")
        else:
            self.text_widget.insert(tk.END, "Data collection was cancelled or failed.")
    
