import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk

import rna_seq_data_collect as rsd

class Collector:
    def __init__(self, root):
        self.root = root
        self.frame = tk.Frame(self.root, borderwidth=1, relief="solid")
        self.frame.pack(side="left", padx=10, pady=10, fill="both", expand=True)

        self.label = ttk.Label(self.frame, text="10X RNA seq Data Collection")
        self.label.pack(padx=10, pady=10)

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
        tk.Button(self.frame, text="Excute Data Collecting", command=self.data_collect).pack(pady=10)

    def browse_folder(self):
        folder_selected = filedialog.askdirectory()
        if folder_selected:
            self.file_path.set(folder_selected)

    def show_datafiles(self):
        # Get option 1
        selected_option1 = [opt for var, opt in zip(self.option1_vars, ["WMB-10Xv2", "WMB-10Xv3", "WMB-10XMulti"]) if var.get() == 1]
        if not selected_option1:
            return
        
        datafiles = rsd.list_file_name(selected_option1)

        # Update listbox
        self.option2_listbox.delete(0, tk.END)
        for file_name in datafiles:
            self.option2_listbox.insert(tk.END, file_name)

    def data_collect(self):
        # User File Path
        path = self.file_path.get()
        if not path:
            messagebox.showerror("Error", "Please select a download path.")
            return
        
        # OPTION 1
        selected_option1 = [opt for var, opt in zip(self.option1_vars, ["WMB-10Xv2", "WMB-10Xv3", "WMB-10XMulti"]) if var.get() == 1]
        if not selected_option1:
            messagebox.showerror("Error", "Please select at least one option in the set.")
            return
        
        # OPTION 2
        selected_option2 = [self.option2_listbox.get(i) for i in self.option2_listbox.curselection()]
        if not selected_option2:
            messagebox.showerror("Error", "Please select at least one file in the set.")
            return
        
        # OPTION 3
        selected_option3 = self.option3_var.get()   # 0 for "log2", 1 for "raw"
        
        # OPTION 4
        start = self.range_start.get()
        end = self.range_end.get()
        if (start == -1 and end == -1):
            start = 0
            end = 33000 # for mis-calculation
        if not (0 <= start <= 32284 and 0 <= end <= 33000 and start <= end):
            messagebox.showerror("Error", "Please enter valid integers in the range 1-32285.")
            return
        
        # import function from rsd
        try:
            print(path)
            print(selected_option1)
            print(selected_option2)
            print(selected_option3)
            print(start)
            print(end)

            rsd.change_address(path)
            # Reflect USER options
            if len(selected_option2) == 1:
                rsd.change_feature_matrix_label(selected_option1[0], selected_option2[0], selected_option3)
            else:
                print("not yet")
            rsd.change_gene_range(start, end)
            rsd.data_collect()
        except Exception as e:
            messagebox.showerror("Error", str(e))

class Tester:
    def __init__(self, root):
        self.root = root
        self.frame = tk.Frame(self.root, borderwidth=1, relief="solid")
        self.frame.pack(side="right", padx=10, pady=10, fill="both", expand=True)

        self.label = tk.Label(self.frame, text="t-test")
        self.label.pack(padx=10, pady=10)
        
if __name__ == "__main__":
    root = tk.Tk()
    root.title("FOS seq")
    root.geometry("640x480")

    collector_app = Collector(root)
    tester_app = Tester(root)

    root.mainloop()