import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import os


# Edit agg and deliver to Tester
class Editor:
    def __init__(self, root, frame):
        self.root = root
        self.agg = None
        self.summary = None
        self.frame = frame
        
        self.label = ttk.Label(self.frame, text="Current Agg Data")
        self.label.pack(padx=10, pady=10)
        
        self.text_widget = tk.Text(self.frame, wrap=tk.WORD)
        self.text_widget.pack(padx=10, pady=(10,0), expand=True, fill=tk.BOTH)
        
        # Buttons
        self.button_frame = tk.Frame(self.frame, bd=0)
        self.button_frame.pack(padx=10, pady=(0, 10))

        self.load_button = ttk.Button(self.button_frame, text="Load", command=self.load_file)
        self.load_button.pack(side=tk.LEFT, padx=(0, 5))

        self.save_button = ttk.Button(self.button_frame, text="Save", command=self.save_file)
        self.save_button.pack(side=tk.RIGHT, padx=(5, 0))
        
        self.update_display()

    # Called by Collector
    def update_agg(self, agg, summary):
        self.agg = agg
        self.summary = summary
        self.update_display()
    
    def update_display(self):
        self.text_widget.delete(1.0, tk.END)
        if self.agg is not None:
            self.text_widget.insert(tk.END, f"{self.summary}\n")
            self.text_widget.insert(tk.END, f"{self.agg}\n")
        else:
            self.text_widget.insert(tk.END, "There is no loaded gene expression data.\n")
    
    def load_file(self):
        file_path = filedialog.askopenfilename(
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if file_path:
            try:
                self.agg = pd.read_csv(file_path, index_col=0, header=[0, 1])
                # Check if proper data format
                if list(self.agg.columns.get_level_values(1)[:3]) != ['mean', 'std', 'count']:
                    messagebox.showerror("Error", "Please check data format.")
                    return
                elif not self.agg.index.to_series().apply(lambda x: isinstance(x, int)).all():
                    messagebox.showerror("Error", "Please check data format.")
                    return

                self.text_widget.delete(1.0, tk.END)
                self.text_widget.insert(tk.END, f"{self.agg}\n")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to laod file:\n{e}")

    def save_file(self):
        if self.agg is None:
            messagebox.showwarning("Save", "No data to save!")
            return

        self.agg.reset_index(inplace=True)

        file_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if file_path:
            try:
                self.agg.to_csv(file_path, index=False)
                if self.summary:
                    txt_file_path = os.path.splitext(file_path)[0] + ".txt"
                    with open(txt_file_path, 'w') as txt_file:
                        txt_file.write(self.summary)

            except Exception as e:
                messagebox.showerror("Error", f"Failed to save file:\n{e}")