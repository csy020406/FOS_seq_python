import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import threading
import queue
import pandas as pd
import numpy as np
import os
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

import merfish as mer


# Excute t-test
# Check data ready
class Tester:
    def __init__(self, root, frame, ttest_frame, text_widget):
        self.root = root
        self.frame = frame
        self.ttest_frame = ttest_frame
        self.text_widget = text_widget
        
        self.agg = None
        self.summary = None
        self.cfos_o = None
        self.cfos_x = None

        # Buttons
        self.button_frame = tk.Frame(self.frame, bd=0)
        self.button_frame.pack(padx=10, pady=(0, 10))

        self.load_button = ttk.Button(self.button_frame, text="Load", command=self.load_file)
        self.load_button.pack(side=tk.LEFT, padx=(0, 5))

        self.save_button = ttk.Button(self.button_frame, text="Save", command=self.save_file)
        self.save_button.pack(side=tk.RIGHT, padx=(5, 0))
        
        self.update_display()

        # T-test frame
        self.ready_agg = 0
        self.ready_cfos = 0
        self.ready_agg_label = tk.Label(self.ttest_frame, text="10X RNA seq Data Collection\tUNREADY", fg="red")
        self.ready_agg_label.pack(padx=10, pady=10)
        self.ready_cfos_label = tk.Label(self.ttest_frame, text="cFOS grouping\tUNREADY", fg="red")
        self.ready_cfos_label.pack(padx=10, pady=10)

        # Excute button
        self.ttest_button = ttk.Button(self.ttest_frame, text="Excute t-test", command=self.start_ttest)
        self.ttest_button.pack(pady=10)
        self.ttest_button.config(state=tk.DISABLED)

        # Matplotlib Figure
        self.fig = Figure(figsize=(10, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.ttest_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def check_ready(self):
        if self.ready_agg:
            self.ready_agg_label.config(text="10X RNA seq Data Collection\tREADY", fg="green")
        else:
            self.ready_agg_label.config(text="10X RNA seq Data Collection\tUNREADY", fg="red")

        if self.ready_cfos:
            self.ready_cfos_label.config(text="cFOS grouping\tREADY", fg="green")
        else:
            self.ready_cfos_label.config(text="cFOS grouping\tUNREADY", fg="red")

        if self.ready_agg and self.ready_cfos:
            self.ttest_button.config(state=tk.NORMAL)
        else:
            self.ttest_button.config(state=tk.DISABLED)

    # Called by Collector
    def update_agg(self, agg, summary):
        self.agg = agg
        self.summary = summary

        self.ready_agg = 1
        self.check_ready()
    
    def update_display(self):
        pd.set_option('display.multi_sparse', True)

        self.text_widget.delete(1.0, tk.END)
        if self.summary is not None:
            self.text_widget.insert(tk.END, f"{self.summary}\n")
        else:
            self.text_widget.insert(tk.END, "There is no loaded summary data.\n")

        if self.agg is not None:
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
                self.ready_agg = 0
                self.check_ready()
                # Check if proper data format
                if list(self.agg.columns.get_level_values(1)[:3]) != ['mean', 'std', 'count']:
                    messagebox.showerror("Error", "Please check data format.")
                    return
                elif not self.agg.index.to_series().apply(lambda x: isinstance(x, int)).all():
                    messagebox.showerror("Error", "Please check data format.")
                    return

                # Check txt file
                txt_file_path = file_path.replace(".csv", ".txt")
                if os.path.exists(txt_file_path):
                    with open(txt_file_path, 'r') as txt_file:
                        self.summary = txt_file.read()
                else:
                    self.summary = None

                self.ready_agg = 1
                self.check_ready()
                self.update_display()
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file:\n{e}")

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

    # Called by Sorter
    def update_cfos_data(self, cfos_o, cfos_x):
        self.cfos_o = cfos_o
        self.cfos_x = cfos_x

        self.ready_cfos = 1
        self.check_ready()

    def start_ttest(self):
        self.ttest_button.config(state=tk.DISABLED)
        self.show_progress_window()

        self.queue = queue.Queue()
        
        self.thread = threading.Thread(target=self.ttest)
        self.thread.start()
        
    def ttest(self):
        try:
            self.ttest_result = mer.t_test(self.agg, self.cfos_o, self.cfos_x, pc=self.update_progress)
        except Exception as e:
            messagebox.showerror("Error", str(e))
        finally:
            self.queue.put("Task Done")

    def show_progress_window(self):
        self.progress_window = tk.Toplevel(self.root)
        self.progress_window.title("Processing")
        self.progress_window.geometry("500x150")
        self.progress_window.resizable(False, False)

        self.progress_window.columnconfigure(0, weight=1)
        self.progress_window.columnconfigure(1, weight=1)
        self.progress_window.rowconfigure(1,weight=1)

        self.progress_label = tk.Label(self.progress_window, text="Doing t-test...", 
                                       font=("Arial", 10),
                                       anchor="w",
                                       justify="left")
        self.progress_label.grid(row=0, column=0, columnspan=2, padx=10, pady=(10,0), sticky="w")

        self.progress_var = tk.DoubleVar()

        self.progress_bar = ttk.Progressbar(self.progress_window, variable=self.progress_var, maximum=100)
        self.progress_bar.grid(row=1, column=0, columnspan=2, padx=10, pady=(0, 10), sticky="ew")

        self.cancel_button = ttk.Button(self.progress_window, text="Cancel", command=self.cancel_task)
        self.cancel_button.grid(row=2, column=1, padx=10, pady=10, sticky="se")

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
            self.ttest_button.config(state=tk.NORMAL)

            # Show completion message
            self.show_results()

    def cancel_task(self):
        # Signal the background task to cancel
        self.queue.put("Cancel")

        # Stop the progress bar and close the window
        self.progress_bar.stop()
        self.progress_window.destroy()
        self.ttest_button.config(state=tk.NORMAL)

    # 'current' represents the current gene
    # 'total' represents total genes
    def update_progress(self, current, total):
        self._gradual_update_progress(self.progress_var.get(), (current / total) * 100)

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

    def show_results(self):
        self.draw_volcano_plot()

    def draw_volcano_plot(self):
        self.ax.clear()

        fold_changes = self.ttest_result['fold_change']
        p_values = self.ttest_result['p_value']

        neg_log_p_values = p_values

        self.ax.scatter(fold_changes, neg_log_p_values, alpha=0.75)
        
        for gene, fold_change, neg_log_p_value in zip(self.ttest_result.index, fold_changes, neg_log_p_values):
            if abs(fold_change) > 1 or neg_log_p_value > 2:
                self.ax.text(fold_change, neg_log_p_value, gene, fontsize=8)
        
        self.ax.set_xlabel('Fold Change')
        self.ax.set_ylabel('-log10(p-value)')
        self.ax.set_title('t-test')
        self.ax.axhline(y=-np.log10(0.05), color='r', linestyle='--')
        self.ax.axvline(x=1, color='r', linestyle='--')
        self.ax.axvline(x=-1, color='r', linestyle='--')

        self.canvas.draw()