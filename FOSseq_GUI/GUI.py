import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox

import rna_seq_data_collect as rsd
import merfish as mer

from collector import Collector
from tester import Tester
from sorter import Sorter


def main():
    root = tk.Tk()
    root.title("FOS seq")
    root.geometry("900x600")

    # Frames
    root.rowconfigure(0, weight=1)
    root.rowconfigure(1, weight=1)
    root.columnconfigure(0, weight=0)
    root.columnconfigure(1, weight=1)
    root.columnconfigure(2, weight=1)
    root.columnconfigure(3, weight=2)

    # Fixed left region
    fixed_left_frame = tk.Frame(root, bd=0, width=300)
    fixed_left_frame.grid(row=0, column=0, rowspan=2, sticky='news')
    fixed_left_frame.grid_propagate(False)

    # ABC address setting region
    abc_frame = tk.Frame(fixed_left_frame, borderwidth=1, relief="solid", height=200)
    abc_frame.pack(fill=tk.X, padx=10, pady=10)
    abc_frame.pack_propagate(False)

    # Collector region
    collector_frame = tk.Frame(fixed_left_frame, borderwidth=1, relief="solid")
    collector_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

    # Tester region
    tester_frame = tk.Frame(root, borderwidth=1, relief="solid")
    tester_frame.grid(row=0, column=1, columnspan=2, padx=10, pady=10, sticky='news')

    # Sorter region
    sorter_frame = tk.Frame(root, borderwidth=1, relief="solid")
    sorter_frame.grid(row=1, column=1, columnspan=2, padx=10, pady=10, sticky='news')
    
    sorter_frame.columnconfigure(0, weight=0)
    sorter_frame.columnconfigure(1, weight=1)

    sorter_option_frame = tk.Frame(sorter_frame, bd=0, width=500)
    sorter_option_frame.grid(row=0, column=0, sticky='news')
    sorter_option_frame.grid_propagate(False)

    sorter_cfos_result_frame = tk.Frame(sorter_frame, bd=0)
    sorter_cfos_result_frame.grid(row=0, column=1, sticky='news')

    ttest_result_frame = tk.Frame(root, borderwidth=1, relief="solid")
    ttest_result_frame.grid(row=0, column=3, rowspan=2, padx=10, pady=10, sticky='news')

    # ABC path setting
    file_abc_path = ''

    def browse_folder(path_label, collector_app, sorter_app):
        nonlocal file_abc_path
        folder_selected = filedialog.askdirectory()
        if folder_selected:
            file_abc_path = folder_selected
            collector_app.set_abc_path(file_abc_path)
            sorter_app.set_abc_path(file_abc_path)
            path_label.config(text=folder_selected)

    tk.Button(abc_frame, text="Select ABC Download Path(*)", command=lambda: browse_folder(path_label, collector_app, sorter_app)).pack(pady=10)
    path_label = tk.Label(abc_frame, text="")
    path_label.pack(pady=10)

    # Text widget
    collection_text_widget = tk.Text(tester_frame, wrap=tk.WORD)
    collection_text_widget.pack(padx=10, pady=10, expand=True, fill=tk.BOTH)

    cfos_text_widget = tk.Text(sorter_cfos_result_frame, wrap=tk.WORD)
    cfos_text_widget.pack(padx=10, pady=10, expand=True, fill=tk.BOTH)

    # Classes
    tester_app = Tester(root, tester_frame, ttest_result_frame, collection_text_widget)
    sorter_app = Sorter(root, sorter_option_frame, sorter_cfos_result_frame, cfos_text_widget, tester_app)
    collector_app = Collector(root, collector_frame, tester_frame, collection_text_widget, tester_app)

    root.mainloop()


if __name__ == "__main__":
    main()