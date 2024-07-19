import tkinter as tk
from tkinter import ttk

import rna_seq_data_collect as rsd
import merfish as mer

from collector import Collector
from editor import Editor
from tester import Tester

 
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

    collector_result_frame = tk.Frame(root, borderwidth=1, relief="solid")
    collector_result_frame.grid(row=0, column=1, padx=10, pady=10, sticky='news')

    editor_frame = tk.Frame(root, borderwidth=1, relief="solid")
    editor_frame.grid(row=1, column=1, padx=10, pady=10, sticky='news')

    tester_frame = tk.Frame(root, borderwidth=1, relief="solid")
    tester_frame.grid(row=0, column=2, rowspan=2, padx=10, pady=10, sticky='news')

    # Classes
    editor_app = Editor(root, editor_frame)
    collector_app = Collector(root, collector_frame, collector_result_frame, editor_app)
    tester_app = Tester(root, tester_frame)
    
    root.mainloop()