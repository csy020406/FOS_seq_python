import tkinter as tk
from tkinter import ttk

import rna_seq_data_collect as rsd
import merfish as mer


class Tester:
    def __init__(self, root, frame):
        self.root = root
        self.frame = frame

        self.label = tk.Label(self.frame, text="t-test")
        self.label.pack(padx=10, pady=10)
  