import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox

from imputer import Imputer


def main():
    root = tk.Tk()
    root.title("FOS seq")
    root.geometry("900x600")

    # Frames
    root.rowconfigure(0, weight=1)
    root.rowconfigure(1, weight=1)
    root.rowconfigure(2, weight=1)
    root.columnconfigure(0, weight=1)
    root.columnconfigure(1, weight=2)

    # ABC download path box
    abc_frame = tk.Frame(root, borderwidth=1, relief="solid")
    abc_frame.grid(row=0, column=0, sticky='news', padx=10, pady=10)

    # Training box
    train_frame = tk.Frame(root, borderwidth=1, relief="solid")
    train_frame.grid(row=1, column=0, sticky='news', padx=10, pady=10)

    # Testing box
    test_frame = tk.Frame(root, borderwidth=1, relief="solid")
    test_frame.grid(row=2, column=0, sticky='news', padx=10, pady=10)

    # Result box
    result_frame = tk.Frame(root, borderwidth=1, relief="solid")
    result_frame.grid(row=0, column=1, rowspan=3, sticky='news', padx=10, pady=10)


    # Classe
    imputer_app = Imputer(
        root=root,
        abc_frame=abc_frame,
        train_frame=train_frame,
        test_frame=test_frame,
        result_frame=result_frame
    )


    # # ABC path setting
    # file_abc_path = ''

    # def browse_folder(path_label, collector_app, sorter_app):
    #     nonlocal file_abc_path
    #     folder_selected = filedialog.askdirectory()
    #     if folder_selected:
    #         file_abc_path = folder_selected
    #         collector_app.set_abc_path(file_abc_path)
    #         sorter_app.set_abc_path(file_abc_path)
    #         path_label.config(text=folder_selected)

    # tk.Button(abc_frame, text="Select ABC Download Path(*)", command=lambda: browse_folder(path_label, collector_app, sorter_app)).pack(pady=10)
    # path_label = tk.Label(abc_frame, text="")
    # path_label.pack(pady=10)

    # # Text widget
    # collection_text_widget = tk.Text(tester_frame, wrap=tk.WORD)
    # collection_text_widget.pack(padx=10, pady=10, expand=True, fill=tk.BOTH)

    # cfos_text_widget = tk.Text(sorter_cfos_result_frame, wrap=tk.WORD)
    # cfos_text_widget.pack(padx=10, pady=10, expand=True, fill=tk.BOTH)


    root.mainloop()


if __name__ == "__main__":
    main()