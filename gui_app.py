from __future__ import annotations

import re
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from analysis_engine import compute_analysis
from export_tab import ExportTab
from fasta_parser import load_fasta
from iupac import count_matches
from ncbi_client import fetch_fasta_by_ids, is_accession_like


class DNAApp(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self._configure_styles()
        self._init_state()
        self.create_menu()
        self.create_layout()
        if self.DEV_MODE:
            self.load_test_data()

    def _configure_styles(self) -> None:
        self.style = ttk.Style(self)
        self.style.theme_use("clam")

        self.UI_FONT = ("Segoe UI", 10)
        self.UI_H1 = ("Segoe UI", 12, "bold")
        self.ACCENT = "#0b5394"

        self.style.configure(".", font=self.UI_FONT)
        self.style.configure("TNotebook.Tab", padding=(12, 6))
        self.style.configure("TLabel", font=self.UI_FONT)
        self.style.configure("Header.TLabel", font=self.UI_H1, foreground=self.ACCENT)
        self.style.configure("Card.TFrame", background="#f7f7f7", borderwidth=1, relief="solid")
        self.style.configure("CardInner.TFrame", background="#f7f7f7")

        self.title("Projekt 1: Analiza motywów sekwencyjnych w DNA")
        self.geometry("1000x600")

    def _init_state(self) -> None:
        self.sequences: dict[str, str] = {}
        self.current_fasta: str | None = None
        self.motifs: list[str] = []
        self.analysis_result = None
        self.figures: dict[str, object] = {}

        self.DEV_MODE = False
        self.analysis_done = False

        self.sort_state: dict[str, bool] = {}
        self.normalization_mode = tk.StringVar(value="raw")

        self.active_label_x = None
        self.active_label_y = None
        self.selected_sequence: str | None = None
        self.selected_motif: str | None = None
        self.current_canvas = None
        self.hover_annotation = None
        self.heatmap_data = None
        self.help_win = None

    def create_menu(self) -> None:
        menubar = tk.Menu(self)

        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Nowa analiza", command=self.new_analysis)
        file_menu.add_separator()
        file_menu.add_command(label="Wczytaj plik", command=self.gui_load_file)
        file_menu.add_separator()
        file_menu.add_command(label="Wyjście", command=self.quit)

        motif_menu = tk.Menu(menubar, tearoff=0)
        motif_menu.add_command(label="Dodaj/usuń motywy", command=self.open_motif_manager)

        ncbi_menu = tk.Menu(menubar, tearoff=0)
        ncbi_menu.add_command(label="Pobierz z NCBI", command=self.download_ncbi)

        menubar.add_cascade(label="Plik", menu=file_menu)
        menubar.add_cascade(label="Motywy", menu=motif_menu)
        menubar.add_cascade(label="NCBI", menu=ncbi_menu)
        menubar.add_command(label="Pomoc", command=self.open_help_window)

        self.config(menu=menubar)
        self.bind_all("<F1>", lambda event: self.open_help_window())

    def create_layout(self) -> None:
        for widget in self.winfo_children():
            widget.pack_forget()
            widget.grid_forget()

        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=0)
        self.grid_columnconfigure(0, weight=1)

        main_frame = tk.Frame(self)
        main_frame.grid(row=0, column=0, sticky="nsew")

        bottom_frame = tk.Frame(self, height=140)
        bottom_frame.grid(row=1, column=0, sticky="ew")
        bottom_frame.grid_propagate(False)

        right = tk.Frame(main_frame)
        right.pack(fill="both", expand=True)

        self.tabs = ttk.Notebook(right)
        self.tabs.pack(fill="both", expand=True)

        self.tab_preview = ttk.Frame(self.tabs)
        self.tab_results = ttk.Frame(self.tabs)
        self.tab_viz = ttk.Frame(self.tabs)
        self.tab_export = ttk.Frame(self.tabs)

        self.tabs.add(self.tab_preview, text="Podgląd sekwencji")
        self.tabs.add(self.tab_results, text="Wyniki analizy")
        self.tabs.add(self.tab_viz, text="Wizualizacja")
        self.tabs.add(self.tab_export, text="Eksport")

        export_widget = ExportTab(
            self.tab_export,
            get_result=lambda: self.analysis_result,
            get_figures=lambda: self.figures,
            log_fn=self.log,
        )
        export_widget.pack(fill="both", expand=True)

        self.preview_box = tk.Text(self.tab_preview)
        self.preview_box.pack(fill="both", expand=True)

        self._build_results_tab()
        self._build_visualization_tab()
        self._build_status_and_logs(bottom_frame)
        self.update_status()

    def _build_results_tab(self) -> None:
        results_frame = tk.Frame(self.tab_results)
        results_frame.pack(fill="both", expand=True)

        mode_frame = tk.Frame(results_frame)
        mode_frame.pack(fill="x", pady=(5, 0))

        tk.Label(mode_frame, text="Tryb wyświetlania:").pack(side="left", padx=(5, 5))
        tk.Radiobutton(mode_frame, text="Surowe liczby", variable=self.normalization_mode, value="raw", command=self.run_analysis).pack(side="left")
        tk.Radiobutton(mode_frame, text="Na 1000 nt", variable=self.normalization_mode, value="norm", command=self.run_analysis).pack(side="left")

        self.results_table = ttk.Treeview(results_frame, show="headings")
        self.results_table.pack(fill="both", expand=True)
        self.results_table.bind("<Motion>", self._on_results_table_motion)

    def _build_visualization_tab(self) -> None:
        self.viz_frame = tk.Frame(self.tab_viz)
        self.viz_frame.pack(fill="both", expand=True)

        viz_controls = tk.Frame(self.viz_frame)
        viz_controls.pack(fill="x", pady=(5, 0))
        tk.Label(viz_controls, text="Tryb wyświetlania:").pack(side="left", padx=(10, 5))
        tk.Radiobutton(viz_controls, text="Surowe", variable=self.normalization_mode, value="raw", command=self.refresh_visualization).pack(side="left")
        tk.Radiobutton(viz_controls, text="Na 1000 nt", variable=self.normalization_mode, value="norm", command=self.refresh_visualization).pack(side="left")

        self.viz_top = tk.Frame(self.viz_frame)
        self.viz_top.pack(fill="both", expand=True)

        self.viz_bottom = tk.Frame(self.viz_frame, height=260)
        self.viz_bottom.pack(side="bottom", fill="x")
        self.viz_bottom.pack_propagate(False)

        self.close_barplot_btn = tk.Button(
            self.viz_bottom,
            text="✖",
            command=self.hide_barplot,
            bd=0,
            highlightthickness=0,
            font=("Segoe UI", 11, "bold"),
            cursor="hand2",
        )
        self.close_barplot_btn.place_forget()

    def _build_status_and_logs(self, parent: tk.Frame) -> None:
        self.status_var = tk.StringVar(value="Gotowe")
        status_bar = tk.Label(parent, textvariable=self.status_var, bd=1, relief="sunken", anchor="w")
        status_bar.pack(side="bottom", fill="x")

        log_frame = tk.Frame(parent)
        log_frame.pack(side="top", fill="both", expand=True)
        tk.Label(log_frame, text="Logi / komunikaty:").pack(anchor="w", padx=6)

        self.log_box = tk.Text(log_frame, height=3)
        self.log_box.pack(fill="both", expand=True, padx=6, pady=(0, 6))

    def _on_results_table_motion(self, event) -> None:
        region = self.results_table.identify_region(event.x, event.y)
        self.results_table.config(cursor="hand2" if region == "heading" else "")

    def log(self, message: str) -> None:
        self.log_box.insert("end", message + "\n")
        self.log_box.see("end")

    def update_status(self) -> None:
        mode = "Surowe" if self.normalization_mode.get() == "raw" else "Na 1000 nt"
        self.status_var.set(f"✔ Sekwencje: {len(self.sequences)} | Motywy: {len(self.motifs)} | Tryb: {mode}")

    def load_test_data(self) -> None:
        self.sequences = {
            "seq1": "ATGCGTATGCGTATGTATAAAGGGCGGATGCGT",
            "seq2": "TTGACATATAATAGGAGGATGCGTTAA",
            "seq10": "GGGAAATGCAAATAATAAATGCGTATG",
            "seq20": "TATAAAGGGCGGCANNTGTGACGTCA",
        }
        self.motifs = ["ATG", "TATAAA", "GGGCGG"]

        self.preview_box.delete("1.0", "end")
        for sequence_id, sequence in self.sequences.items():
            self.preview_box.insert("end", f">{sequence_id}\n{sequence}\n\n")

        self.run_analysis()
        self.log("Załadowano dane testowe (DEV_MODE)")
        self.update_status()

    def _confirm_replace_current_analysis(self, parent=None, source_label: str = "nowy plik") -> bool:
        if not self.analysis_done:
            return True
        return messagebox.askyesno(
            "Nowa analiza",
            f"Masz już wykonaną analizę.\nCzy chcesz wczytać {source_label} i usunąć poprzednie wyniki oraz motywy?",
            parent=parent,
        )

    def gui_load_file(self) -> None:
        if not self._confirm_replace_current_analysis(source_label="nowy plik"):
            return
        if self.analysis_done:
            self.clear_results()

        path = filedialog.askopenfilename(
            title="Wybierz FASTA",
            filetypes=[("FASTA", "*.fasta *.fa *.fna"), ("All", "*.*")],
        )
        if not path:
            return

        try:
            self.load_sequences_from_file(path)
            self.log(f"Wczytano plik: {path}")
            self.update_status()
        except Exception as exc:
            messagebox.showerror("Błąd FASTA", str(exc))

    def load_sequences_from_file(self, path: str) -> None:
        self.sequences = load_fasta(path)
        self.current_fasta = path
        self.analysis_result = None
        self.figures.clear()
        self.analysis_done = False
        self._render_preview()
        self._clear_results_views_only()

    def _render_preview(self) -> None:
        self.preview_box.delete("1.0", "end")
        for idx, (sequence_id, sequence) in enumerate(self.sequences.items()):
            if idx == 20:
                self.preview_box.insert("end", "\n... (reszta ukryta)\n")
                break
            suffix = "..." if len(sequence) > 200 else ""
            self.preview_box.insert("end", f">{sequence_id}\n{sequence[:200]}{suffix}\n\n")

    def open_motif_manager(self) -> None:
        IUPAC_ALLOWED = set("ACGTURYSWKMBDHVN")

        def normalize_motif(text: str) -> str:
            return "".join(text.split()).upper()

        def validate_motif(motif: str) -> tuple[bool, str]:
            if not motif:
                return False, "Motyw jest pusty."
            bad = sorted({char for char in motif if char not in IUPAC_ALLOWED})
            if bad:
                allowed = "".join(sorted(IUPAC_ALLOWED))
                return False, f"Nieprawidłowe znaki: {', '.join(bad)}.\nDozwolone: {allowed}"
            if len(motif) < 2:
                return False, "Motyw jest za krótki (min. 2 znaki)."
            return True, ""

        def unique_preserve_order(items: list[str]) -> list[str]:
            seen = set()
            ordered: list[str] = []
            for item in items:
                if item not in seen:
                    seen.add(item)
                    ordered.append(item)
            return ordered

        def plural_motif(count: int) -> str:
            if count == 1:
                return "motyw"
            if 12 <= count % 100 <= 14:
                return "motywów"
            if 2 <= count % 10 <= 4:
                return "motywy"
            return "motywów"

        win = tk.Toplevel(self)
        win.title("Motywy")
        win.transient(self)
        win.grab_set()
        win.resizable(False, True)

        initial_width = 460
        win.geometry(f"{initial_width}x520")

        def autosize() -> None:
            win.update_idletasks()
            req_h = win.winfo_reqheight()
            max_h = max(360, win.winfo_screenheight() - 120)
            win.geometry(f"{initial_width}x{min(req_h, max_h)}")

        root = ttk.Frame(win, padding=10)
        root.pack(fill="both", expand=True)
        root.columnconfigure(0, weight=1)
        root.rowconfigure(5, weight=1)

        top_row = ttk.Frame(root)
        top_row.grid(row=0, column=0, sticky="ew", pady=(0, 2))
        top_row.columnconfigure(0, weight=1)
        ttk.Label(top_row, text="Dozwolone symbole:", font=("Segoe UI", 9, "bold")).grid(row=0, column=0, sticky="w")

        def run_and_close() -> None:
            win.destroy()
            self.tabs.select(self.tab_results)
            self.after(50, self.run_analysis)

        btn_run = ttk.Button(top_row, text="Uruchom analizę", command=run_and_close)
        btn_run.grid(row=0, column=1, sticky="e")

        iupac_lbl = ttk.Label(root, text="A C G T + IUPAC (R Y S W K M B D H V N)", foreground="#0b5394")
        iupac_lbl.grid(row=1, column=0, sticky="w", pady=(0, 6))

        tooltip_text = (
            "R = A lub G\nY = C lub T\nS = G lub C\nW = A lub T\nK = G lub T\nM = A lub C\n"
            "B = C G T\nD = A G T\nH = A C T\nV = A C G\nN = dowolna baza"
        )
        tooltip = None

        def show_tooltip(event) -> None:
            nonlocal tooltip
            if tooltip:
                return
            tooltip = tk.Toplevel(win)
            tooltip.wm_overrideredirect(True)
            tooltip.geometry(f"+{event.x_root + 12}+{event.y_root + 12}")
            tk.Label(tooltip, text=tooltip_text, justify="left", bg="#ffffe0", relief="solid", borderwidth=1, padx=6, pady=4).pack()

        def hide_tooltip(_event) -> None:
            nonlocal tooltip
            if tooltip:
                tooltip.destroy()
                tooltip = None

        iupac_lbl.bind("<Enter>", show_tooltip)
        iupac_lbl.bind("<Leave>", hide_tooltip)

        ttk.Label(root, text="Dodaj motyw:", font=("Segoe UI", 9, "bold")).grid(row=2, column=0, sticky="w", pady=(4, 2))
        add_row = ttk.Frame(root)
        add_row.grid(row=3, column=0, sticky="ew")
        add_row.columnconfigure(0, weight=1)
        motif_entry = ttk.Entry(add_row)
        motif_entry.grid(row=0, column=0, sticky="ew")

        ttk.Label(root, text="Aktualne motywy:", font=("Segoe UI", 9, "bold")).grid(row=4, column=0, sticky="w", pady=(10, 2))
        list_frame = ttk.Frame(root)
        list_frame.grid(row=5, column=0, sticky="nsew")
        list_frame.columnconfigure(0, weight=1)
        list_frame.rowconfigure(0, weight=1)

        motifs_lb = tk.Listbox(list_frame, selectmode="extended", height=8, activestyle="none")
        motifs_lb.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(list_frame, orient="vertical", command=motifs_lb.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        motifs_lb.configure(yscrollcommand=scroll.set)

        actions = ttk.Frame(root)
        actions.grid(row=6, column=0, sticky="ew", pady=(10, 6))
        for col in range(4):
            actions.columnconfigure(col, weight=1)

        badge_var = tk.StringVar(value="0 motywów")
        ttk.Label(actions, textvariable=badge_var, anchor="e").grid(row=1, column=0, columnspan=4, sticky="e", pady=(6, 0))

        def invalidate_analysis() -> None:
            self.analysis_result = None
            self.figures.clear()
            self.analysis_done = False
            self._clear_results_views_only()

        def refresh_lists() -> None:
            motifs_lb.delete(0, "end")
            for motif in self.motifs:
                motifs_lb.insert("end", motif)
            count = len(self.motifs)
            badge_var.set(f"{count} {plural_motif(count)}")
            has_any = bool(self.motifs)
            btn_export.configure(state=("normal" if has_any else "disabled"))
            btn_clear.configure(state=("normal" if has_any else "disabled"))
            btn_remove.configure(state=("normal" if motifs_lb.curselection() else "disabled"))
            btn_run.configure(state=("normal" if self.sequences and self.motifs else "disabled"))
            autosize()

        def add_motif(raw_text: str) -> None:
            motif = normalize_motif(raw_text)
            ok, reason = validate_motif(motif)
            if not ok:
                messagebox.showwarning("Nieprawidłowy motyw", reason, parent=win)
                return
            if motif in self.motifs:
                self.log(f"Motyw już istnieje: {motif}")
                motif_entry.delete(0, "end")
                return
            self.motifs.append(motif)
            self.motifs = unique_preserve_order(self.motifs)
            invalidate_analysis()
            motif_entry.delete(0, "end")
            self.log(f"Dodano motyw: {motif}")
            refresh_lists()
            refresh_accordion_checks()
            self.update_status()

        def remove_selected() -> None:
            indices = sorted(motifs_lb.curselection(), reverse=True)
            if not indices:
                return
            removed = [self.motifs.pop(index) for index in indices]
            self.log(f"Usunięto: {', '.join(removed)}")
            invalidate_analysis()
            refresh_lists()
            refresh_accordion_checks()
            self.update_status()

        def clear_all() -> None:
            if not self.motifs:
                return
            if not messagebox.askyesno("Wyczyścić motywy?", "Na pewno usunąć wszystkie motywy?", parent=win):
                return
            self.motifs.clear()
            invalidate_analysis()
            self.log("Wyczyszczono wszystkie motywy.")
            refresh_lists()
            refresh_accordion_checks()
            self.update_status()

        def add_imported_tokens(tokens: list[str]) -> None:
            added = 0
            for token in tokens:
                motif = normalize_motif(token)
                ok, _ = validate_motif(motif)
                if ok and motif not in self.motifs:
                    self.motifs.append(motif)
                    added += 1
            self.motifs = unique_preserve_order(self.motifs)
            if added:
                invalidate_analysis()
            self.log(f"Import CSV/TXT (dodano {added}, jest {len(self.motifs)})")
            refresh_lists()
            refresh_accordion_checks()
            self.update_status()

        def import_motifs() -> None:
            path = filedialog.askopenfilename(parent=win, title="Import motywów", filetypes=[("CSV", "*.csv"), ("Tekst", "*.txt"), ("Wszystkie", "*.*")])
            if not path:
                return
            file_path = Path(path)
            if file_path.suffix.lower() == ".csv":
                import csv

                with file_path.open(encoding="utf-8", errors="ignore") as handle:
                    reader = csv.DictReader(handle)
                    headers = reader.fieldnames
                    if not headers:
                        messagebox.showerror("Błąd", "Nie znaleziono nagłówków w pliku CSV.", parent=win)
                        return

                col_win = tk.Toplevel(win)
                col_win.title("Wybierz kolumnę z motywami")
                col_win.transient(win)
                col_win.grab_set()
                col_win.resizable(False, False)
                ttk.Label(col_win, text="Wybierz kolumnę:").pack(padx=10, pady=(10, 4))
                selected_col = tk.StringVar(value=headers[0])
                ttk.Combobox(col_win, values=headers, textvariable=selected_col, state="readonly").pack(padx=10, pady=4)

                def confirm_column() -> None:
                    column = selected_col.get()
                    col_win.destroy()
                    with file_path.open(encoding="utf-8", errors="ignore") as handle:
                        reader = csv.DictReader(handle)
                        add_imported_tokens([row[column] for row in reader if row.get(column)])

                ttk.Button(col_win, text="Importuj", command=confirm_column).pack(pady=(4, 10))
                return

            text = file_path.read_text(encoding="utf-8", errors="ignore")
            tokens: list[str] = []
            for line in text.splitlines():
                line = line.strip()
                if not line:
                    continue
                for part in line.replace(";", ",").replace("\t", ",").split(","):
                    part = part.strip()
                    if part:
                        tokens.append(part)
            add_imported_tokens(tokens)

        def export_motifs() -> None:
            path = filedialog.asksaveasfilename(parent=win, title="Eksport motywów", defaultextension=".txt", filetypes=[("Tekst", "*.txt"), ("CSV", "*.csv")])
            if not path:
                return
            file_path = Path(path)
            content = "motif\n" + "\n".join(self.motifs) + "\n" if file_path.suffix.lower() == ".csv" else "\n".join(self.motifs) + "\n"
            file_path.write_text(content, encoding="utf-8")
            self.log(f"Eksport: {file_path.name} ({len(self.motifs)} {plural_motif(len(self.motifs))})")

        motif_entry.bind("<Return>", lambda _event: add_motif(motif_entry.get()))
        ttk.Button(add_row, text="Dodaj", command=lambda: add_motif(motif_entry.get())).grid(row=0, column=1, padx=(8, 0))

        btn_remove = ttk.Button(actions, text="Usuń zaznaczone", command=remove_selected)
        btn_remove.grid(row=0, column=0, sticky="ew", padx=(0, 6))
        btn_clear = ttk.Button(actions, text="Wyczyść motywy", command=clear_all)
        btn_clear.grid(row=0, column=1, sticky="ew", padx=(0, 6))
        ttk.Button(actions, text="Import", command=import_motifs).grid(row=0, column=2, sticky="ew", padx=(0, 6))
        btn_export = ttk.Button(actions, text="Eksport", command=export_motifs)
        btn_export.grid(row=0, column=3, sticky="ew")
        motifs_lb.bind("<<ListboxSelect>>", lambda _event: btn_remove.configure(state=("normal" if motifs_lb.curselection() else "disabled")))

        ttk.Separator(root).grid(row=7, column=0, sticky="ew", pady=(8, 8))
        ttk.Label(root, text="Motywy referencyjne (opcjonalnie):", font=("Segoe UI", 9, "bold")).grid(row=8, column=0, sticky="w")

        accordion_container = ttk.Frame(root)
        accordion_container.grid(row=9, column=0, sticky="ew", pady=(6, 0))
        accordion_container.columnconfigure(0, weight=1)

        eukaryotic_motifs = [
            ("ATG", "start kodon"), ("TAA", "stop kodon"), ("TAG", "stop kodon"), ("TGA", "stop kodon"),
            ("TATAAA", "TATA-box"), ("CCAAT", "CAAT-box"), ("GGGCGG", "GC-box"), ("CANNTG", "E-box"),
            ("TGACGTCA", "CRE"), ("GGGAAA", "NF-κB-like"), ("ATGCAAAT", "Octamer"), ("AATAAA", "Poly-A"),
            ("GTAG", "splice site"),
        ]
        prokaryotic_motifs = [
            ("ATG", "start kodon"), ("TAA", "stop kodon"), ("TAG", "stop kodon"), ("TGA", "stop kodon"),
            ("TATAAT", "Pribnow box"), ("TTGACA", "−35 promotor"), ("AGGAGG", "Shine-Dalgarno"),
            ("TGTGAG", "lac operator"), ("AATGAG", "trp operator"),
        ]

        accordion_frames: dict[str, ttk.Frame] = {}
        icon_vars = {"Eukariota": tk.StringVar(value="▸  Eukariota"), "Prokariota": tk.StringVar(value="▸  Prokariota")}
        section_vars: dict[str, list[tk.BooleanVar]] = {"Eukariota": [], "Prokariota": []}

        def toggle_section(name: str) -> None:
            frame = accordion_frames[name]
            if frame.winfo_ismapped():
                frame.grid_remove()
                icon_vars[name].set(f"▸  {name}")
                autosize()
                return
            for other_name, other_frame in accordion_frames.items():
                other_frame.grid_remove()
                icon_vars[other_name].set(f"▸  {other_name}")
            frame.grid()
            icon_vars[name].set(f"▾  {name}")
            autosize()

        def on_toggle_reference_motif(seq: str, checked: bool) -> None:
            if checked and seq not in self.motifs:
                self.motifs.append(seq)
                self.motifs = unique_preserve_order(self.motifs)
                self.log(f"Dodano motyw: {seq}")
                invalidate_analysis()
            elif not checked and seq in self.motifs:
                self.motifs = [motif for motif in self.motifs if motif != seq]
                self.log(f"Usunięto motyw: {seq}")
                invalidate_analysis()
            refresh_lists()
            refresh_accordion_checks()
            self.update_status()

        def build_section(row: int, name: str, motifs: list[tuple[str, str]]) -> None:
            ttk.Button(accordion_container, textvariable=icon_vars[name], command=lambda: toggle_section(name)).grid(row=row, column=0, sticky="ew", pady=(0, 4))
            frame = ttk.Frame(accordion_container)
            frame.grid(row=row + 1, column=0, sticky="ew")
            frame.grid_remove()
            accordion_frames[name] = frame
            for seq, desc in motifs:
                var = tk.BooleanVar(master=win, value=(seq in self.motifs))
                section_vars[name].append(var)
                ttk.Checkbutton(frame, text=f"{seq} — {desc}", variable=var, command=lambda s=seq, v=var: on_toggle_reference_motif(s, v.get())).pack(anchor="w")

        def refresh_accordion_checks() -> None:
            for name, motifs in (("Eukariota", eukaryotic_motifs), ("Prokariota", prokaryotic_motifs)):
                for (seq, _), var in zip(motifs, section_vars[name]):
                    var.set(seq in self.motifs)

        build_section(0, "Eukariota", eukaryotic_motifs)
        build_section(2, "Prokariota", prokaryotic_motifs)

        self.motifs = unique_preserve_order([normalize_motif(motif) for motif in self.motifs if normalize_motif(motif)])
        refresh_lists()
        refresh_accordion_checks()
        autosize()

    def run_analysis(self) -> None:
        if not self.sequences or not self.motifs:
            messagebox.showwarning("Błąd", "Wczytaj FASTA i dodaj motywy")
            return
        try:
            self.analysis_result = compute_analysis(
                sequences=self.sequences,
                motifs=self.motifs,
                iupac_count_fn=count_matches,
                fasta_path=self.current_fasta,
            )
            self.render_results_table()
            self.draw_visualization()
            self.update_status()
            self.log("Analiza zakończona")
            self.analysis_done = True
        except Exception as exc:
            messagebox.showerror("Błąd analizy", str(exc))

    def sort_column(self, col: str) -> None:
        reverse = self.sort_state.get(col, False)
        rows = [(self.results_table.item(child)["values"], child) for child in self.results_table.get_children()]
        col_index = self.results_table["columns"].index(col)

        if col == "Sekwencja":
            rows.sort(key=lambda item: self._extract_sequence_number(item[0][col_index]), reverse=reverse)
        else:
            try:
                rows.sort(key=lambda item: float(item[0][col_index]), reverse=reverse)
            except ValueError:
                rows.sort(key=lambda item: item[0][col_index], reverse=reverse)

        for index, (_, child) in enumerate(rows):
            self.results_table.move(child, "", index)

        for column in self.results_table["columns"]:
            self.results_table.heading(column, text=column, command=lambda c=column: self.sort_column(c))

        arrow = " ▲" if not reverse else " ▼"
        self.results_table.heading(col, text=col + arrow, command=lambda: self.sort_column(col))
        self.sort_state[col] = not reverse

    @staticmethod
    def _extract_sequence_number(text: str) -> int | float:
        match = re.search(r"\d+", str(text))
        return int(match.group()) if match else float("inf")

    def draw_visualization(self) -> None:
        for widget in self.viz_top.winfo_children():
            widget.destroy()
        if not self.analysis_result:
            return

        data = self.analysis_result.matrix(self.normalization_mode.get())
        self.heatmap_data = data

        fig, ax = plt.subplots(figsize=(8, 5))
        self.figures["heatmap"] = fig
        ax.text(0.5, 1.15, "Kliknij nazwę sekwencji (oś Y) lub motywu (oś X), aby wyświetlić wykres słupkowy.", transform=ax.transAxes, ha="center", va="bottom", fontsize=9, color="dimgray")

        from matplotlib import colors
        norm = colors.Normalize(vmin=0, vmax=np.max(data) if np.max(data) > 0 else 1)
        image = ax.imshow(data, aspect="auto", cmap="viridis_r", norm=norm)

        self.hover_annotation = ax.annotate("", xy=(0, 0), xytext=(10, 10), textcoords="offset points", bbox=dict(boxstyle="round", fc="white", ec="gray"), fontsize=9)
        self.hover_annotation.set_visible(False)

        for row_idx in range(data.shape[0]):
            for col_idx in range(data.shape[1]):
                value = data[row_idx, col_idx]
                text = f"{value:.1f}" if self.normalization_mode.get() == "norm" else f"{int(value)}"
                ax.text(col_idx, row_idx, text, ha="center", va="center", color="black", fontsize=8)

        motifs = self.analysis_result.motifs
        seq_ids = self.analysis_result.seq_ids
        ax.set_xticks(np.arange(len(motifs)))
        ax.set_xticklabels(motifs, rotation=45, ha="right")
        ax.set_yticks(np.arange(len(seq_ids)))
        ax.set_yticklabels(seq_ids)

        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_picker(5)

        colorbar = fig.colorbar(image)
        colorbar.set_label("Liczba motywów" if self.normalization_mode.get() == "raw" else "Na 1000 nt")
        ax.set_title("Heatmapa wystąpień motywów", fontsize=14, fontweight="bold", pad=15)
        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=self.viz_top)
        canvas.draw()
        self.current_canvas = canvas
        canvas.mpl_connect("pick_event", self.on_pick)
        canvas.mpl_connect("motion_notify_event", self.on_hover)
        canvas.get_tk_widget().pack(fill="both", expand=True)

    def on_pick(self, event) -> None:
        if not hasattr(event.artist, "get_text"):
            return
        label = event.artist
        text = label.get_text()

        self._reset_axis_highlights()

        if text in self.sequences:
            self.selected_sequence = text
            self.selected_motif = None
            self.draw_barplot_sequence()
            label.set_color("red")
            label.set_fontweight("bold")
            self.active_label_y = label
        elif text in self.motifs:
            self.selected_motif = text
            self.selected_sequence = None
            self.draw_barplot_motif()
            label.set_color("red")
            label.set_fontweight("bold")
            self.active_label_x = label

        event.canvas.draw_idle()

    def _reset_axis_highlights(self) -> None:
        if self.active_label_x:
            self.active_label_x.set_color("black")
            self.active_label_x.set_fontweight("normal")
            self.active_label_x = None
        if self.active_label_y:
            self.active_label_y.set_color("black")
            self.active_label_y.set_fontweight("normal")
            self.active_label_y = None

    def draw_barplot_sequence(self) -> None:
        if not self.selected_sequence or not self.analysis_result:
            return
        seq_index = self.analysis_result.seq_ids.index(self.selected_sequence)
        values = self.analysis_result.matrix(self.normalization_mode.get())[seq_index, :]
        self._draw_barplot(
            figure_key=f"barplot_seq_{self.selected_sequence}",
            title="Rozkład motywów w sekwencji",
            subtitle=f"Sekwencja: {self.selected_sequence}",
            labels=self.analysis_result.motifs,
            values=values,
        )

    def draw_barplot_motif(self) -> None:
        if not self.selected_motif or not self.analysis_result:
            return
        motif_index = self.analysis_result.motifs.index(self.selected_motif)
        values = self.analysis_result.matrix(self.normalization_mode.get())[:, motif_index]
        self._draw_barplot(
            figure_key=f"barplot_motif_{self.selected_motif}",
            title="Rozkład motywu w sekwencjach",
            subtitle=f"Motyw: {self.selected_motif}",
            labels=self.analysis_result.seq_ids,
            values=values,
        )

    def _draw_barplot(self, figure_key: str, title: str, subtitle: str, labels: list[str], values) -> None:
        self._clear_barplot_area()
        fig, ax = plt.subplots(figsize=(6, 4))
        self.figures[figure_key] = fig
        fig.suptitle(title, fontsize=14, fontweight="bold")

        x_positions = np.arange(len(labels))
        ax.bar(x_positions, values)
        ax.set_xticks(x_positions)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_title(subtitle, fontsize=11, pad=5)
        ax.set_ylabel("Na 1000 nt" if self.normalization_mode.get() == "norm" else "Liczba")
        fig.tight_layout(rect=[0, 0, 1, 0.88])
        fig.subplots_adjust(bottom=0.30)

        self.viz_bottom.configure(height=260)
        canvas = FigureCanvasTkAgg(fig, master=self.viz_bottom)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=10)
        self.close_barplot_btn.lift()
        self.close_barplot_btn.place(relx=1.0, rely=0.0, anchor="ne", x=-6, y=6)
        self.autosize_viz_window()

    def _clear_barplot_area(self) -> None:
        for widget in self.viz_bottom.winfo_children():
            if widget is not self.close_barplot_btn:
                widget.destroy()

    def on_hover(self, event) -> None:
        if event.inaxes is None:
            if self.hover_annotation:
                self.hover_annotation.set_visible(False)
                event.canvas.draw_idle()
            return
        if self.heatmap_data is None or event.xdata is None or event.ydata is None:
            return

        x = int(round(event.xdata))
        y = int(round(event.ydata))
        if 0 <= x < len(self.motifs) and 0 <= y < len(self.sequences):
            sequence_id = list(self.sequences.keys())[y]
            motif = self.motifs[x]
            self.hover_annotation.xy = (x, y)
            self.hover_annotation.set_text(f"{sequence_id}\n{motif}")
            self.hover_annotation.set_visible(True)
            event.canvas.draw_idle()
            return

        self.hover_annotation.set_visible(False)
        event.canvas.draw_idle()

    def refresh_visualization(self) -> None:
        if not self.analysis_result:
            return
        self.render_results_table()
        self.draw_visualization()
        if self.selected_sequence:
            self.draw_barplot_sequence()
        elif self.selected_motif:
            self.draw_barplot_motif()
        self.update_status()

    def hide_barplot(self) -> None:
        self._clear_barplot_area()
        self.viz_bottom.configure(height=0)
        self.close_barplot_btn.place_forget()
        self._reset_axis_highlights()
        self.selected_sequence = None
        self.selected_motif = None
        if self.current_canvas:
            self.current_canvas.draw_idle()
        self.autosize_viz_window()

    def download_ncbi(self) -> None:
        win = tk.Toplevel(self)
        win.title("Pobierz FASTA z NCBI")
        win.resizable(False, False)
        win.transient(self)
        win.grab_set()

        frm = ttk.Frame(win, padding=12)
        frm.pack(fill="both", expand=True)
        frm.columnconfigure(0, weight=1)

        ttk.Label(frm, text="Wpisz accession/UID (jeden lub kilka, oddzielone przecinkami lub spacją):", font=("Segoe UI", 10, "bold")).grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 8))
        entry = ttk.Entry(frm, width=56)
        entry.grid(row=1, column=0, sticky="ew")
        entry.focus_set()

        ttk.Button(
            frm,
            text="?",
            width=3,
            command=lambda: messagebox.showinfo(
                "Podpowiedź",
                "Podaj accession/UID, np.:\n• NC_001416.1\n• NM_000546.6\n• MN908947.3\n• AY274119.3",
                parent=win,
            ),
        ).grid(row=1, column=1, padx=(8, 0))

        def do_download() -> None:
            if not self._confirm_replace_current_analysis(parent=win, source_label="nową sekwencję"):
                return
            if self.analysis_done:
                self.clear_results()

            query = entry.get().strip()
            if not query:
                messagebox.showwarning("Brak danych", "Wpisz accession/UID.", parent=win)
                return

            ids = [item for item in re.split(r"[,\s]+", query) if item]
            if not ids or not all(is_accession_like(item) for item in ids):
                messagebox.showwarning("Nieprawidłowe dane", "Podaj accession number lub UID, bez wyszukiwania tekstowego.", parent=win)
                return

            try:
                self.status_var.set("Pobieram FASTA z NCBI...")
                win.update_idletasks()
                fasta_text = fetch_fasta_by_ids(ids)
            except Exception as exc:
                messagebox.showerror("Błąd pobierania", str(exc), parent=win)
                self.update_status()
                return

            out_path = filedialog.asksaveasfilename(parent=win, title="Zapisz pobrane FASTA", defaultextension=".fasta", filetypes=[("FASTA (DNA/RNA)", "*.fasta *.fa *.fna"), ("Wszystkie", "*.*")])
            if not out_path:
                self.update_status()
                return

            Path(out_path).write_text(fasta_text, encoding="utf-8")
            self.log(f"Pobrano FASTA z NCBI: {out_path}")

            try:
                self.load_sequences_from_file(out_path)
                self.log(f"Wczytano {len(self.sequences)} sekwencji")
                self.tabs.select(self.tab_preview)
                self.update_status()
                win.destroy()
            except Exception as exc:
                messagebox.showerror("Błąd FASTA", str(exc), parent=win)
                self.update_status()

        btns = ttk.Frame(win, padding=(12, 0, 12, 12))
        btns.pack(side="bottom", fill="x")
        ttk.Button(btns, text="Anuluj", command=win.destroy).pack(side="right")
        ttk.Button(btns, text="Pobierz", command=do_download).pack(side="right", padx=(0, 8))
        entry.bind("<Return>", lambda _event: do_download())

    def render_results_table(self) -> None:
        self.results_table.delete(*self.results_table.get_children())
        if not self.analysis_result:
            return

        mode = self.normalization_mode.get()
        mat = self.analysis_result.matrix(mode)
        columns = ["Sekwencja"] + self.analysis_result.motifs + ["SUMA"]
        self.results_table["columns"] = columns

        for column in columns:
            self.results_table.heading(column, text=f"{column} ⇅", command=lambda c=column: self.sort_column(c))
            self.results_table.column(column, width=90, anchor="center")

        for row_idx, sequence_id in enumerate(self.analysis_result.seq_ids):
            row_vals = mat[row_idx, :]
            cells = [str(int(value)) for value in row_vals] if mode == "raw" else [f"{value:.1f}" for value in row_vals]
            total = str(int(row_vals.sum())) if mode == "raw" else f"{row_vals.sum():.1f}"
            self.results_table.insert("", "end", values=[sequence_id] + cells + [total])

    def new_analysis(self) -> None:
        if not self.analysis_done and not self.sequences and not self.motifs:
            return
        ok = messagebox.askyesno(
            "Nowa analiza",
            "Czy na pewno chcesz rozpocząć nową analizę?\nPoprzednie wyniki, wizualizacje, FASTA i motywy zostaną usunięte.",
        )
        if ok:
            self.clear_results()
            self.log("Rozpoczęto nową analizę.")

    def clear_results(self) -> None:
        self.analysis_result = None
        self.figures.clear()
        self.sequences.clear()
        self.current_fasta = None
        self.motifs.clear()
        self.analysis_done = False
        self.preview_box.delete("1.0", "end")
        self.log_box.delete("1.0", "end")
        self._clear_results_views_only()
        self.status_var.set("Gotowe")
        self.tabs.select(self.tab_preview)
        self.update_status()

    def _clear_results_views_only(self) -> None:
        self.results_table.delete(*self.results_table.get_children())
        for widget in self.viz_top.winfo_children():
            widget.destroy()
        self._clear_barplot_area()
        self.viz_bottom.configure(height=0)
        self.close_barplot_btn.place_forget()
        self.current_canvas = None
        self.hover_annotation = None
        self.heatmap_data = None
        self.selected_sequence = None
        self.selected_motif = None
        self._reset_axis_highlights()

    def autosize_viz_window(self, pad: int = 16) -> None:
        win = self.viz_frame.winfo_toplevel()
        win.update_idletasks()
        win.geometry(f"{win.winfo_reqwidth() + pad}x{win.winfo_reqheight() + pad}")

    def open_help_window(self) -> None:
        if self.help_win is not None and self.help_win.winfo_exists():
            self.help_win.lift()
            self.help_win.focus_force()
            return

        self.help_win = tk.Toplevel(self)
        self.help_win.title("Pomoc")
        self.help_win.geometry("780x520")
        self.help_win.minsize(520, 360)

        top = ttk.Frame(self.help_win)
        top.pack(fill="x", padx=10, pady=(10, 0))
        txt = ScrolledText(self.help_win, wrap="word")
        txt.pack(fill="both", expand=True, padx=10, pady=10)

        help_text = """O aplikacji
Program służy do wczytywania plików FASTA i wyszukiwania w sekwencjach zdefiniowanych motywów nukleotydowych z obsługą kodów IUPAC (np. N, R, Y). Wyniki prezentowane są w formie tabelarycznej oraz na wykresach.

Szybki start
1) Plik → Wczytaj plik (FASTA)
2) Motywy → Dodaj/usuń motywy
3) Uruchom analizę
4) Zakładki: Podgląd / Wyniki / Wizualizacja / Eksport

Format FASTA
- Nagłówek zaczyna się od '>'
- Sekwencja może być w wielu liniach
- Program usuwa spacje, taby i cyfry, zamienia '?' na 'N' i waliduje znaki

Dozwolone znaki
A C G T U R Y S W K M B D H V N

Uwaga: wyszukiwanie uwzględnia nakładające się dopasowania.

Tryb wyników
- Surowe: liczba wystąpień
- Na 1000 nt: normalizacja względem długości sekwencji
"""

        txt.insert("1.0", help_text)
        txt.configure(state="disabled")

        def copy_help() -> None:
            self.clipboard_clear()
            self.clipboard_append(help_text)
            self.update_idletasks()
            messagebox.showinfo("Pomoc", "Skopiowano treść pomocy do schowka.")

        def save_help() -> None:
            path = filedialog.asksaveasfilename(title="Zapisz pomoc", defaultextension=".txt", filetypes=[("Plik tekstowy", "*.txt"), ("Wszystkie pliki", "*.*")])
            if not path:
                return
            try:
                with open(path, "w", encoding="utf-8") as handle:
                    handle.write(help_text)
                messagebox.showinfo("Pomoc", "Zapisano plik pomocy.")
            except Exception as exc:
                messagebox.showerror("Błąd zapisu", str(exc))

        ttk.Button(top, text="Kopiuj do schowka", command=copy_help).pack(side="left")
        ttk.Button(top, text="Zapisz do pliku…", command=save_help).pack(side="left", padx=(8, 0))
        ttk.Button(top, text="Zamknij", command=self.help_win.destroy).pack(side="right")
