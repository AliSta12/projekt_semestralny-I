"""
Projekt 1 – Analiza motywów sekwencyjnych DNA

Opis:
Program umożliwia wczytanie sekwencji DNA w formacie FASTA,
dodawanie motywów nukleotydowych, analizę ich liczby oraz pozycji,
a także wizualizację wyników i eksport do plików CSV.

Architektura programu:
1. LOGIKA
   - funkcje niezależne od GUI (FASTA, motywy)
2. GUI
   - interfejs użytkownika (Tkinter)
   - obsługa zdarzeń
   - wizualizacja i eksport
"""

# ==================================================
# IMPORTY
# ==================================================

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog
import csv
import os
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np


# ==================================================
# LOGIKA – PRZETWARZANIE SEKWENCJI DNA
# ==================================================

def load_fasta(path):
    """
    Wczytuje plik FASTA i zwraca słownik:
    {ID_sekwencji : sekwencja_DNA}
    """

    sequences = {}
    current_id = None

    with open(path, "r") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()

            # pomijamy puste linie
            if not line:
                continue

            # nagłówek FASTA
            if line.startswith(">"):
                current_id = line[1:].strip()

                if not current_id:
                    raise ValueError(f"Pusty nagłówek w linii {line_num}")

                if current_id in sequences:
                    raise ValueError(f"Duplikat ID: {current_id}")

                sequences[current_id] = ""

            # linia sekwencji
            else:
                if current_id is None:
                    raise ValueError("Plik nie zaczyna się od nagłówka FASTA")

                if any(c not in "ACGTacgt" for c in line):
                    raise ValueError(f"Niepoprawne znaki w linii {line_num}")

                sequences[current_id] += line.upper()

    if not sequences:
        raise ValueError("Nie znaleziono sekwencji")

    return sequences


def find_motif_positions(sequence, motif):
    """
    Zwraca listę pozycji (indeksów),
    w których dany motyw występuje w sekwencji
    """

    positions = []
    start = 0

    while True:
        pos = sequence.find(motif, start)
        if pos == -1:
            break

        positions.append(pos)
        start = pos + 1

    return positions


def count_motif(sequence, motif):
    """
    Liczy liczbę wystąpień motywu w sekwencji
    (z uwzględnieniem nachodzących na siebie motywów)
    """

    count = 0
    start = 0

    while True:
        pos = sequence.find(motif, start)
        if pos == -1:
            break

        count += 1
        start = pos + 1

    return count
# ==================================================
# IUPAC – MAPA KODÓW I FUNKCJE DOPASOWANIA
# ==================================================

IUPAC_MAP = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"},
    "S": {"G", "C"}, "W": {"A", "T"},
    "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

def matches_iupac(seq_char: str, motif_char: str) -> bool:
    """
    Sprawdza, czy nukleotyd seq_char pasuje
    do symbolu IUPAC motif_char.
    """
    seq_char = seq_char.upper()
    motif_char = motif_char.upper()
    return seq_char in IUPAC_MAP.get(motif_char, set())


def iupac_find_positions(sequence: str, motif: str) -> list[int]:
    """
    Znajduje pozycje motywu z IUPAC w sekwencji.
    Zwraca listę pozycji (0-based).
    """
    positions = []
    L = len(sequence)
    mL = len(motif)

    for i in range(L - mL + 1):
        match = True
        for j in range(mL):
            if not matches_iupac(sequence[i + j], motif[j]):
                match = False
                break
        if match:
            positions.append(i)

    return positions


def iupac_count(sequence: str, motif: str) -> int:
    """
    Liczy liczbę dopasowań motywu IUPAC w sekwencji,
    z allow overlapping.
    """
    return len(iupac_find_positions(sequence, motif))



# ==================================================
# GUI – APLIKACJA TKINTER
# ==================================================

class DNAApp(tk.Tk):

    # ==================================================
    # INICJALIZACJA APLIKACJI
    # ==================================================

    def __init__(self):
        super().__init__()

        self.title("Projekt 1: Analiza motywów sekwencyjnych w DNA")
        self.geometry("1000x600")

        # dane aplikacji
        self.sequences = {}
        self.current_fasta = None
        self.motifs = []
        self.DEV_MODE = True  # <-- zmień na False gdy niepotrzebne

        # budowa interfejsu
        self.sort_state = {}  # zapamiętuje kierunek sortowania kolumn
        self.normalization_mode = tk.StringVar(value="raw")

        self.create_menu()
        self.create_layout()
        if self.DEV_MODE:
            self.load_test_data()

    # ==================================================
    # MENU GÓRNE
    # ==================================================

    def create_menu(self):
        """Tworzy menu główne aplikacji"""

        menubar = tk.Menu(self)

        plik = tk.Menu(menubar, tearoff=0)
        plik.add_command(label="Wczytaj plik", command=self.gui_load_file)
        plik.add_separator()
        plik.add_command(label="Wyjście", command=self.quit)

        motywy = tk.Menu(menubar, tearoff=0)
        motywy.add_command(label="Dodaj / usuń motywy", command=self.open_motif_manager)

        ncbi = tk.Menu(menubar, tearoff=0)
        ncbi.add_command(label="Pobierz z NCBI", command=self.download_ncbi)

        eksport = tk.Menu(menubar, tearoff=0)
        eksport.add_command(label="Eksportuj CSV/PDF", command=self.export_data)

        pomoc = tk.Menu(menubar, tearoff=0)
        pomoc.add_command(label="O programie", command=self.about)

        menubar.add_cascade(label="Plik", menu=plik)
        menubar.add_cascade(label="Motywy", menu=motywy)
        menubar.add_cascade(label="NCBI", menu=ncbi)
        menubar.add_cascade(label="Eksport", menu=eksport)
        menubar.add_cascade(label="Pomoc", menu=pomoc)

        self.config(menu=menubar)

    # ==================================================
    # UKŁAD OKNA I WIDŻETY
    # ==================================================

    def create_layout(self):
        """
        Buduje cały interfejs graficzny:
        panele, zakładki, tabele, wizualizację i logi
        """

        # --- główny kontener ---
        main_frame = tk.Frame(self)
        main_frame.pack(fill="both", expand=True)

        # --- panel boczny ---
        side = tk.Frame(main_frame, width=200)
        side.pack(side="left", fill="y")

        tk.Button(side, text="Wczytaj plik", command=self.gui_load_file).pack(fill="x", pady=2)
        tk.Button(side, text="Pobierz z NCBI", command=self.download_ncbi).pack(fill="x", pady=2)
        tk.Button(
            side,
            text="Dodaj / usuń motywy",
            command=self.open_motif_manager
        ).pack(fill="x", pady=2)
        tk.Button(side, text="Uruchom analizę", command=self.run_analysis).pack(fill="x", pady=2)
        tk.Button(side, text="Eksportuj CSV/PDF", command=self.export_data).pack(fill="x", pady=2)

        # --- część prawa ---
        right = tk.Frame(main_frame)
        right.pack(side="right", fill="both", expand=True)

        self.tabs = ttk.Notebook(right)
        self.tabs.pack(fill="both", expand=True)

        # zakładki
        self.tab_preview = ttk.Frame(self.tabs)
        self.tab_motifs = ttk.Frame(self.tabs)
        self.tab_results = ttk.Frame(self.tabs)
        self.tab_viz = ttk.Frame(self.tabs)
        self.tab_export = ttk.Frame(self.tabs)

        self.tabs.add(self.tab_preview, text="Podgląd sekwencji")
        self.tabs.add(self.tab_motifs, text="Szukane motywy")
        self.tabs.add(self.tab_results, text="Wyniki analizy")
        self.tabs.add(self.tab_viz, text="Wizualizacja")
        self.tabs.add(self.tab_export, text="Eksport")

        # --- podgląd FASTA ---
        self.preview_box = tk.Text(self.tab_preview)
        self.preview_box.pack(fill="both", expand=True)

        # --- lista motywów ---
        self.motif_listbox = tk.Listbox(self.tab_motifs)
        self.motif_listbox.pack(fill="both", expand=True)

        # --- tabela wyników (macierzowa) ---
        results_frame = tk.Frame(self.tab_results)
        results_frame.pack(fill="both", expand=True)

        # ===== TRYB WYŚWIETLANIA =====
        mode_frame = tk.Frame(results_frame)
        mode_frame.pack(fill="x", pady=(5, 0))

        tk.Label(
            mode_frame,
            text="Tryb wyświetlania:"
        ).pack(side="left", padx=(5, 5))

        tk.Radiobutton(
            mode_frame,
            text="Surowe liczby",
            variable=self.normalization_mode,
            value="raw",
            command=self.run_analysis
        ).pack(side="left")

        tk.Radiobutton(
            mode_frame,
            text="Na 1000 nt",
            variable=self.normalization_mode,
            value="norm",
            command=self.run_analysis
        ).pack(side="left")

        # ===== TABELA =====
        self.results_table = ttk.Treeview(
            results_frame,
            show="headings"
        )
        self.results_table.pack(fill="both", expand=True)

        # zmiana kursora na rączkę przy najechaniu na nagłówki
        def on_heading_enter(event):
            region = self.results_table.identify_region(event.x, event.y)
            if region == "heading":
                self.results_table.config(cursor="hand2")
            else:
                self.results_table.config(cursor="")

        self.results_table.bind("<Motion>", on_heading_enter)

        # --- wizualizacja ---
        self.viz_frame = tk.Frame(self.tab_viz)
        self.viz_frame.pack(fill="both", expand=True)

        self.viz_canvas = tk.Canvas(self.viz_frame, bg="white")
        self.viz_canvas.pack(side="left", fill="both", expand=True)

        self.viz_scroll = ttk.Scrollbar(
            self.viz_frame,
            orient="vertical",
            command=self.viz_canvas.yview
        )
        self.viz_scroll.pack(side="right", fill="y")

        self.viz_canvas.configure(yscrollcommand=self.viz_scroll.set)

        self.viz_inner = tk.Frame(self.viz_canvas)
        self.viz_canvas.create_window((0, 0), window=self.viz_inner, anchor="nw")

        self.viz_inner.bind(
            "<Configure>",
            lambda e: self.viz_canvas.configure(
                scrollregion=self.viz_canvas.bbox("all")
            )
        )

        # --- eksport ---
        tk.Label(self.tab_export, text="Pliki wynikowe CSV:").pack(anchor="w")
        self.export_listbox = tk.Listbox(self.tab_export)
        self.export_listbox.pack(fill="both", expand=True)

        tk.Button(
            self.tab_export,
            text="Odśwież listę",
            command=self.refresh_csv_files
        ).pack(pady=5)

        # --- logi ---
        log_frame = tk.Frame(self)
        log_frame.pack(fill="x")

        tk.Label(log_frame, text="Logi / komunikaty:").pack(anchor="w")
        self.log_box = tk.Text(log_frame, height=6)
        self.log_box.pack(fill="x")

    # ==================================================
    # FUNKCJE POMOCNICZE
    # ==================================================

    def log(self, msg):
        """Dodaje komunikat do pola logów"""
        self.log_box.insert("end", msg + "\n")
        self.log_box.see("end")

    def refresh_csv_files(self):
        """Odświeża listę plików CSV w zakładce Eksport"""

        self.export_listbox.delete(0, "end")
        project_dir = os.path.dirname(os.path.abspath(__file__))

        for f in os.listdir(project_dir):
            if f.lower().endswith(".csv"):
                self.export_listbox.insert("end", f)

        if self.export_listbox.size() == 0:
            self.export_listbox.insert("end", "(brak plików CSV)")

    def load_test_data(self):
        """Ładuje przykładowe dane testowe (tryb developerski)"""

        # --- przykładowe sekwencje ---
        self.sequences = {
            "seq1": "ATGCGTATGCGTATGTATAAAGGGCGGATGCGT",
            "seq2": "TTGACATATAATAGGAGGATGCGTTAA",
            "seq10": "GGGAAATGCAAATAATAAATGCGTATG",
            "seq20": "TATAAAGGGCGGCANNTGTGACGTCA"
        }

        # --- przykładowe motywy ---
        self.motifs = ["ATG", "TATAAA", "GGGCGG"]

        # wyczyść GUI
        self.preview_box.delete("1.0", "end")
        self.motif_listbox.delete(0, "end")

        # wstaw podgląd sekwencji
        for k, v in self.sequences.items():
            self.preview_box.insert("end", f">{k}\n{v}\n\n")

        # wstaw motywy
        for m in self.motifs:
            self.motif_listbox.insert("end", m)

        self.log("Załadowano dane testowe (DEV_MODE)")

    # ==================================================
    # AKCJE UŻYTKOWNIKA
    # ==================================================

    def gui_load_file(self):
        """Wczytuje plik FASTA wybrany przez użytkownika"""

        path = filedialog.askopenfilename(
            title="Wybierz FASTA",
            filetypes=[("FASTA", "*.fasta *.fa *.fna"), ("All", "*.*")]
        )

        if not path:
            return

        try:
            self.sequences.clear()
            self.preview_box.delete("1.0", "end")

            self.sequences = load_fasta(path)
            self.current_fasta = path

            for i, (k, v) in enumerate(self.sequences.items()):
                if i == 20:
                    self.preview_box.insert("end", "\n... (reszta ukryta)\n")
                    break

                self.preview_box.insert("end", f">{k}\n{v[:200]}...\n\n")

            self.log(f"Wczytano plik: {path}")

        except Exception as e:
            messagebox.showerror("Błąd FASTA", str(e))

    def open_motif_manager(self):
        """
        Okno do zarządzania motywami:
        - wpisywanie własnego
        - lista szukanych
        - akordeon z predefiniowanymi
        """

        win = tk.Toplevel(self)
        win.title("Dodaj/usuń motywy")
        # Ustal stałą szerokość, dynamiczną wysokość
        initial_width = 380

        win.geometry(f"{initial_width}x500")
        win.resizable(False, True)
        win.transient(self)
        win.grab_set()

        # ===== LEGENDA =====
        tk.Label(win, text="Dozwolone symbole:").pack(pady=(6, 0))
        tk.Label(win, text="A  C  G  T").pack()

        iupac_label = tk.Label(win, text="R Y S W K M B D H V N (IUPAC)", fg="blue")
        iupac_label.pack(pady=(0, 8))

        tooltip_text = (
            "R = A lub G\n"
            "Y = C lub T\n"
            "S = G lub C\n"
            "W = A lub T\n"
            "K = G lub T\n"
            "M = A lub C\n"
            "B = C G T\n"
            "D = A G T\n"
            "H = A C T\n"
            "V = A C G\n"
            "N = dowolna baza"
        )

        tooltip = None

        def show_tooltip(e):
            nonlocal tooltip
            if tooltip: return
            x, y = e.x_root + 10, e.y_root + 10
            tooltip = tk.Toplevel(win)
            tooltip.wm_overrideredirect(True)
            tooltip.geometry(f"+{x}+{y}")
            tk.Label(
                tooltip, text=tooltip_text,
                justify="left", bg="#ffffe0",
                relief="solid", borderwidth=1,
                padx=4, pady=2
            ).pack()

        def hide_tooltip(e):
            nonlocal tooltip
            if tooltip:
                tooltip.destroy()
                tooltip = None

        iupac_label.bind("<Enter>", show_tooltip)
        iupac_label.bind("<Leave>", hide_tooltip)

        # ===== Własny motyw =====
        tk.Label(win, text="Wpisz swój motyw:", anchor="w").pack(fill="x", padx=10)
        motif_entry = tk.Entry(win)
        motif_entry.pack(fill="x", padx=10, pady=(2, 6))

        def add_manual(event):
            m = motif_entry.get().strip().upper()
            if not m:
                return
            if any(c not in IUPAC_MAP for c in m):
                messagebox.showerror("Błąd", "Motyw może zawierać tylko symbole ACGT lub IUPAC")
                return
            if m in self.motifs:
                messagebox.showinfo("Info", "Taki motyw już istnieje")
                return
            self.motifs.append(m)
            self.motif_listbox.insert("end", m)
            searched_listbox.insert("end", m)
            motif_entry.delete(0, "end")
            self.log(f"Dodano motyw: {m}")

        motif_entry.bind("<Return>", add_manual)

        # ===== Szukane motywy =====
        tk.Label(win, text="Szukane motywy:", anchor="w").pack(fill="x", padx=10)
        searched_listbox = tk.Listbox(win, selectmode="extended")
        searched_listbox.pack(fill="both", expand=True, padx=10, pady=(0, 8))

        # ===== Akordeon =====
        tk.Label(win, text="Wybierz z listy motywów:", anchor="w").pack(
            fill="x", padx=10, pady=(4, 6)
        )

        # kontener akordeonu
        accordion_container = tk.Frame(win)
        accordion_container.pack(fill="x", padx=10, pady=(0, 6))

        EU_MOTIFS = [
            ("ATG", "start kodon"),
            ("TAA", "stop kodon"),
            ("TAG", "stop kodon"),
            ("TGA", "stop kodon"),
            ("TATAAA", "TATA-box"),
            ("CCAAT", "CAAT-box"),
            ("GGGCGG", "GC-box"),
            ("CANNTG", "E-box"),
            ("TGACGTCA", "CRE"),
            ("GGGAAA", "NF-κB-like"),
            ("ATGCAAAT", "Octamer"),
            ("AATAAA", "Poly-A"),
            ("GTAG", "splice site"),
        ]

        PRO_MOTIFS = [
            ("ATG", "start kodon"),
            ("TAA", "stop kodon"),
            ("TAG", "stop kodon"),
            ("TGA", "stop kodon"),
            ("TATAAT", "Pribnow box"),
            ("TTGACA", "−35 promotor"),
            ("AGGAGG", "Shine-Dalgarno"),
            ("TGTGAG", "lac operator"),
            ("AATGAG", "trp operator"),
        ]

        accordion_frames = {}
        icons = {
            "Eukariota": tk.StringVar(value="▸  Eukariota"),
            "Prokariota": tk.StringVar(value="▸  Prokariota"),
        }
        checkbox_vars = {}

        def toggle_section(name):
            # jeśli ta sekcja już jest widoczna → chowaj
            if accordion_frames[name].winfo_ismapped():
                accordion_frames[name].pack_forget()
                icons[name].set(f"▸  {name}")

                win.update_idletasks()
                h = win.winfo_reqheight()
                win.geometry(f"{initial_width}x{h}")
                return

            # chowamy wszystkie i reset ikon
            for nm, frame in accordion_frames.items():
                frame.pack_forget()
                icons[nm].set(f"▸  {nm}")

            # pokaz sekcję
            accordion_frames[name].pack(fill="x", padx=20, pady=(2, 6))
            icons[name].set(f"▾  {name}")

            win.update_idletasks()
            h = win.winfo_reqheight()
            win.geometry(f"{initial_width}x{h}")

        # — Eukariota — #
        btn_euk = tk.Button(
            accordion_container,
            textvariable=icons["Eukariota"],
            anchor="w",
            fg="#0b5394",
            font=("Arial", 10, "bold"),
            relief="flat",
            command=lambda: toggle_section("Eukariota")
        )
        btn_euk.pack(fill="x")

        frm_euk = tk.Frame(accordion_container)
        accordion_frames["Eukariota"] = frm_euk
        checkbox_vars["Eukariota"] = []

        for seq, desc in EU_MOTIFS:
            var = tk.BooleanVar(master=win, value=(seq in self.motifs))

            def on_check(v=var, s=seq):
                if v.get():
                    if s not in self.motifs:
                        self.motifs.append(s)
                        self.motif_listbox.insert("end", s)
                        searched_listbox.insert("end", s)
                else:
                    if s in self.motifs:
                        self.motifs.remove(s)
                        try:
                            idx_main = self.motifs.index(s)
                            self.motif_listbox.delete(idx_main)
                        except ValueError:
                            pass
                    items = list(searched_listbox.get(0, "end"))
                    if s in items:
                        searched_listbox.delete(items.index(s))

            cb = tk.Checkbutton(
                frm_euk,
                text=f"{seq} — {desc}",
                variable=var,
                command=on_check
            )
            cb.pack(anchor="w")
            checkbox_vars["Eukariota"].append((var, seq))

        # chowamy na start
        frm_euk.pack_forget()

        # — Prokariota — #
        btn_pro = tk.Button(
            accordion_container,
            textvariable=icons["Prokariota"],
            anchor="w",
            fg="#0b5394",
            font=("Arial", 10, "bold"),
            relief="flat",
            command=lambda: toggle_section("Prokariota")
        )
        btn_pro.pack(fill="x", pady=(4, 0))

        frm_pro = tk.Frame(accordion_container)
        accordion_frames["Prokariota"] = frm_pro
        checkbox_vars["Prokariota"] = []

        for seq, desc in PRO_MOTIFS:
            var = tk.BooleanVar(master=win, value=(seq in self.motifs))

            def on_check_p(v=var, s=seq):
                if v.get():
                    if s not in self.motifs:
                        self.motifs.append(s)
                        self.motif_listbox.insert("end", s)
                        searched_listbox.insert("end", s)
                else:
                    if s in self.motifs:
                        self.motifs.remove(s)
                        try:
                            idx_main = self.motifs.index(s)
                            self.motif_listbox.delete(idx_main)
                        except ValueError:
                            pass
                    items = list(searched_listbox.get(0, "end"))
                    if s in items:
                        searched_listbox.delete(items.index(s))

            cb = tk.Checkbutton(
                frm_pro,
                text=f"{seq} — {desc}",
                variable=var,
                command=on_check_p
            )
            cb.pack(anchor="w")
            checkbox_vars["Prokariota"].append((var, seq))

        frm_pro.pack_forget()

    def run_analysis(self):
        """Buduje macierzową tabelę pivot"""

        if not self.sequences or not self.motifs:
            messagebox.showwarning("Błąd", "Wczytaj FASTA i dodaj motywy")
            return

        # wyczyść tabelę
        for row in self.results_table.get_children():
            self.results_table.delete(row)

        # dynamiczne kolumny
        columns = ["Sekwencja"] + self.motifs + ["SUMA"]
        self.results_table["columns"] = columns

        for col in columns:
            self.results_table.heading(
                col,
                text=f"{col} ⇅",
                command=lambda c=col: self.sort_column(c)
            )
            self.results_table.column(col, width=90, anchor="center")

        # wypełnianie tabeli
        for seq_id, seq in self.sequences.items():
            row = [seq_id]
            total = 0
            seq_length = len(seq)

            for motif in self.motifs:
                count = iupac_count(seq, motif)

                if self.normalization_mode.get() == "norm" and seq_length > 0:
                    value = round((count / seq_length) * 1000, 1)
                else:
                    value = count

                row.append(value)
                total += value

            row.append(round(total, 1))
            self.results_table.insert("", "end", values=row)

        self.log("Analiza zakończona")
        self.draw_visualization()

    def sort_column(self, col):
        """Sortowanie kolumny + strzałki w nagłówkach"""

        # ustalamy kierunek
        reverse = self.sort_state.get(col, False)

        data = []
        for child in self.results_table.get_children():
            values = self.results_table.item(child)["values"]
            data.append((values, child))

        col_index = self.results_table["columns"].index(col)

        # specjalne sortowanie dla kolumny Sekwencja (numeryczne)
        if col == "Sekwencja":
            def extract_number(text):
                match = re.search(r"\d+", str(text))
                return int(match.group()) if match else float("inf")

            data.sort(
                key=lambda x: extract_number(x[0][col_index]),
                reverse=reverse
            )

        else:
            try:
                data.sort(key=lambda x: float(x[0][col_index]), reverse=reverse)
            except ValueError:
                data.sort(key=lambda x: x[0][col_index], reverse=reverse)

        # przestaw wiersze
        for index, (values, child) in enumerate(data):
            self.results_table.move(child, "", index)

        # reset wszystkich nagłówków
        for c in self.results_table["columns"]:
            self.results_table.heading(
                c,
                text=c,
                command=lambda col=c: self.sort_column(col)
            )

        # ustaw strzałkę tylko dla aktywnej kolumny
        arrow = " ▲" if not reverse else " ▼"
        self.results_table.heading(
            col,
            text=col + arrow,
            command=lambda: self.sort_column(col)
        )

        # zapamiętaj nowy stan
        self.sort_state[col] = not reverse

    # ==================================================
    # WIZUALIZACJA I PODSUMOWANIA
    # ==================================================

    def draw_visualization(self):
        """Rysuje heatmapę motywów"""

        # wyczyść poprzednią zawartość
        for w in self.viz_inner.winfo_children():
            w.destroy()

        if not self.sequences or not self.motifs:
            return

        # budujemy macierz danych
        data = []

        for seq_id, seq in self.sequences.items():
            row = []
            seq_length = len(seq)

            for motif in self.motifs:
                count = iupac_count(seq, motif)

                if self.normalization_mode.get() == "norm" and seq_length > 0:
                    value = (count / seq_length) * 1000
                else:
                    value = count

                row.append(value)

            data.append(row)

        data = np.array(data)

        # tworzymy wykres
        fig, ax = plt.subplots(figsize=(8, 5))

        from matplotlib import colors

        norm = colors.Normalize(
            vmin=0,
            vmax=np.max(data) if np.max(data) > 0 else 1
        )

        im = ax.imshow(
            data,
            aspect="auto",
            cmap="viridis_r",
            norm=norm
        )

        # wpisanie wartości do komórek
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                value = data[i, j]

                if self.normalization_mode.get() == "norm":
                    text = f"{value:.1f}"
                else:
                    text = f"{int(value)}"

                ax.text(
                    j,
                    i,
                    text,
                    ha="center",
                    va="center",
                    color="black",
                    fontsize=8
                )

        # etykiety osi
        ax.set_xticks(np.arange(len(self.motifs)))
        ax.set_xticklabels(self.motifs, rotation=45, ha="right")

        ax.set_yticks(np.arange(len(self.sequences)))
        ax.set_yticklabels(list(self.sequences.keys()))

        # colorbar
        cbar = fig.colorbar(im)
        cbar.set_label("Liczba motywów" if self.normalization_mode.get() == "raw" else "Na 1000 nt")

        ax.set_title("Heatmapa wystąpień motywów")

        fig.tight_layout()

        # osadzamy w Tkinter
        canvas = FigureCanvasTkAgg(fig, master=self.viz_inner)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

    # ==================================================
    # EKSPORT I INFORMACJE
    # ==================================================

    def export_data(self):
        """Eksportuje wyniki analizy do pliku CSV"""

        if not self.results_table.get_children():
            messagebox.showinfo("Eksport", "Brak wyników do eksportu.")
            return

        path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv")]
        )

        if not path:
            return

        with open(path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Sekwencja", "Motyw", "Liczba"])

            for row in self.results_table.get_children():
                writer.writerow(self.results_table.item(row)["values"])

        self.log(f"Wyeksportowano wyniki do: {path}")

    def download_ncbi(self):
        """Placeholder – pobieranie danych z NCBI"""
        self.log("NCBI (placeholder)")

    def about(self):
        """Informacje o programie"""
        messagebox.showinfo(
            "O programie",
            "Analiza motywów DNA\nProjekt zaliczeniowy"
        )


# ==================================================
# URUCHOMIENIE PROGRAMU
# ==================================================

if __name__ == "__main__":
    app = DNAApp()
    app.mainloop()
