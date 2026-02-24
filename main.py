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

        # ===== WYGLĄD / STYLE (TTK) =====
        self.style = ttk.Style(self)
        self.style.theme_use("clam")  # estetyczniejszy od defaultu na Windows

        # Stałe
        self.UI_FONT = ("Segoe UI", 10)
        self.UI_FONT_BOLD = ("Segoe UI", 10, "bold")
        self.UI_H1 = ("Segoe UI", 12, "bold")
        self.ACCENT = "#0b5394"  # ten granat co już lubisz

        # Ustawienia ogólne stylu
        self.style.configure(".", font=self.UI_FONT)
        self.style.configure("TNotebook.Tab", padding=(12, 6))
        self.style.configure("TLabel", font=self.UI_FONT)
        self.style.configure("Header.TLabel", font=self.UI_H1, foreground=self.ACCENT)

        # „Karta” – ramka z delikatnym tłem i obramowaniem
        self.style.configure("Card.TFrame", background="#f7f7f7", borderwidth=1, relief="solid")
        self.style.configure("CardInner.TFrame", background="#f7f7f7")

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

        self.active_label_x = None
        self.active_label_y = None
        self.selected_sequence = None
        self.selected_motif = None

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
        - zakładki (preview / wyniki / wizualizacja)
        - tabela wyników + tryb raw/norm
        - heatmapa + barplot + przycisk zamykania barplotu
        - logi (z przewijaniem)
        - pasek statusu (na samym dole)
        """

        # --- GÓRNA CZĘŚĆ (tabs) ---
        main_frame = tk.Frame(self)
        main_frame.pack(side="top", fill="both", expand=True)

        # ==================================================
        # PRAWA CZĘŚĆ Z ZAKŁADKAMI
        # ==================================================
        right = tk.Frame(main_frame)
        right.pack(fill="both", expand=True)

        self.tabs = ttk.Notebook(right)
        self.tabs.pack(fill="both", expand=True)

        # --- zakładki ---
        self.tab_preview = ttk.Frame(self.tabs)
        self.tab_results = ttk.Frame(self.tabs)
        self.tab_viz = ttk.Frame(self.tabs)

        self.tabs.add(self.tab_preview, text="Podgląd sekwencji")
        self.tabs.add(self.tab_results, text="Wyniki analizy")
        self.tabs.add(self.tab_viz, text="Wizualizacja")

        # ==================================================
        # PODGLĄD FASTA
        # ==================================================
        self.preview_box = tk.Text(self.tab_preview)
        self.preview_box.pack(fill="both", expand=True)

        # ==================================================
        # WYNIKI ANALIZY: RADIOBUTTONY + TABELA
        # ==================================================
        results_frame = tk.Frame(self.tab_results)
        results_frame.pack(fill="both", expand=True)

        # --- tryb wyświetlania (raw/norm) nad tabelą ---
        mode_frame = tk.Frame(results_frame)
        mode_frame.pack(fill="x", pady=(5, 0))

        tk.Label(mode_frame, text="Tryb wyświetlania:").pack(side="left", padx=(5, 5))

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

        # --- tabela wyników (Treeview) ---
        self.results_table = ttk.Treeview(results_frame, show="headings")
        self.results_table.pack(fill="both", expand=True)

        # --- UX: kursor "rączka" na nagłówkach tabeli (sortowanie) ---
        def on_heading_enter(event):
            region = self.results_table.identify_region(event.x, event.y)
            if region == "heading":
                self.results_table.config(cursor="hand2")
            else:
                self.results_table.config(cursor="")

        self.results_table.bind("<Motion>", on_heading_enter)

        # ==================================================
        # WIZUALIZACJA: RADIOBUTTONY NAD HEATMAPĄ + HEATMAPA + BARPLOT
        # ==================================================
        self.viz_frame = tk.Frame(self.tab_viz)
        self.viz_frame.pack(fill="both", expand=True)

        # --- panel sterowania (MUSI być w viz_frame, żeby był na górze zakładki) ---
        viz_controls = tk.Frame(self.viz_frame)
        viz_controls.pack(fill="x", pady=(5, 0))

        tk.Label(viz_controls, text="Tryb wyświetlania:").pack(side="left", padx=(10, 5))

        tk.Radiobutton(
            viz_controls,
            text="Surowe",
            variable=self.normalization_mode,
            value="raw",
            command=self.refresh_visualization
        ).pack(side="left")

        tk.Radiobutton(
            viz_controls,
            text="Na 1000 nt",
            variable=self.normalization_mode,
            value="norm",
            command=self.refresh_visualization
        ).pack(side="left")

        # --- górna część: heatmapa ---
        self.viz_top = tk.Frame(self.viz_frame)
        self.viz_top.pack(fill="both", expand=True)

        # --- dolna część: barplot (stała wysokość, żeby heatmapa nie znikała) ---
        self.viz_bottom = tk.Frame(self.viz_frame, height=260)
        self.viz_bottom.pack(fill="x")
        self.viz_bottom.pack_propagate(False)

        # --- przycisk zamykania barplotu (ukryty na starcie) ---
        self.close_barplot_btn = tk.Button(
            self.viz_frame,
            text="✖ Zamknij wykres",
            command=self.hide_barplot
        )
        self.close_barplot_btn.pack_forget()

        # ==================================================
        # DOLNY PANEL: LOGI + STATUS BAR
        # ==================================================

        # kontener dolny
        bottom_frame = tk.Frame(self)
        bottom_frame.pack(side="bottom", fill="x")

        # --- LOGI ---
        self.log_box = tk.Text(
            bottom_frame,
            height=6,          # stała wysokość
            bg="#f5f5f5"
        )
        self.log_box.pack(fill="x", padx=5, pady=(5, 0))

        # --- STATUS BAR ---
        self.status_var = tk.StringVar(value="Gotowe")

        status_bar = tk.Label(
            bottom_frame,
            textvariable=self.status_var,
            bd=1,
            relief="sunken",
            anchor="w"
        )
        status_bar.pack(fill="x")

    # ==================================================
    # FUNKCJE POMOCNICZE
    # ==================================================

    def log(self, message):
        """
        Dodaje komunikat do okna logów
        oraz aktualizuje pasek statusu.
        """

        # wpis do okna logów
        self.log_box.insert("end", message + "\n")
        self.log_box.see("end")  # auto-scroll na dół

        # wpis do status bara
        self.status_var.set(message)

    def update_status(self):
        seq_n = len(self.sequences) if self.sequences else 0
        mot_n = len(self.motifs) if self.motifs else 0
        mode = "Surowe" if self.normalization_mode.get() == "raw" else "Na 1000 nt"
        self.status_var.set(f"✔ Sekwencje: {seq_n} | Motywy: {mot_n} | Tryb: {mode}")

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

        # wyczyść podgląd
        self.preview_box.delete("1.0", "end")

        # wstaw podgląd sekwencji
        for k, v in self.sequences.items():
            self.preview_box.insert("end", f">{k}\n{v}\n\n")

        # automatycznie przelicz analizę
        self.run_analysis()

        self.log("Załadowano dane testowe (DEV_MODE)")
        self.update_status()

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
            self.update_status()

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
        self.update_status()

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

        # Czyścimy tylko górną część (heatmapę)
        for w in self.viz_top.winfo_children():
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
        self.heatmap_data = data

        # tworzymy wykres
        fig, ax = plt.subplots(figsize=(8, 5))

        # krótka instrukcja dla użytkownika
        ax.text(
            0.5,
            1.15,
            "Kliknij nazwę sekwencji (oś Y) lub motywu (oś X), aby wyświetlić wykres słupkowy.",
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=9,
            color="dimgray"
        )

        # ustalamy maksymalną skalę kolorów
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

        # ===== Tooltip (hover) =====
        self.hover_annotation = ax.annotate(
            "",
            xy=(0, 0),
            xytext=(10, 10),
            textcoords="offset points",
            bbox=dict(boxstyle="round", fc="white", ec="gray"),
            fontsize=9
        )
        self.hover_annotation.set_visible(False)

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

        # etykiety osi x
        ax.set_xticks(np.arange(len(self.motifs)))
        ax.set_xticklabels(self.motifs, rotation=45, ha="right")

        # etykiety osi y
        ax.set_yticks(np.arange(len(self.sequences)))
        ax.set_yticklabels(list(self.sequences.keys()))

        # Umożliwiamy klikanie w etykiety osi
        for label in ax.get_xticklabels():
            label.set_picker(5)

        for label in ax.get_yticklabels():
            label.set_picker(5)

        # colorbar
        cbar = fig.colorbar(im)
        cbar.set_label("Liczba motywów" if self.normalization_mode.get() == "raw" else "Na 1000 nt")

        ax.set_title(
            "Heatmapa wystąpień motywów",
            fontsize=14,
            fontweight="bold",
            pad=15
        )

        fig.tight_layout()

        # osadzamy w Tkinter
        canvas = FigureCanvasTkAgg(fig, master=self.viz_top)
        canvas.draw()
        self.current_canvas = canvas
        canvas.mpl_connect("pick_event", self.on_pick)
        canvas.mpl_connect("motion_notify_event", self.on_hover)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack(fill="both", expand=True)

    def on_pick(self, event):
        if not hasattr(event.artist, "get_text"):
            return

        label = event.artist
        text = label.get_text()

        # Reset poprzednich podświetleń
        if self.active_label_x:
            self.active_label_x.set_color("black")
            self.active_label_x.set_fontweight("normal")

        if self.active_label_y:
            self.active_label_y.set_color("black")
            self.active_label_y.set_fontweight("normal")

        # Kliknięto sekwencję (oś Y)
        if text in self.sequences:
            self.selected_sequence = text
            self.draw_barplot_sequence()

            label.set_color("red")
            label.set_fontweight("bold")
            self.active_label_y = label
            self.active_label_x = None

        # Kliknięto motyw (oś X)
        elif text in self.motifs:
            self.selected_motif = text
            self.draw_barplot_motif()

            label.set_color("red")
            label.set_fontweight("bold")
            self.active_label_x = label
            self.active_label_y = None

        # 🔥 kluczowe do podświetlenia
        event.canvas.draw_idle()

    def draw_barplot_sequence(self):
        """
        Rysuje wykres słupkowy dla aktualnie wybranej sekwencji.

        Wykres pokazuje:
        - liczbę wystąpień każdego motywu
        - lub wartości znormalizowane (na 1000 nt),
          w zależności od wybranego trybu (raw / norm).
        """

        # Jeśli nie wybrano sekwencji (np. brak kliknięcia w heatmapę)
        if not self.selected_sequence:
            return

        # Pobieramy sekwencję DNA na podstawie zapisanego ID
        seq = self.sequences[self.selected_sequence]
        seq_length = len(seq)

        # Lista wartości dla osi Y (liczby motywów)
        values = []

        # Dla każdego motywu liczymy jego wystąpienia
        for motif in self.motifs:
            count = iupac_count(seq, motif)

            # Jeśli wybrano tryb normalizacji
            # przeliczamy na 1000 nukleotydów
            if self.normalization_mode.get() == "norm" and seq_length > 0:
                value = (count / seq_length) * 1000
            else:
                value = count

            values.append(value)

        # Tworzymy nową figurę matplotlib
        fig, ax = plt.subplots(figsize=(6, 4))

        # Główny tytuł wykresu (spójny z heatmapą)
        fig.suptitle(
            "Rozkład motywów w sekwencji",
            fontsize=14,
            fontweight="bold"
        )

        # Pozycje słupków
        x_positions = np.arange(len(self.motifs))

        # Rysowanie słupków
        ax.bar(x_positions, values)

        ax.set_xticks(x_positions)
        ax.set_xticklabels(self.motifs, rotation=45, ha="right")

        # Podtytuł – konkretna sekwencja
        ax.set_title(
            f"Sekwencja: {self.selected_sequence}",
            fontsize=11,
            fontweight="normal",
            pad=5
        )

        # Opis osi Y
        ax.set_ylabel(
            "Na 1000 nt"
            if self.normalization_mode.get() == "norm"
            else "Liczba"
        )

        fig.tight_layout(rect=[0, 0, 1, 0.92])
        fig.subplots_adjust(bottom=0.30)

        # Czyścimy poprzedni wykres słupkowy
        for w in self.viz_bottom.winfo_children():
            w.destroy()

        self.viz_bottom.configure(height=260)
        self.viz_bottom.pack_propagate(False)

        # Osadzamy wykres w Tkinter
        canvas = FigureCanvasTkAgg(fig, master=self.viz_bottom)
        canvas.draw()

        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack(fill="both", expand=True, padx=10, pady=10)

        # pokaż panel i przycisk
        self.viz_bottom.pack(fill="x")
        self.close_barplot_btn.pack(pady=(5, 0))

    def draw_barplot_motif(self):
        """Rysuje wykres słupkowy dla wybranego motywu (rozkład w sekwencjach)."""

        if not self.selected_motif:
            return

        values = []

        for seq_id, seq in self.sequences.items():
            seq_length = len(seq)
            count = iupac_count(seq, self.selected_motif)

            if self.normalization_mode.get() == "norm" and seq_length > 0:
                value = (count / seq_length) * 1000
            else:
                value = count

            values.append(value)

        # Czyścimy dolny panel
        for w in self.viz_bottom.winfo_children():
            w.destroy()

        fig, ax = plt.subplots(figsize=(6, 4))

        # Główny tytuł wykresu (taki sam styl jak heatmapa)
        fig.suptitle(
            "Rozkład motywu w sekwencjach",
            fontsize=14,
            fontweight="bold"
        )

        x_positions = np.arange(len(self.sequences))
        ax.bar(x_positions, values)

        ax.set_xticks(x_positions)
        ax.set_xticklabels(list(self.sequences.keys()), rotation=45, ha="right")

        # Podtytuł – konkretna nazwa motywu
        ax.set_title(
            f"Motyw: {self.selected_motif}",
            fontsize=11,
            fontweight="normal",
            pad=5
        )

        ax.set_ylabel(
            "Na 1000 nt"
            if self.normalization_mode.get() == "norm"
            else "Liczba"
        )

        fig.tight_layout(rect=[0, 0, 1, 0.92])
        fig.subplots_adjust(bottom=0.30)

        self.viz_bottom.configure(height=260)
        self.viz_bottom.pack_propagate(False)

        canvas = FigureCanvasTkAgg(fig, master=self.viz_bottom)
        canvas.draw()

        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack(fill="both", expand=True, padx=10, pady=10)

        # pokaż przycisk (panel zostaje w layoucie cały czas)
        self.close_barplot_btn.pack(pady=(5, 0))

    def on_hover(self, event):
        # jeśli poza wykresem → nic nie rób
        if event.inaxes is None:
            if hasattr(self, "hover_annotation"):
                self.hover_annotation.set_visible(False)
                event.canvas.draw_idle()
            return

        if not hasattr(self, "heatmap_data"):
            return

        x = int(round(event.xdata))
        y = int(round(event.ydata))

        if 0 <= x < len(self.motifs) and 0 <= y < len(self.sequences):
            seq_id = list(self.sequences.keys())[y]
            motif = self.motifs[x]
            value = self.heatmap_data[y, x]

            text = f"{seq_id}\n{motif}"

            self.hover_annotation.xy = (x, y)
            self.hover_annotation.set_text(text)
            self.hover_annotation.set_visible(True)
            event.canvas.draw_idle()
        else:
            self.hover_annotation.set_visible(False)
            event.canvas.draw_idle()

    def refresh_visualization(self):
        """
        Odświeża heatmapę i barplot po zmianie trybu raw / norm.
        """

        self.draw_visualization()

        # jeśli coś było wybrane – przerysuj barplot
        if hasattr(self, "selected_sequence") and self.selected_sequence:
            self.draw_barplot_sequence()

        if hasattr(self, "selected_motif") and self.selected_motif:
            self.draw_barplot_motif()

        self.update_status()

    def hide_barplot(self):
        """
        Ukrywa wykres słupkowy (czyści zawartość),
        resetuje zaznaczenia etykiet.
        """

        # 1) wyczyść dolny panel
        for w in self.viz_bottom.winfo_children():
            w.destroy()

        # 2) ZAMIANA KLUCZOWA:
        #    NIE pack_forget() -> zostawiamy viz_bottom w layoucie
        #    żeby nie rozjeżdżało układu (logi/status).
        self.viz_bottom.configure(height=0)
        self.viz_bottom.pack_propagate(False)

        # 3) ukryj przycisk
        self.close_barplot_btn.pack_forget()

        # 4) reset podświetleń etykiet
        if getattr(self, "active_label_x", None):
            self.active_label_x.set_color("black")
            self.active_label_x.set_fontweight("normal")
            self.active_label_x = None

        if getattr(self, "active_label_y", None):
            self.active_label_y.set_color("black")
            self.active_label_y.set_fontweight("normal")
            self.active_label_y = None

        # 5) reset wyboru
        self.selected_sequence = None
        self.selected_motif = None

        # 6) odśwież canvas (żeby zmiana koloru była widoczna)
        if hasattr(self, "current_canvas"):
            self.current_canvas.draw_idle()

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
