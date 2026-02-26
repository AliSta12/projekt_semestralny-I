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
from analysis_engine import compute_analysis
from export_tab import ExportTab


# ==================================================
# LOGIKA – PRZETWARZANIE SEKWENCJI DNA
# ==================================================

def load_fasta(path):
    """
    Wczytuje plik FASTA i zwraca słownik {ID: sekwencja}.

    Polityka zgodna z rekomendacją:
    - FULL IUPAC DNA: A C G T R Y S W K M B D H V N
    - U -> T
    - usuwa alignment gapy/placeholdery: '-' i '.'
    - '?' (nieznana baza) -> 'N'  (zachowuje "dziurę" bez sklejania)
    - usuwa spacje, taby, cyfry (formatowanie)
    - inne znaki -> błąd
    """

    sequences = {}
    current_id = None

    allowed = set("ACGTRYSWKMBDHVN")

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                current_id = line[1:].strip()
                if not current_id:
                    raise ValueError(f"Pusty nagłówek w linii {line_num}")
                if current_id in sequences:
                    raise ValueError(f"Duplikat ID: {current_id}")
                sequences[current_id] = ""
                continue

            if current_id is None:
                raise ValueError("Plik nie zaczyna się od nagłówka FASTA")

            seq_line = line.upper()

            # RNA -> DNA
            seq_line = seq_line.replace("U", "T")

            # alignment gapy/placeholdery: usuwamy
            seq_line = seq_line.replace("-", "").replace(".", "")

            # nieznane znaki: mapujemy na N (uczciwsze niż usuwanie)
            seq_line = seq_line.replace("?", "N")

            # formatowanie: usuwamy
            seq_line = seq_line.replace(" ", "").replace("\t", "")
            seq_line = "".join(c for c in seq_line if not c.isdigit())

            bad = {c for c in seq_line if c not in allowed}
            if bad:
                raise ValueError(f"Niepoprawne znaki {sorted(bad)} w linii {line_num}")

            sequences[current_id] += seq_line

    if not sequences:
        raise ValueError("Nie znaleziono sekwencji")

    return sequences

# ==================================================
# IUPAC – MAPA (bitmask) + DOPASOWANIE
# ==================================================

_IUPAC_MASK = {
    "A": 0b0001,
    "C": 0b0010,
    "G": 0b0100,
    "T": 0b1000,
    "U": 0b1000,
    "R": 0b0101,
    "Y": 0b1010,
    "S": 0b0110,
    "W": 0b1001,
    "K": 0b1100,
    "M": 0b0011,
    "B": 0b1110,
    "D": 0b1101,
    "H": 0b1011,
    "V": 0b0111,
    "N": 0b1111,
}


def matches_iupac(seq_char: str, motif_char: str) -> bool:
    """
    Sprawdza, czy znak sekwencji i motywu mają
    część wspólną w sensie IUPAC.
    """
    try:
        return (_IUPAC_MASK[seq_char.upper()] &
                _IUPAC_MASK[motif_char.upper()]) != 0
    except KeyError:
        return False


def iupac_find_positions(sequence: str, motif: str) -> list[int]:
    """
    Znajduje pozycje motywu (IUPAC vs IUPAC).
    Zwraca listę pozycji (0-based).
    """
    seq = sequence.upper().replace("U", "T")
    mot = motif.upper().replace("U", "T")

    positions = []
    L = len(seq)
    mL = len(mot)

    for i in range(L - mL + 1):
        match = True
        for j in range(mL):
            if not matches_iupac(seq[i + j], mot[j]):
                match = False
                break
        if match:
            positions.append(i)

    return positions

def iupac_count(sequence: str, motif: str) -> int:
    """
    Liczy wystąpienia motywu w sekwencji z uwzględnieniem IUPAC po obu stronach.
    Liczy z nakładaniem (overlapping).
    """
    if not sequence or not motif:
        return 0

    seq = sequence.upper().replace("U", "T")
    mot = motif.upper().replace("U", "T")

    # przygotuj maski; jeśli trafisz na nieznany znak -> błąd (łatwiej debugować)
    try:
        seq_masks = [_IUPAC_MASK[ch] for ch in seq]
        mot_masks = [_IUPAC_MASK[ch] for ch in mot]
    except KeyError as e:
        raise ValueError(f"Nieznany symbol IUPAC: {e.args[0]}")

    m = len(mot_masks)
    n = len(seq_masks)
    if m > n:
        return 0

    count = 0
    for i in range(n - m + 1):
        ok = True
        for j in range(m):
            # pasuje, jeśli zbiory możliwych zasad mają część wspólną
            if (seq_masks[i + j] & mot_masks[j]) == 0:
                ok = False
                break
        if ok:
            count += 1

    return count

def find_motif_positions(sequence: str, motif: str) -> list[int]:
    """Pozycje motywu z obsługą IUPAC (overlap)."""
    return iupac_find_positions(sequence, motif)

def count_motif(sequence: str, motif: str) -> int:
    """Liczba wystąpień motywu z obsługą IUPAC (overlap)."""
    return iupac_count(sequence, motif)

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
        self.DEV_MODE = False  # <-- zmień na False gdy niepotrzebne
        self.analysis_result = None
        self.figures = {}

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
        motywy.add_command(label="Dodaj/usuń motywy", command=self.open_motif_manager)

        ncbi = tk.Menu(menubar, tearoff=0)
        ncbi.add_command(label="Pobierz z NCBI", command=self.download_ncbi)

        pomoc = tk.Menu(menubar, tearoff=0)
        pomoc.add_command(label="O programie", command=self.about)

        menubar.add_cascade(label="Plik", menu=plik)
        menubar.add_cascade(label="Motywy", menu=motywy)
        menubar.add_cascade(label="NCBI", menu=ncbi)
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

        # ==================================================
        # UKŁAD: 2 WIERSZE (GRID) -> dół zawsze widoczny
        # ==================================================

        # kasujemy ustawienia pack/grid z root (bezpieczne przy restartach layoutu)
        for w in self.winfo_children():
            w.pack_forget()
            w.grid_forget()

        self.grid_rowconfigure(0, weight=1)  # zakładki rosną
        self.grid_rowconfigure(1, weight=0)  # dół ma stałe miejsce
        self.grid_columnconfigure(0, weight=1)

        # --- GÓRA: zakładki ---
        main_frame = tk.Frame(self)
        main_frame.grid(row=0, column=0, sticky="nsew")

        # --- DÓŁ: logi + status (STAŁA wysokość) ---
        bottom_frame = tk.Frame(self, height=140)
        bottom_frame.grid(row=1, column=0, sticky="ew")
        bottom_frame.grid_propagate(False)  # nie pozwól, żeby dół się “zwinął”

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
        self.tab_export = ttk.Frame(self.tabs)
        self.tabs.add(self.tab_export, text="Eksport")

        export_widget = ExportTab(
            self.tab_export,
            get_result=lambda: self.analysis_result,
            get_figures=lambda: self.figures,
            log_fn=self.log,
        )
        export_widget.pack(fill="both", expand=True)

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

        # --- STATUS BAR (zawsze na dole) ---
        self.status_var = tk.StringVar(value="Gotowe")
        status_bar = tk.Label(
            bottom_frame,
            textvariable=self.status_var,
            bd=1,
            relief="sunken",
            anchor="w"
        )
        status_bar.pack(side="bottom", fill="x")

        # --- LOGI (nad statusem) ---
        log_frame = tk.Frame(bottom_frame)
        log_frame.pack(side="top", fill="both", expand=True)

        tk.Label(log_frame, text="Logi / komunikaty:").pack(anchor="w", padx=6)

        self.log_box = tk.Text(log_frame, height=3)
        self.log_box.pack(fill="both", expand=True, padx=6, pady=(0, 6))

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

            self.analysis_result = None
            self.figures.clear()

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
        - dodawanie ręczne (Enter + przycisk)
        - lista aktualnych motywów z usuwaniem
        - import / eksport / wyczyść
        - akordeon (opcjonalne motywy referencyjne) z checkboxami
        - okno automatycznie dopasowuje wysokość do zawartości (żeby nic nie ucinać)
        """

        import tkinter as tk
        from tkinter import ttk, filedialog, messagebox
        from pathlib import Path

        IUPAC_ALLOWED = set("ACGTURYSWKMBDHVN")

        def normalize_motif(s: str) -> str:
            return "".join(s.split()).upper()

        def validate_motif(m: str) -> tuple[bool, str]:
            if not m:
                return False, "Motyw jest pusty."
            bad = sorted({ch for ch in m if ch not in IUPAC_ALLOWED})
            if bad:
                return False, f"Nieprawidłowe znaki: {', '.join(bad)}.\nDozwolone: {''.join(sorted(IUPAC_ALLOWED))}"
            if len(m) < 2:
                return False, "Motyw jest za krótki (min. 2 znaki)."
            return True, ""

        def unique_preserve_order(items: list[str]) -> list[str]:
            seen = set()
            out = []
            for x in items:
                if x not in seen:
                    seen.add(x)
                    out.append(x)
            return out

        def pl_motyw(n: int) -> str:
            # 1 motyw, 2-4 motywy, 5+ motywów; wyjątek 12-14 motywów
            if n == 1:
                return "motyw"
            if 12 <= (n % 100) <= 14:
                return "motywów"
            if 2 <= (n % 10) <= 4:
                return "motywy"
            return "motywów"

        def sync_to_main_listbox():
            if hasattr(self, "motif_listbox") and self.motif_listbox:
                try:
                    self.motif_listbox.delete(0, "end")
                    for m in self.motifs:
                        self.motif_listbox.insert("end", m)
                except Exception:
                    pass

        # ===== Window =====
        win = tk.Toplevel(self)
        win.title("Motywy")
        win.transient(self)
        win.grab_set()
        win.resizable(False, True)  # stała szerokość, wysokość dopasowujemy

        initial_width = 460
        win.geometry(f"{initial_width}x520")

        # autosize: dopasuj wysokość do reqheight (z limitem ekranu)
        def autosize():
            win.update_idletasks()
            req_h = win.winfo_reqheight()
            screen_h = win.winfo_screenheight()
            max_h = max(360, screen_h - 120)
            new_h = min(req_h, max_h)
            win.geometry(f"{initial_width}x{new_h}")

        def on_close():
            # nie czyścimy motywów przy zamknięciu okna (Twoje wymaganie)
            win.destroy()

        win.protocol("WM_DELETE_WINDOW", on_close)

        # ===== Layout root =====
        root = ttk.Frame(win, padding=10)
        root.pack(fill="both", expand=True)
        root.columnconfigure(0, weight=1)
        root.rowconfigure(5, weight=1)  # lista motywów rośnie

        # ===== Legenda + tooltip =====
        iupac_text = "A C G T + IUPAC (R Y S W K M B D H V N)"
        iupac_lbl = ttk.Label(root, text=iupac_text, foreground="#0b5394")
        iupac_lbl.grid(row=1, column=0, sticky="w", pady=(0, 6))

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

        top_row = ttk.Frame(root)
        top_row.grid(row=0, column=0, sticky="ew", pady=(0, 2))
        top_row.columnconfigure(0, weight=1)

        ttk.Label(top_row, text="Dozwolone symbole:", font=("Segoe UI", 9, "bold")).grid(row=0, column=0, sticky="w")

        btn_run = ttk.Button(top_row, text="Uruchom analizę", command=self.run_analysis)
        btn_run.grid(row=0, column=1, sticky="e")

        def run_and_close():
            # zamknij okno motywów
            try:
                win.destroy()
            except Exception:
                pass

            # przełącz na wyniki
            try:
                self.tabs.select(self.tab_results)
            except Exception:
                pass

            # uruchom analizę po chwili, żeby UI zdążyło się przełączyć
            self.after(50, self.run_analysis)

        btn_run.configure(command=run_and_close)

        def show_tooltip(event):
            nonlocal tooltip
            if tooltip:
                return
            x, y = event.x_root + 12, event.y_root + 12
            tooltip = tk.Toplevel(win)
            tooltip.wm_overrideredirect(True)
            tooltip.geometry(f"+{x}+{y}")
            tk.Label(
                tooltip, text=tooltip_text, justify="left", bg="#ffffe0",
                relief="solid", borderwidth=1, padx=6, pady=4
            ).pack()

        def hide_tooltip(event):
            nonlocal tooltip
            if tooltip:
                tooltip.destroy()
                tooltip = None

        iupac_lbl.bind("<Enter>", show_tooltip)
        iupac_lbl.bind("<Leave>", hide_tooltip)

        # ===== Add motif =====
        ttk.Label(root, text="Dodaj motyw:", font=("Segoe UI", 9, "bold")).grid(row=2, column=0, sticky="w",
                                                                                pady=(4, 2))
        add_row = ttk.Frame(root)
        add_row.grid(row=3, column=0, sticky="ew")
        add_row.columnconfigure(0, weight=1)

        motif_entry = ttk.Entry(add_row)
        motif_entry.grid(row=0, column=0, sticky="ew")

        # ===== Motifs list =====
        ttk.Label(root, text="Aktualne motywy:", font=("Segoe UI", 9, "bold")).grid(row=4, column=0, sticky="w",
                                                                                    pady=(10, 2))

        list_frame = ttk.Frame(root)
        list_frame.grid(row=5, column=0, sticky="nsew")
        list_frame.columnconfigure(0, weight=1)
        list_frame.rowconfigure(0, weight=1)

        motifs_lb = tk.Listbox(list_frame, selectmode="extended", height=8, activestyle="none")
        motifs_lb.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(list_frame, orient="vertical", command=motifs_lb.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        motifs_lb.configure(yscrollcommand=scroll.set)

        # ===== Actions (przyciski + licznik w osobnym wierszu) =====
        actions = ttk.Frame(root)
        actions.grid(row=6, column=0, sticky="ew", pady=(10, 6))
        for c in range(4):
            actions.columnconfigure(c, weight=1)

        badge_var = tk.StringVar(value="0 motywów")
        badge = ttk.Label(actions, textvariable=badge_var, anchor="e")
        badge.grid(row=1, column=0, columnspan=4, sticky="e", pady=(6, 0))

        # ---- state/refresh ----
        def refresh_lists():
            motifs_lb.delete(0, "end")
            for m in self.motifs:
                motifs_lb.insert("end", m)

            n = len(self.motifs)
            badge_var.set(f"{n} {pl_motyw(n)}")

            has_any = bool(self.motifs)
            btn_export.configure(state=("normal" if has_any else "disabled"))
            btn_clear.configure(state=("normal" if has_any else "disabled"))

            # stan "usuń" zależy od selekcji
            btn_remove.configure(state=("normal" if motifs_lb.curselection() else "disabled"))

            sync_to_main_listbox()

            # przycisk analizy w oknie motywów
            btn_run.configure(state=("normal" if (self.sequences and self.motifs) else "disabled"))

            autosize()

        def add_motif(m: str):
            m = normalize_motif(m)
            ok, reason = validate_motif(m)
            if not ok:
                messagebox.showwarning("Nieprawidłowy motyw", reason, parent=win)
                return

            if m in self.motifs:
                self.log(f"Motyw już istnieje: {m}")
                motif_entry.delete(0, "end")
                return

            self.motifs.append(m)
            self.motifs = unique_preserve_order(self.motifs)
            self.analysis_result = None
            self.figures.clear()
            motif_entry.delete(0, "end")
            self.log(f"Dodano motyw: {m}")
            refresh_lists()
            refresh_accordion_checks()

        def remove_selected():
            sel = motifs_lb.curselection()
            if not sel:
                return
            indices = sorted(sel, reverse=True)
            removed = []
            for idx in indices:
                try:
                    removed.append(self.motifs.pop(idx))
                except Exception:
                    pass
            if removed:
                self.log(f"Usunięto: {', '.join(removed)}")
            refresh_lists()
            refresh_accordion_checks()

        def clear_all():
            if not self.motifs:
                return
            if not messagebox.askyesno("Wyczyścić motywy?", "Na pewno usunąć wszystkie motywy?", parent=win):
                return
            self.motifs.clear()
            self.analysis_result = None
            self.figures.clear()
            self.log("Wyczyszczono wszystkie motywy.")
            refresh_lists()
            refresh_accordion_checks()

        def import_motifs():
            path = filedialog.askopenfilename(
                parent=win,
                title="Import motywów",
                filetypes=[("CSV", "*.csv"), ("Tekst", "*.txt"), ("Wszystkie", "*.*")]
            )
            if not path:
                return

            p = Path(path)

            # ========================
            # IMPORT CSV Z WYBOREM KOLUMNY
            # ========================
            if p.suffix.lower() == ".csv":
                import csv

                with p.open(encoding="utf-8", errors="ignore") as f:
                    reader = csv.DictReader(f)
                    headers = reader.fieldnames

                    if not headers:
                        messagebox.showerror("Błąd", "Nie znaleziono nagłówków w pliku CSV.", parent=win)
                        return

                    # --- okno wyboru kolumny ---
                    col_win = tk.Toplevel(win)
                    col_win.title("Wybierz kolumnę z motywami")
                    col_win.transient(win)
                    col_win.grab_set()
                    col_win.resizable(False, False)

                    ttk.Label(col_win, text="Wybierz kolumnę:").pack(padx=10, pady=(10, 4))

                    selected_col = tk.StringVar(value=headers[0])
                    combo = ttk.Combobox(col_win, values=headers, textvariable=selected_col, state="readonly")
                    combo.pack(padx=10, pady=4)

                    def confirm_column():
                        col = selected_col.get()
                        col_win.destroy()

                        # wczytaj ponownie plik
                        with p.open(encoding="utf-8", errors="ignore") as f2:
                            reader2 = csv.DictReader(f2)
                            tokens = []
                            for row in reader2:
                                value = row.get(col)
                                if value:
                                    tokens.append(value)

                        add_imported_tokens(tokens)

                    ttk.Button(col_win, text="Importuj", command=confirm_column).pack(pady=(4, 10))

                return

            # ========================
            # IMPORT TXT (jak wcześniej)
            # ========================
            text = p.read_text(encoding="utf-8", errors="ignore")

            tokens = []
            for line in text.splitlines():
                line = line.strip()
                if not line:
                    continue
                for part in line.replace(";", ",").replace("\t", ",").split(","):
                    part = part.strip()
                    if part:
                        tokens.append(part)

            add_imported_tokens(tokens)

        def add_imported_tokens(tokens: list[str]):
            before = len(self.motifs)
            added = 0

            for t in tokens:
                m = normalize_motif(t)
                ok, _ = validate_motif(m)
                if ok and m not in self.motifs:
                    self.motifs.append(m)
                    added += 1

            self.motifs = unique_preserve_order(self.motifs)
            if added > 0:
                self.analysis_result = None
                self.figures.clear()

            self.log(f"Import CSV/TXT (dodano {added}, było {before}, jest {len(self.motifs)})")

            refresh_lists()
            refresh_accordion_checks()

        def export_motifs():
            path = filedialog.asksaveasfilename(
                parent=win,
                title="Eksport motywów",
                defaultextension=".txt",
                filetypes=[("Tekst", "*.txt"), ("CSV", "*.csv")]
            )
            if not path:
                return
            p = Path(path)
            if p.suffix.lower() == ".csv":
                content = "motif\n" + "\n".join(self.motifs) + "\n"
            else:
                content = "\n".join(self.motifs) + "\n"
            p.write_text(content, encoding="utf-8")
            self.log(f"Eksport: {p.name} ({len(self.motifs)} {pl_motyw(len(self.motifs))})")

        # ===== Podpinanie akcji do widgetów (po definicjach funkcji!) =====
        motif_entry.bind("<Return>", lambda e: add_motif(motif_entry.get()))

        btn_add = ttk.Button(add_row, text="Dodaj", command=lambda: add_motif(motif_entry.get()))
        btn_add.grid(row=0, column=1, padx=(8, 0))

        btn_remove = ttk.Button(actions, text="Usuń zaznaczone", command=remove_selected)
        btn_remove.grid(row=0, column=0, sticky="ew", padx=(0, 6))

        btn_clear = ttk.Button(actions, text="Wyczyść motywy", command=clear_all)
        btn_clear.grid(row=0, column=1, sticky="ew", padx=(0, 6))

        btn_import = ttk.Button(actions, text="Import", command=import_motifs)
        btn_import.grid(row=0, column=2, sticky="ew", padx=(0, 6))

        btn_export = ttk.Button(actions, text="Eksport", command=export_motifs)
        btn_export.grid(row=0, column=3, sticky="ew")

        def on_list_select(event=None):
            btn_remove.configure(state=("normal" if motifs_lb.curselection() else "disabled"))

        motifs_lb.bind("<<ListboxSelect>>", on_list_select)

        # ===== Accordion =====
        ttk.Separator(root).grid(row=7, column=0, sticky="ew", pady=(8, 8))
        ttk.Label(root, text="Motywy referencyjne (opcjonalnie):", font=("Segoe UI", 9, "bold")).grid(
            row=8, column=0, sticky="w"
        )

        accordion_container = ttk.Frame(root)
        accordion_container.grid(row=9, column=0, sticky="ew", pady=(6, 0))
        accordion_container.columnconfigure(0, weight=1)

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

        accordion_frames: dict[str, ttk.Frame] = {}
        icon_vars = {
            "Eukariota": tk.StringVar(value="▸  Eukariota"),
            "Prokariota": tk.StringVar(value="▸  Prokariota"),
        }
        section_vars: dict[str, list[tk.BooleanVar]] = {"Eukariota": [], "Prokariota": []}

        def toggle_section(name: str):
            frm = accordion_frames[name]
            if frm.winfo_ismapped():
                frm.grid_remove()
                icon_vars[name].set(f"▸  {name}")
                autosize()
                return

            for nm, fr in accordion_frames.items():
                fr.grid_remove()
                icon_vars[nm].set(f"▸  {nm}")

            frm.grid()
            icon_vars[name].set(f"▾  {name}")
            autosize()

        def on_toggle_motif(seq: str, checked: bool):
            changed = False

            if checked:
                if seq not in self.motifs:
                    self.motifs.append(seq)
                    self.motifs = unique_preserve_order(self.motifs)
                    self.log(f"Dodano motyw: {seq}")
                    changed = True
            else:
                if seq in self.motifs:
                    self.motifs = [m for m in self.motifs if m != seq]
                    self.log(f"Usunięto motyw: {seq}")
                    changed = True

            if changed:
                self.analysis_result = None
                self.figures.clear()

            refresh_lists()
            refresh_accordion_checks()

        def build_section(row: int, name: str, motifs: list[tuple[str, str]]):
            btn = ttk.Button(
                accordion_container,
                textvariable=icon_vars[name],
                command=lambda: toggle_section(name)
            )
            btn.grid(row=row, column=0, sticky="ew", pady=(0, 4))

            frm = ttk.Frame(accordion_container)
            frm.grid(row=row + 1, column=0, sticky="ew")
            frm.grid_remove()
            accordion_frames[name] = frm

            for seq, desc in motifs:
                var = tk.BooleanVar(master=win, value=(seq in self.motifs))
                section_vars[name].append(var)
                cb = ttk.Checkbutton(
                    frm,
                    text=f"{seq} — {desc}",
                    variable=var,
                    command=lambda s=seq, v=var: on_toggle_motif(s, v.get())
                )
                cb.pack(anchor="w")

        def refresh_accordion_checks():
            for name, motifs in [("Eukariota", EU_MOTIFS), ("Prokariota", PRO_MOTIFS)]:
                vars_list = section_vars[name]
                for (seq, _), var in zip(motifs, vars_list):
                    var.set(seq in self.motifs)

        build_section(0, "Eukariota", EU_MOTIFS)
        build_section(2, "Prokariota", PRO_MOTIFS)

        # ===== Init =====
        self.motifs = [normalize_motif(m) for m in getattr(self, "motifs", []) if normalize_motif(m)]
        self.motifs = unique_preserve_order(self.motifs)

        refresh_lists()
        refresh_accordion_checks()
        on_list_select()
        autosize()

    def run_analysis(self):
        if not self.sequences or not self.motifs:
            messagebox.showwarning("Błąd", "Wczytaj FASTA i dodaj motywy")
            return

        try:
            # 🔬 liczenie tylko raz
            self.analysis_result = compute_analysis(
                sequences=self.sequences,
                motifs=self.motifs,
                iupac_count_fn=iupac_count,
                fasta_path=self.current_fasta,
            )

            # 🔁 render widoków
            self.render_results_table()
            self.draw_visualization()

            self.update_status()
            self.log("Analiza zakończona")

        except Exception as e:
            messagebox.showerror("Błąd analizy", str(e))

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

        if self.analysis_result is None:
            return

        data = self.analysis_result.matrix(self.normalization_mode.get())
        self.heatmap_data = data

        # tworzymy wykres
        fig, ax = plt.subplots(figsize=(8, 5))
        self.figures["heatmap"] = fig

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

        motifs = self.analysis_result.motifs
        seq_ids = self.analysis_result.seq_ids

        ax.set_xticks(np.arange(len(motifs)))
        ax.set_xticklabels(motifs, rotation=45, ha="right")

        ax.set_yticks(np.arange(len(seq_ids)))
        ax.set_yticklabels(seq_ids)

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

        if self.analysis_result is None:
            return

        mode = self.normalization_mode.get()
        mat = self.analysis_result.matrix(mode)

        # znajdź indeks wybranej sekwencji
        i = self.analysis_result.seq_ids.index(self.selected_sequence)

        # pobierz wartości dla wszystkich motywów
        values = mat[i, :]

        # Tworzymy nową figurę matplotlib
        fig, ax = plt.subplots(figsize=(6, 4))

        self.figures[f"barplot_seq_{self.selected_sequence}"] = fig

        # Główny tytuł wykresu (spójny z heatmapą)
        fig.suptitle(
            "Rozkład motywów w sekwencji",
            fontsize=14,
            fontweight="bold"
        )

        motifs = self.analysis_result.motifs

        x_positions = np.arange(len(motifs))
        ax.bar(x_positions, values)

        ax.set_xticks(x_positions)
        ax.set_xticklabels(motifs, rotation=45, ha="right")

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

        if self.analysis_result is None:
            return

        mode = self.normalization_mode.get()
        mat = self.analysis_result.matrix(mode)

        # indeks wybranego motywu
        j = self.analysis_result.motifs.index(self.selected_motif)

        # wartości dla wszystkich sekwencji
        values = mat[:, j]

        # Czyścimy dolny panel
        for w in self.viz_bottom.winfo_children():
            w.destroy()

        fig, ax = plt.subplots(figsize=(6, 4))

        self.figures[f"barplot_motif_{self.selected_motif}"] = fig

        # Główny tytuł wykresu (taki sam styl jak heatmapa)
        fig.suptitle(
            "Rozkład motywu w sekwencjach",
            fontsize=14,
            fontweight="bold"
        )

        seq_ids = self.analysis_result.seq_ids

        x_positions = np.arange(len(seq_ids))
        ax.bar(x_positions, values)

        ax.set_xticks(x_positions)
        ax.set_xticklabels(seq_ids, rotation=45, ha="right")

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
        self.render_results_table()
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


    def download_ncbi(self):
        """Pobierz FASTA z NCBI (E-utilities) i wczytaj do aplikacji."""
        import tkinter as tk
        from tkinter import ttk, messagebox, filedialog
        from urllib.parse import urlencode
        from urllib.request import urlopen, Request
        import re
        from pathlib import Path

        def http_get(url: str) -> str:
            req = Request(url, headers={"User-Agent": "DNAApp/1.0 (NCBI E-utilities)"})
            with urlopen(req, timeout=30) as r:
                return r.read().decode("utf-8", errors="replace")

        def is_accession_like(s: str) -> bool:
            # liberalnie: bez spacji, tylko znaki typowe dla accession/UID
            return bool(re.fullmatch(r"[A-Za-z0-9_.]+", s))

        def parse_esearch_ids(xml_text: str) -> list[str]:
            return re.findall(r"<Id>(\d+)</Id>", xml_text)

        win = tk.Toplevel(self)
        win.title("Pobierz FASTA z NCBI")
        win.transient(self)
        win.grab_set()
        win.resizable(False, False)

        frm = ttk.Frame(win, padding=12)
        frm.pack(fill="both", expand=True)

        ttk.Label(frm, text="Wpisz accession/UID (np. NC_000001.11) lub zapytanie (ESearch):",
                  font=("Segoe UI", 9, "bold")).grid(row=0, column=0, columnspan=3, sticky="w")

        entry = ttk.Entry(frm, width=56)
        entry.grid(row=1, column=0, columnspan=3, sticky="ew", pady=(4, 8))
        entry.focus_set()

        ttk.Label(frm, text="Baza:").grid(row=2, column=0, sticky="w")
        db_var = tk.StringVar(value="nuccore")
        db = ttk.Combobox(frm, textvariable=db_var, values=["nuccore", "protein"], state="readonly", width=12)
        db.grid(row=2, column=1, sticky="w", padx=(6, 0))

        ttk.Label(frm, text="Max rekordów (dla zapytania):").grid(row=3, column=0, sticky="w", pady=(6, 0))
        retmax_var = tk.IntVar(value=20)
        retmax = ttk.Spinbox(frm, from_=1, to=500, textvariable=retmax_var, width=8)
        retmax.grid(row=3, column=1, sticky="w", padx=(6, 0), pady=(6, 0))

        status_var = tk.StringVar(value="")
        ttk.Label(frm, textvariable=status_var).grid(row=4, column=0, columnspan=3, sticky="w", pady=(10, 0))

        frm.columnconfigure(2, weight=1)

        def do_download():
            q = entry.get().strip()
            if not q:
                messagebox.showwarning("Brak danych", "Wpisz accession/UID albo zapytanie.", parent=win)
                return

            db_name = db_var.get().strip()
            base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

            # 1) ustal listę ID
            ids: list[str] = []
            parts = [x for x in re.split(r"[,\s]+", q) if x]

            if parts and all(is_accession_like(x) for x in parts):
                # traktujemy jako accession/uid listę
                ids = parts
            else:
                # ESearch -> IDs
                params = {
                    "db": db_name,
                    "term": q,
                    "retmax": int(retmax_var.get()),
                    "retmode": "xml",
                }
                status_var.set("Szukam w NCBI (ESearch)...")
                win.update_idletasks()
                xml = http_get(base + "esearch.fcgi?" + urlencode(params))
                ids = parse_esearch_ids(xml)

                if not ids:
                    messagebox.showinfo("Brak wyników", "Nie znaleziono rekordów dla podanego zapytania.", parent=win)
                    status_var.set("")
                    return

            # 2) EFetch FASTA
            params = {
                "db": db_name,
                "id": ",".join(ids),
                "rettype": "fasta",
                "retmode": "text",
            }
            status_var.set("Pobieram FASTA (EFetch)...")
            win.update_idletasks()
            fasta_text = http_get(base + "efetch.fcgi?" + urlencode(params))

            if not fasta_text.strip().startswith(">"):
                messagebox.showerror("Błąd", "NCBI nie zwróciło FASTA (sprawdź zapytanie/bazę).", parent=win)
                status_var.set("")
                return

            # 3) zapisz
            out_path = filedialog.asksaveasfilename(
                parent=win,
                title="Zapisz pobrane FASTA",
                defaultextension=".fasta",
                filetypes=[("FASTA", "*.fasta *.fa *.fna *.faa"), ("Wszystkie", "*.*")]
            )
            if not out_path:
                status_var.set("")
                return

            Path(out_path).write_text(fasta_text, encoding="utf-8")
            self.log(f"Pobrano FASTA z NCBI: {out_path}")

            # 4) wczytaj do aplikacji przez Twoją funkcję load_fasta(path)
            try:
                seqs = load_fasta(out_path)  # <- używa Twojej funkcji globalnej
                self.sequences = seqs
                self.log(f"Wczytano {len(seqs)} sekwencji")
            except Exception as e:
                messagebox.showerror("Błąd FASTA", str(e), parent=win)
                return

            # 5) przełącz na podgląd
            try:
                self.tabs.select(self.tab_preview)
            except Exception:
                pass

            win.destroy()

        btns = ttk.Frame(frm)
        btns.grid(row=5, column=0, columnspan=3, sticky="e", pady=(12, 0))
        ttk.Button(btns, text="Anuluj", command=win.destroy).pack(side="right")
        ttk.Button(btns, text="Pobierz", command=do_download).pack(side="right", padx=(0, 8))

        entry.bind("<Return>", lambda e: do_download())

    def about(self):
        """Informacje o programie"""
        messagebox.showinfo(
            "O programie",
            "Analiza motywów DNA\nProjekt zaliczeniowy"
        )

    def render_results_table(self):
        # wyczyść tabelę
        for row in self.results_table.get_children():
            self.results_table.delete(row)

        if self.analysis_result is None:
            return

        mode = self.normalization_mode.get()
        mat = self.analysis_result.matrix(mode)
        seq_ids = self.analysis_result.seq_ids
        motifs = self.analysis_result.motifs

        # dynamiczne kolumny
        columns = ["Sekwencja"] + motifs + ["SUMA"]
        self.results_table["columns"] = columns

        for col in columns:
            self.results_table.heading(col, text=f"{col} ⇅", command=lambda c=col: self.sort_column(c))
            self.results_table.column(col, width=90, anchor="center")

        for i, sid in enumerate(seq_ids):
            row_vals = mat[i, :]
            if mode == "raw":
                cells = [str(int(v)) for v in row_vals]
                total = str(int(row_vals.sum()))
            else:
                cells = [f"{v:.1f}" for v in row_vals]
                total = f"{row_vals.sum():.1f}"

            self.results_table.insert("", "end", values=[sid] + cells + [total])

# ==================================================
# URUCHOMIENIE PROGRAMU
# ==================================================

if __name__ == "__main__":
    app = DNAApp()
    app.mainloop()
