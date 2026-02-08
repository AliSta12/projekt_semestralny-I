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

        # budowa interfejsu
        self.create_menu()
        self.create_layout()

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

        # --- tabela wyników ---
        results_frame = tk.Frame(self.tab_results)
        results_frame.pack(fill="both", expand=True)

        self.results_table = ttk.Treeview(
            results_frame,
            columns=("seq", "motif", "count"),
            show="headings"
        )
        self.results_table.heading("seq", text="Sekwencja")
        self.results_table.heading("motif", text="Motyw")
        self.results_table.heading("count", text="Liczba")
        self.results_table.pack(fill="both", expand=True)

        # --- tabela podsumowań ---
        self.summary_table = ttk.Treeview(results_frame, show="headings")
        self.summary_table.pack(fill="x")

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
        Otwiera okno do zarządzania motywami:
        - dodawanie przez Enter
        - usuwanie wielu zaznaczonych
        - wybór predefiniowanych motywów
        """

        win = tk.Toplevel(self)
        win.title("Dodaj/usuń motywy")
        win.geometry("350x550")
        win.transient(self)
        win.grab_set()

        # ===== LEGENDA =====
        tk.Label(win, text="Dozwolone symbole:").pack(pady=(6, 0))
        tk.Label(win, text="A C G T").pack()

        # IUPAC tooltip
        iupac_label = tk.Label(win, text="R Y S W K M B D H V N (IUPAC)", fg="blue")
        iupac_label.pack(pady=(0, 8))

        tooltip_text = (
            "R = A lub G\n"
            "Y = C lub T\n"
            "S = G lub C\n"
            "W = A lub T\n"
            "K = G lub T\n"
            "M = A lub C\n"
            "B = C lub G lub T\n"
            "D = A lub G lub T\n"
            "H = A lub C lub T\n"
            "V = A lub C lub G\n"
            "N = dowolna baza"
        )
        tooltip = None

        def on_enter_label(e):
            nonlocal tooltip
            if tooltip:
                return
            x = e.x_root + 10
            y = e.y_root + 10
            tooltip = tk.Toplevel(win)
            tooltip.wm_overrideredirect(True)
            tooltip.geometry(f"+{x}+{y}")
            tk.Label(
                tooltip,
                text=tooltip_text,
                justify="left",
                bg="#ffffe0",
                relief="solid",
                borderwidth=1,
                padx=4,
                pady=2
            ).pack()

        def on_leave_label(e):
            nonlocal tooltip
            if tooltip:
                tooltip.destroy()
                tooltip = None

        iupac_label.bind("<Enter>", on_enter_label)
        iupac_label.bind("<Leave>", on_leave_label)

        # ===== Pole wpisywania własnych =====
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
            motif_listbox.insert("end", m)
            motif_entry.delete(0, "end")
            self.log(f"Dodano motyw: {m}")

        motif_entry.bind("<Return>", add_manual)

        # ===== Lista aktualnych motywów =====
        tk.Label(win, text="Aktualne motywy:").pack(pady=(6, 0))
        motif_listbox = tk.Listbox(win, selectmode="extended")
        motif_listbox.pack(fill="both", expand=True, padx=10)

        # === AKORDEON PREDEFINIOWANYCH ===
        tk.Label(win, text="Wybierz z listy motywów:").pack(anchor="w", padx=10, pady=(8, 4))

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
        checkbox_vars = {}

        def toggle_section(sec):
            f = accordion_frames[sec]
            if f.winfo_ismapped():
                f.pack_forget()
            else:
                f.pack(fill="x", padx=20, pady=(2, 6))

        # — Eukariota —
        btn_euk = tk.Button(win, text="Eukariota", anchor="w", relief="flat", fg="blue")
        btn_euk.pack(fill="x", padx=10)

        frm_euk = tk.Frame(win)
        accordion_frames["Eukariota"] = frm_euk
        checkbox_vars["Eukariota"] = []

        for seq, desc in EU_MOTIFS:
            var = tk.BooleanVar(master=win, value=False)
            cb = tk.Checkbutton(frm_euk, text=f"{seq} — {desc}", variable=var)
            cb.pack(anchor="w")
            checkbox_vars["Eukariota"].append((var, seq))

        btn_euk.config(command=lambda: toggle_section("Eukariota"))

        # — Prokariota —
        btn_pro = tk.Button(win, text="Prokariota", anchor="w", relief="flat", fg="blue")
        btn_pro.pack(fill="x", padx=10, pady=(4, 0))

        frm_pro = tk.Frame(win)
        accordion_frames["Prokariota"] = frm_pro
        checkbox_vars["Prokariota"] = []

        for seq, desc in PRO_MOTIFS:
            var = tk.BooleanVar(master=win, value=False)
            cb = tk.Checkbutton(frm_pro, text=f"{seq} — {desc}", variable=var)
            cb.pack(anchor="w")
            checkbox_vars["Prokariota"].append((var, seq))

        btn_pro.config(command=lambda: toggle_section("Prokariota"))

        # — Dodaj zaznaczone motywy —
        def add_selected():
            for section in checkbox_vars:
                for var, seq in checkbox_vars[section]:
                    if var.get() and seq not in self.motifs:
                        self.motifs.append(seq)
                        self.motif_listbox.insert("end", seq)
                        self.log(f"Dodano predefiniowany motyw: {seq}")

        tk.Button(win, text="Dodaj zaznaczone motywy", command=add_selected).pack(pady=(8, 10))

        # — Usuwanie —
        def remove_motifs_local():
            sel = motif_listbox.curselection()
            if not sel:
                return
            for i in sorted(sel, reverse=True):
                motif = motif_listbox.get(i)
                motif_listbox.delete(i)
                self.motif_listbox.delete(i)
                self.motifs.remove(motif)
                self.log(f"Usunięto motyw: {motif}")

        tk.Button(win, text="Usuń zaznaczone", command=remove_motifs_local).pack(pady=8)

        motif_entry.focus_set()

    def run_analysis(self):
        """Uruchamia analizę motywów dla wszystkich sekwencji"""

        if not self.sequences or not self.motifs:
            messagebox.showwarning("Błąd", "Wczytaj FASTA i dodaj motywy")
            return

        for row in self.results_table.get_children():
            self.results_table.delete(row)

        for seq_id, seq in self.sequences.items():
            for motif in self.motifs:
                c = iupac_count(seq, motif)
                self.results_table.insert("", "end", values=(seq_id, motif, c))

        self.log("Analiza zakończona")
        self.draw_visualization()
        self.build_summary()

    # ==================================================
    # WIZUALIZACJA I PODSUMOWANIA
    # ==================================================

    def draw_visualization(self):
        """Rysuje mapę wystąpień motywów na sekwencjach DNA"""

        for w in self.viz_inner.winfo_children():
            w.destroy()

        if not self.sequences or not self.motifs:
            return

        left = 120
        width = 700
        line_h = 40

        palette = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628"]
        motif_colors = {m: palette[i % len(palette)] for i, m in enumerate(self.motifs)}

        legend = tk.Frame(self.viz_inner)
        legend.pack(anchor="w", pady=5)

        for m, c in motif_colors.items():
            tk.Label(legend, text=m, fg=c).pack(side="left", padx=10)

        for row, (seq_id, seq) in enumerate(self.sequences.items()):
            canvas = tk.Canvas(
                self.viz_inner,
                width=left + width + 20,
                height=line_h,
                bg="white",
                highlightthickness=0
            )
            canvas.pack(anchor="w")

            canvas.create_text(10, line_h // 2, anchor="w", text=seq_id)
            canvas.create_line(left, line_h // 2, left + width, line_h // 2, width=2)

            L = len(seq)
            if L == 0:
                continue

            for motif in self.motifs:
                color = motif_colors[motif]
                positions = iupac_find_positions(seq, motif)

                for p in positions:
                    x = left + (p / L) * width
                    line_id = canvas.create_line(
                        x,
                        line_h // 2 - 8,
                        x,
                        line_h // 2 + 8,
                        fill=color,
                        width=3
                    )

                    canvas.tag_bind(
                        line_id,
                        "<Button-1>",
                        lambda e, s=seq_id, m=motif, pos=p:
                        messagebox.showinfo(
                            "Motyw",
                            f"Sekwencja: {s}\nMotyw: {m}\nPozycja: {pos}"
                        )
                    )

    def build_summary(self):
        """Buduje tabelę podsumowań liczby motywów"""

        for c in self.summary_table.get_children():
            self.summary_table.delete(c)

        cols = ["Sekwencja"] + self.motifs + ["SUMA"]
        self.summary_table["columns"] = cols

        for c in cols:
            self.summary_table.heading(c, text=c)
            self.summary_table.column(c, width=80)

        for seq_id in self.sequences:
            row = [seq_id]
            total = 0

            for motif in self.motifs:
                s = 0
                for item in self.results_table.get_children():
                    v = self.results_table.item(item)["values"]
                    if v[0] == seq_id and v[1] == motif:
                        s += v[2]

                row.append(s)
                total += s

            row.append(total)
            self.summary_table.insert("", "end", values=row)

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
