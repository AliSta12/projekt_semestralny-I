import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog
import csv
import os




# ==================================================
# LOGIKA (funkcje)
# ==================================================

def load_fasta(path):
    sequences = {}
    current_id = None

    with open(path, "r") as f:
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
            else:
                if current_id is None:
                    raise ValueError("Plik nie zaczyna się od nagłówka FASTA")

                if any(c not in "ACGTacgt" for c in line):
                    raise ValueError(f"Niepoprawne znaki w linii {line_num}")

                sequences[current_id] += line.upper()

    if not sequences:
        raise ValueError("Nie znaleziono sekwencji")

    return sequences
#funkcja licząca motywy
def count_motif(sequence, motif):
    count = 0
    start = 0

    while True:
        pos = sequence.find(motif, start)
        if pos == -1:
            break

        count += 1
        start = pos + 1  # pozwala na overlapping

    return count
#szukanie pozycji motywów
def find_motif_positions(sequence, motif):
    positions = []
    start = 0

    while True:
        pos = sequence.find(motif, start)
        if pos == -1:
            break

        positions.append(pos)
        start = pos + 1

    return positions




# ==================================================
# GUI
# ==================================================

class DNAApp(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("Projekt 1: Analiza motywów sekwencyjnych w DNA")
        self.geometry("1000x600")

        self.sequences = {}
        self.current_fasta = None
        self.motifs = []

        self.create_menu()
        self.create_layout()

    # ================= MENU =================

    def create_menu(self):
        menubar = tk.Menu(self)

        plik = tk.Menu(menubar, tearoff=0)
        plik.add_command(label="Wczytaj plik", command=self.gui_load_file)
        plik.add_separator()
        plik.add_command(label="Wyjście", command=self.quit)

        motywy = tk.Menu(menubar, tearoff=0)
        motywy.add_command(label="Dodaj motyw", command=self.add_motif)

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

    # ================= LAYOUT =================

    def create_layout(self):
        main_frame = tk.Frame(self)
        main_frame.pack(fill="both", expand=True)

        side = tk.Frame(main_frame, width=200)
        side.pack(side="left", fill="y")

        tk.Button(side, text="Wczytaj plik", command=self.gui_load_file).pack(fill="x", pady=2)
        tk.Button(side, text="Pobierz z NCBI", command=self.download_ncbi).pack(fill="x", pady=2)
        tk.Button(side, text="Dodaj motyw", command=self.add_motif).pack(fill="x", pady=2)
        tk.Button(side, text="Usuń zaznaczone motywy", command=self.remove_motif).pack(fill="x", pady=2)
        tk.Button(side, text="Uruchom analizę", command=self.run_analysis).pack(fill="x", pady=2)
        tk.Button(side, text="Eksportuj CSV/PDF", command=self.export_data).pack(fill="x", pady=2)

        right = tk.Frame(main_frame)
        right.pack(side="right", fill="both", expand=True)

        self.tabs = ttk.Notebook(right)
        self.tabs.pack(fill="both", expand=True)

        self.tab_preview = ttk.Frame(self.tabs)
        self.tab_motifs = ttk.Frame(self.tabs)
        self.tab_results = ttk.Frame(self.tabs)
        self.tab_viz = ttk.Frame(self.tabs)
        self.tab_export = ttk.Frame(self.tabs)

        self.tabs.add(self.tab_preview, text="Podgląd sekwencji")
        self.tabs.add(self.tab_motifs, text="Wybór motywów")
        self.tabs.add(self.tab_results, text="Wyniki analizy")
        self.tabs.add(self.tab_viz, text="Wizualizacja")
        self.tabs.add(self.tab_export, text="Eksport")

        # Podgląd FASTA
        self.preview_box = tk.Text(self.tab_preview)
        self.preview_box.pack(fill="both", expand=True)

        # === Wizualizacja (Canvas + Scrollbar) ===
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

        def _on_configure(event):
            self.viz_canvas.configure(scrollregion=self.viz_canvas.bbox("all"))

        self.viz_inner.bind("<Configure>", _on_configure)

        # === Widok plików CSV ===
        tk.Label(self.tab_export, text="Pliki wynikowe CSV:").pack(anchor="w")

        self.export_listbox = tk.Listbox(self.tab_export)
        self.export_listbox.pack(fill="both", expand=True)

        tk.Button(
            self.tab_export,
            text="Odśwież listę",
            command=self.refresh_csv_files
        ).pack(pady=5)

        self.motif_listbox = tk.Listbox(self.tab_motifs)
        self.motif_listbox.pack(fill="both", expand=True)

        # TABELA WYNIKÓW – TYLKO RAZ
        self.results_table = ttk.Treeview(
            self.tab_results,
            columns=("seq", "motif", "count"),
            show="headings"
        )

        self.results_table.heading("seq", text="Sekwencja")
        self.results_table.heading("motif", text="Motyw")
        self.results_table.heading("count", text="Liczba")

        # kontener na wyniki
        results_frame = tk.Frame(self.tab_results)
        results_frame.pack(fill="both", expand=True)

        # tabela szczegółowa (góra)
        self.results_table = ttk.Treeview(
            results_frame,
            columns=("seq", "motif", "count"),
            show="headings"
        )
        self.results_table.heading("seq", text="Sekwencja")
        self.results_table.heading("motif", text="Motyw")
        self.results_table.heading("count", text="Liczba")
        self.results_table.pack(fill="both", expand=True)

        # tabela agregacji (dół)
        self.summary_table = ttk.Treeview(results_frame, show="headings")
        self.summary_table.pack(fill="x")

        # LOGI
        log_frame = tk.Frame(self)
        log_frame.pack(fill="x")

        tk.Label(log_frame, text="Logi / komunikaty:").pack(anchor="w")

        self.log_box = tk.Text(log_frame, height=6)
        self.log_box.pack(fill="x")

    def log(self, msg):
        self.log_box.insert("end", msg + "\n")
        self.log_box.see("end")

    def refresh_csv_files(self):
        self.export_listbox.delete(0, "end")

        project_dir = os.path.dirname(os.path.abspath(__file__))

        for f in os.listdir(project_dir):
            if f.lower().endswith(".csv"):
                self.export_listbox.insert("end", f)

        if self.export_listbox.size() == 0:
            self.export_listbox.insert("end", "(brak plików CSV)")

    # ================= GUI CALLBACKS =================

    def gui_load_file(self):
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

            total_len = sum(len(v) for v in self.sequences.values())

            for i, (k, v) in enumerate(self.sequences.items()):
                if i == 20:
                    self.preview_box.insert("end", "\n... (reszta ukryta)\n")
                    break

                self.preview_box.insert("end", f">{k}\n")
                self.preview_box.insert(
                    "end",
                    v[:200] + ("..." if len(v) > 200 else "") + "\n\n"
                )

            self.log(f"Wczytano plik: {path}")
            self.log(f"Sekwencje: {len(self.sequences)}")
            self.log(f"Łączna długość: {total_len}")

            self.tabs.select(self.tab_preview)

        except Exception as e:
            messagebox.showerror("Błąd FASTA", str(e))

    def add_motif(self):
        text = simpledialog.askstring(
            "Dodaj motywy",
            "Podaj motyw lub kilka motywów.\n"
            "Jeśli kilka – oddziel przecinkami.\n\n"
            "Przykład:\nATG, CGT, AAA"
        )

        if not text:
            return

        motifs = [m.strip().upper() for m in text.split(",") if m.strip()]

        for motif in motifs:
            if any(c not in "ACGT" for c in motif):
                messagebox.showerror(
                    "Błąd",
                    f"Niepoprawny motyw: {motif}\nDozwolone tylko A C G T"
                )
                continue

            if motif not in self.motifs:
                self.motifs.append(motif)
                self.motif_listbox.insert("end", motif)
                self.log(f"Dodano motyw: {motif}")

        self.tabs.select(self.tab_motifs)

    def remove_motif(self):
        selected = self.motif_listbox.curselection()

        if not selected:
            messagebox.showinfo(
                "Usuń motyw",
                "Zaznacz motyw lub kilka motywów na liście."
            )
            return

        for i in reversed(selected):
            motif = self.motif_listbox.get(i)
            self.motif_listbox.delete(i)
            self.motifs.remove(motif)
            self.log(f"Usunięto motyw: {motif}")

    def run_analysis(self):
        if not self.sequences:
            messagebox.showwarning("Brak danych", "Najpierw wczytaj FASTA")
            return

        if not self.motifs:
            messagebox.showwarning("Brak motywów", "Dodaj przynajmniej jeden motyw")
            return

        for row in self.results_table.get_children():
            self.results_table.delete(row)

        for seq_id, seq in self.sequences.items():
            for motif in self.motifs:
                c = count_motif(seq, motif)

                self.results_table.insert(
                    "",
                    "end",
                    values=(seq_id, motif, c)
                )

        self.log("Analiza zakończona")
        self.tabs.select(self.tab_results)
        self.draw_visualization()
        self.build_summary()

    def draw_visualization(self):
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
                positions = find_motif_positions(seq, motif)

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

                    def make_handler(seq_id=seq_id, motif=motif, pos=p):
                        return lambda e: messagebox.showinfo(
                            "Motyw",
                            f"Sekwencja: {seq_id}\nMotyw: {motif}\nPozycja: {pos}"
                        )

                    canvas.tag_bind(line_id, "<Button-1>", make_handler())

    def build_summary(self):
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

    def download_ncbi(self):
        self.log("NCBI (placeholder)")

    def export_data(self):
        if not self.results_table.get_children():
            messagebox.showinfo("Eksport", "Brak wyników do eksportu.")
            return

        path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv")],
            title="Zapisz wyniki jako CSV"
        )

        if not path:
            return

        try:
            with open(path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["Sekwencja", "Motyw", "Liczba"])

                for row in self.results_table.get_children():
                    writer.writerow(self.results_table.item(row)["values"])

            self.log(f"Wyeksportowano wyniki do: {path}")

        except Exception as e:
            messagebox.showerror("Błąd eksportu", str(e))

    def about(self):
        messagebox.showinfo("O programie", "Analiza motywów DNA\nSzkielet projektu")


# ==================================================

if __name__ == "__main__":
    app = DNAApp()
    app.mainloop()
