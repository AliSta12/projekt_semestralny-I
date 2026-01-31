import tkinter as tk
from tkinter import ttk, filedialog, messagebox


# ==================================================
# LOGIKA (funkcje)
# ==================================================

def load_fasta(path):
    sequences = {}
    current_id = None

    with open(path) as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                current_id = line[1:]
                sequences[current_id] = ""
            else:
                sequences[current_id] += line

    return sequences


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

        # Placeholder reszta
        for tab in [
            self.tab_motifs,
            self.tab_results,
            self.tab_viz,
            self.tab_export,
        ]:
            tk.Label(tab, text="(puste – do implementacji)").pack(pady=20)

        log_frame = tk.Frame(self)
        log_frame.pack(fill="x")

        tk.Label(log_frame, text="Logi / komunikaty:").pack(anchor="w")

        self.log_box = tk.Text(log_frame, height=6)
        self.log_box.pack(fill="x")

    # ================= LOG =================

    def log(self, msg):
        self.log_box.insert("end", msg + "\n")
        self.log_box.see("end")

    # ================= GUI CALLBACKS =================

    def gui_load_file(self):
        path = filedialog.askopenfilename(
            title="Wybierz FASTA",
            filetypes=[("FASTA", "*.fasta *.fa *.fna"), ("All", "*.*")]
        )

        if not path:
            return

        try:
            self.sequences = load_fasta(path)
            self.current_fasta = path

            self.preview_box.delete("1.0", "end")

            for k, v in self.sequences.items():
                self.preview_box.insert("end", f">{k}\n{v}\n\n")

            self.log(f"Plik wczytany: {path}")
            self.log(f"Liczba sekwencji: {len(self.sequences)}")

            self.tabs.select(self.tab_preview)

        except Exception as e:
            messagebox.showerror("Błąd", str(e))

    def add_motif(self):
        self.log("Dodaj motyw (placeholder)")

    def run_analysis(self):
        self.log("Analiza (placeholder)")

    def download_ncbi(self):
        self.log("NCBI (placeholder)")

    def export_data(self):
        self.log("Eksport (placeholder)")

    def about(self):
        messagebox.showinfo("O programie", "Analiza motywów DNA\nSzkielet projektu")


# ==================================================

if __name__ == "__main__":
    app = DNAApp()
    app.mainloop()
