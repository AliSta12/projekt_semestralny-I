import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from export_manager import export_results_csv, export_report, export_figure, export_all_figures


class ExportTab(ttk.Frame):
    """
    GUI eksportu. Nie liczy niczego. Pobiera:
      - analysis_result (model)
      - figures (dict nazw -> matplotlib Figure)
    """

    def __init__(self, parent, get_result, get_figures, log_fn):
        super().__init__(parent)
        self.get_result = get_result
        self.get_figures = get_figures
        self.log = log_fn

        self._build_ui()

    def _build_ui(self):
        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True, padx=8, pady=8)

        t1 = ttk.Frame(nb)
        t2 = ttk.Frame(nb)
        t3 = ttk.Frame(nb)

        nb.add(t1, text="Wyniki (CSV/TSV)")
        nb.add(t2, text="Raport (TXT/HTML)")
        nb.add(t3, text="Wykresy")

        # ---------- Wyniki ----------
        self.sep = tk.StringVar(value=";")
        self.mode = tk.StringVar(value="raw")
        self.layout = tk.StringVar(value="wide")

        box = ttk.LabelFrame(t1, text="Eksport wyników")
        box.pack(fill="x", padx=8, pady=8)

        r = ttk.Frame(box); r.pack(fill="x", padx=8, pady=6)
        ttk.Label(r, text="Format:").pack(side="left")
        ttk.Radiobutton(r, text="CSV (;)", value=";", variable=self.sep).pack(side="left", padx=6)
        ttk.Radiobutton(r, text="CSV (,)", value=",", variable=self.sep).pack(side="left", padx=6)
        ttk.Radiobutton(r, text="TSV", value="\t", variable=self.sep).pack(side="left", padx=6)

        r = ttk.Frame(box); r.pack(fill="x", padx=8, pady=6)
        ttk.Label(r, text="Tryb:").pack(side="left")
        ttk.Radiobutton(r, text="Surowe", value="raw", variable=self.mode).pack(side="left", padx=6)
        ttk.Radiobutton(r, text="Na 1000 nt", value="norm", variable=self.mode).pack(side="left", padx=6)

        r = ttk.Frame(box); r.pack(fill="x", padx=8, pady=6)
        ttk.Label(r, text="Układ:").pack(side="left")
        ttk.Radiobutton(r, text="Wide (macierz)", value="wide", variable=self.layout).pack(side="left", padx=6)
        ttk.Radiobutton(r, text="Long (tidy)", value="long", variable=self.layout).pack(side="left", padx=6)

        ttk.Button(t1, text="Zapisz wyniki…", command=self._export_results).pack(anchor="e", padx=16, pady=12)

        # ---------- Raport ----------
        self.rep_fmt = tk.StringVar(value="txt")
        self.rep_mode = tk.StringVar(value="raw")
        self.rep_table = tk.BooleanVar(value=True)

        box = ttk.LabelFrame(t2, text="Eksport raportu")
        box.pack(fill="x", padx=8, pady=8)

        r = ttk.Frame(box); r.pack(fill="x", padx=8, pady=6)
        ttk.Label(r, text="Format:").pack(side="left")
        ttk.Radiobutton(r, text="TXT", value="txt", variable=self.rep_fmt).pack(side="left", padx=6)
        ttk.Radiobutton(r, text="HTML", value="html", variable=self.rep_fmt).pack(side="left", padx=6)

        r = ttk.Frame(box); r.pack(fill="x", padx=8, pady=6)
        ttk.Label(r, text="Tryb:").pack(side="left")
        ttk.Radiobutton(r, text="Surowe", value="raw", variable=self.rep_mode).pack(side="left", padx=6)
        ttk.Radiobutton(r, text="Na 1000 nt", value="norm", variable=self.rep_mode).pack(side="left", padx=6)

        ttk.Checkbutton(box, text="Dołącz tabelę wyników", variable=self.rep_table).pack(anchor="w", padx=12, pady=(0, 6))
        ttk.Button(t2, text="Zapisz raport…", command=self._export_report).pack(anchor="e", padx=16, pady=12)

        # ---------- Wykresy ----------
        self.plot_fmt = tk.StringVar(value="png")
        self.plot_dpi = tk.IntVar(value=200)

        box = ttk.LabelFrame(t3, text="Eksport wykresów")
        box.pack(fill="both", expand=True, padx=8, pady=8)

        r = ttk.Frame(box); r.pack(fill="x", padx=8, pady=6)
        ttk.Label(r, text="Format:").pack(side="left")
        ttk.Radiobutton(r, text="PNG", value="png", variable=self.plot_fmt).pack(side="left", padx=6)
        ttk.Radiobutton(r, text="SVG", value="svg", variable=self.plot_fmt).pack(side="left", padx=6)
        ttk.Radiobutton(r, text="PDF", value="pdf", variable=self.plot_fmt).pack(side="left", padx=6)

        ttk.Label(r, text="DPI (PNG):").pack(side="left", padx=(18, 4))
        ttk.Spinbox(r, from_=72, to=600, textvariable=self.plot_dpi, width=6).pack(side="left")

        ttk.Label(box, text="Dostępne wykresy:").pack(anchor="w", padx=10, pady=(6, 0))
        self.lb = tk.Listbox(box, height=8)
        self.lb.pack(fill="both", expand=True, padx=10, pady=(4, 8))

        b = ttk.Frame(box); b.pack(fill="x", padx=10, pady=(0, 10))
        ttk.Button(b, text="Odśwież", command=self._refresh_plots).pack(side="left")
        ttk.Button(b, text="Zapisz wybrany…", command=self._export_selected_plot).pack(side="right", padx=6)
        ttk.Button(b, text="Zapisz wszystkie…", command=self._export_all_plots).pack(side="right")

        self._refresh_plots()

    def _require_result(self):
        result = self.get_result()
        if result is None:
            messagebox.showwarning("Brak danych", "Najpierw uruchom analizę.")
            return None
        return result

    def _export_results(self):
        result = self._require_result()
        if result is None:
            return

        sep = self.sep.get()
        ext = "tsv" if sep == "\t" else "csv"
        path = filedialog.asksaveasfilename(
            title="Zapisz wyniki",
            defaultextension=f".{ext}",
            filetypes=[("CSV/TSV", f"*.{ext}"), ("All files", "*.*")]
        )
        if not path:
            return

        export_results_csv(
            result=result,
            path=path,
            sep=sep,
            mode=self.mode.get(),
            layout=self.layout.get(),
        )
        self.log(f"Eksport wyników: {path}")
        messagebox.showinfo("OK", "Zapisano wyniki.")

    def _export_report(self):
        result = self._require_result()
        if result is None:
            return

        fmt = self.rep_fmt.get()
        ext = "txt" if fmt == "txt" else "html"
        path = filedialog.asksaveasfilename(
            title="Zapisz raport",
            defaultextension=f".{ext}",
            filetypes=[(ext.upper(), f"*.{ext}"), ("All files", "*.*")]
        )
        if not path:
            return

        export_report(
            result=result,
            path=path,
            mode=self.rep_mode.get(),
            fmt=fmt,
            include_table=bool(self.rep_table.get()),
        )
        self.log(f"Eksport raportu: {path}")
        messagebox.showinfo("OK", "Zapisano raport.")

    def _refresh_plots(self):
        self.lb.delete(0, tk.END)
        figs = self.get_figures() or {}
        for name in sorted(figs.keys()):
            self.lb.insert(tk.END, name)

    def _export_selected_plot(self):
        figs = self.get_figures() or {}
        if not figs:
            messagebox.showwarning("Brak wykresów", "Brak wykresów do eksportu.")
            return

        sel = self.lb.curselection()
        if not sel:
            messagebox.showwarning("Brak wyboru", "Wybierz wykres z listy.")
            return

        name = self.lb.get(sel[0])
        fig = figs.get(name)
        if fig is None:
            messagebox.showerror("Błąd", "Nie znaleziono figury.")
            return

        fmt = self.plot_fmt.get()
        path = filedialog.asksaveasfilename(
            title="Zapisz wykres",
            defaultextension=f".{fmt}",
            filetypes=[(fmt.upper(), f"*.{fmt}"), ("All files", "*.*")]
        )
        if not path:
            return

        export_figure(fig, path, fmt=fmt, dpi=int(self.plot_dpi.get()))
        self.log(f"Eksport wykresu '{name}': {path}")
        messagebox.showinfo("OK", "Zapisano wykres.")

    def _export_all_plots(self):
        figs = self.get_figures() or {}
        if not figs:
            messagebox.showwarning("Brak wykresów", "Brak wykresów do eksportu.")
            return

        directory = filedialog.askdirectory(title="Wybierz folder")
        if not directory:
            return

        fmt = self.plot_fmt.get()
        n = export_all_figures(figs, directory, fmt=fmt, dpi=int(self.plot_dpi.get()))
        self.log(f"Eksport {n} wykresów do: {directory}")
        messagebox.showinfo("OK", f"Zapisano {n} wykresów.")