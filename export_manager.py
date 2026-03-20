from __future__ import annotations

from dataclasses import asdict
from typing import Dict, Optional, Literal
import os

from analysis_engine import AnalysisResult


Layout = Literal["long", "wide"]
ReportFmt = Literal["txt", "html"]


def export_results_csv(
    result: AnalysisResult,
    path: str,
    sep: str = ";",
    mode: str = "raw",
    layout: Layout = "wide",
    include_sum: bool = True,
) -> None:
    """
    Eksport wyników do CSV/TSV.
    - layout='wide': wiersz=sekwencja, kolumny=motywy (+ SUMA)
    - layout='long': (sequence_id, motif, value)
    """
    mat = result.matrix(mode)
    seq_ids = result.seq_ids
    motifs = result.motifs

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    with open(path, "w", encoding="utf-8", newline="") as f:
        if layout == "wide":
            header = ["Sekwencja"] + motifs + (["SUMA"] if include_sum else [])
            f.write(sep.join(header) + "\n")

            for i, sid in enumerate(seq_ids):
                row_vals = mat[i, :]

                # format: raw -> int, norm -> 1 miejsce
                if mode == "raw":
                    cells = [str(int(v)) for v in row_vals]
                    total = str(int(row_vals.sum()))
                else:
                    cells = [f"{v:.1f}" for v in row_vals]
                    total = f"{row_vals.sum():.1f}"

                out = [sid] + cells + ([total] if include_sum else [])
                f.write(sep.join(out) + "\n")

        elif layout == "long":
            f.write(sep.join(["sequence_id", "motif", "value"]) + "\n")
            for i, sid in enumerate(seq_ids):
                for j, mot in enumerate(motifs):
                    v = mat[i, j]
                    if mode == "raw":
                        val = str(int(v))
                    else:
                        val = f"{v:.1f}"
                    f.write(sep.join([sid, mot, val]) + "\n")
        else:
            raise ValueError("layout must be 'wide' or 'long'")


def export_report(
    result: AnalysisResult,
    path: str,
    mode: str = "raw",
    fmt: ReportFmt = "txt",
    include_table: bool = True,
    table_max_rows: int = 200,
    figures: Optional[Dict[str, object]] = None,
    selected_figure_names: Optional[list[str]] = None,
) -> None:
    """
    Prosty raport TXT/HTML: metadane + (opcjonalnie) tabela wyników (wide).
    """
    mat = result.matrix(mode)
    motifs = result.motifs
    seq_ids = result.seq_ids

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    def h(t: str) -> str:
        return f"<h2>{t}</h2>\n" if fmt == "html" else f"\n== {t} ==\n"

    def p(t: str) -> str:
        return f"<p>{t}</p>\n" if fmt == "html" else t + "\n"

    out = ""
    if fmt == "html":
        out += "<!doctype html><html><meta charset='utf-8'><body>\n"

    out += h("Raport analizy motywów DNA")
    out += p(f"Data utworzenia wyniku: {result.created_at}")
    if result.fasta_path:
        out += p(f"Plik FASTA: {result.fasta_path}")
    out += p(f"Tryb: {mode} ({'surowe' if mode=='raw' else 'na 1000 nt'})")
    out += p(f"IUPAC: {'enabled' if result.iupac_enabled else 'disabled'}")
    out += p(f"Dopasowania: {'overlapping' if result.overlapping else 'non-overlapping'}")
    out += p(f"Liczba sekwencji: {len(seq_ids)}")
    out += p(f"Liczba motywów: {len(motifs)}")

    # proste statystyki
    out += h("Podsumowanie")
    total_sum = mat.sum()
    out += p(f"Suma wszystkich wartości w macierzy: {total_sum:.1f}" if mode != "raw" else f"Suma wszystkich zliczeń: {int(total_sum)}")

    if include_table:
        out += h("Tabela wyników (wide)")
        # budujemy TSV (łatwe do wklejenia do Excela)
        lines = []
        header = ["Sekwencja"] + motifs + ["SUMA"]
        lines.append("\t".join(header))

        n = min(len(seq_ids), table_max_rows)
        for i in range(n):
            sid = seq_ids[i]
            row = mat[i, :]
            if mode == "raw":
                cells = [str(int(v)) for v in row]
                total = str(int(row.sum()))
            else:
                cells = [f"{v:.1f}" for v in row]
                total = f"{row.sum():.1f}"
            lines.append("\t".join([sid] + cells + [total]))

        if len(seq_ids) > table_max_rows:
            lines.append(f"... (ucięto do {table_max_rows} wierszy)")

        table_txt = "\n".join(lines)

        if fmt == "html":
            out += "<pre>\n" + table_txt + "\n</pre>\n"
        else:
            out += table_txt + "\n"

    # --- WYKRESY (tylko dla HTML) ---
    if fmt == "html" and figures:
        out += h("Wykresy")

        names = selected_figure_names or []

        import base64
        from io import BytesIO

        for name in names:
            fig = figures.get(name)
            if fig is None:
                continue

            buf = BytesIO()
            fig.savefig(buf, format="png", bbox_inches="tight", dpi=150)
            buf.seek(0)

            img_base64 = base64.b64encode(buf.read()).decode("utf-8")

            out += f"<h3>{name}</h3>\n"
            out += f"<img src='data:image/png;base64,{img_base64}' style='max-width:100%;'><br><br>\n"

    if fmt == "html":
        out += "</body></html>"

    with open(path, "w", encoding="utf-8") as f:
        f.write(out)


def export_figure(fig, path: str, fmt: str = "png", dpi: int = 200) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    kwargs = {"bbox_inches": "tight"}
    if fmt.lower() == "png":
        kwargs["dpi"] = int(dpi)
    fig.savefig(path, **kwargs)


def export_all_figures(figures: Dict[str, object], directory: str, fmt: str = "png", dpi: int = 200) -> int:
    os.makedirs(directory, exist_ok=True)
    n = 0
    for name, fig in figures.items():
        safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in name)
        path = os.path.join(directory, f"{safe}.{fmt}")
        export_figure(fig, path, fmt=fmt, dpi=dpi)
        n += 1
    return n
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import textwrap


def export_report_pdf(
    result: AnalysisResult,
    path: str,
    figures: Dict[str, object],
    selected_figure_names: list[str],
    mode: str = "raw",
) -> None:
    """
    Eksport jednego PDF zawierającego:
    1) stronę z raportem tekstowym
    2) wybrane wykresy (po jednej stronie na wykres)
    """
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    mat = result.matrix(mode)
    motifs = result.motifs
    seq_ids = result.seq_ids

    report_lines = [
        "Raport analizy motywów DNA",
        "",
        f"Data utworzenia wyniku: {result.created_at}",
        f"Plik FASTA: {result.fasta_path or '-'}",
        f"Tryb: {mode} ({'surowe' if mode == 'raw' else 'na 1000 nt'})",
        f"IUPAC: {'enabled' if result.iupac_enabled else 'disabled'}",
        f"Dopasowania: {'overlapping' if result.overlapping else 'non-overlapping'}",
        f"Liczba sekwencji: {len(seq_ids)}",
        f"Liczba motywów: {len(motifs)}",
        "",
        "Podsumowanie",
        f"Suma wszystkich wartości w macierzy: {mat.sum():.1f}" if mode != "raw"
        else f"Suma wszystkich zliczeń: {int(mat.sum())}",
        "",
        "Motywy:",
        ", ".join(motifs) if motifs else "-",
        "",
        "Sekwencje:",
        ", ".join(seq_ids[:20]) + (" ..." if len(seq_ids) > 20 else ""),
    ]

    with PdfPages(path) as pdf:
        # --- strona 1: raport tekstowy ---
        fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait in inches
        fig.patch.set_facecolor("white")

        y = 0.97
        line_height = 0.028

        for raw_line in report_lines:
            wrapped = textwrap.wrap(str(raw_line), width=95) or [""]
            for line in wrapped:
                fig.text(0.06, y, line, fontsize=10, va="top", ha="left")
                y -= line_height
                if y < 0.05:
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig)
                    fig = plt.figure(figsize=(8.27, 11.69))
                    fig.patch.set_facecolor("white")
                    y = 0.97

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # --- kolejne strony: wybrane wykresy ---
        for name in selected_figure_names:
            fig = figures.get(name)
            if fig is not None:
                pdf.savefig(fig, bbox_inches="tight")