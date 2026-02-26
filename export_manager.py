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