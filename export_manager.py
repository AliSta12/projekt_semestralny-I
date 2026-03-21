from __future__ import annotations

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
    mat = result.matrix(mode)
    seq_ids = result.seq_ids
    motifs = result.motifs

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    with open(path, "w", encoding="utf-8", newline="") as handle:
        if layout == "wide":
            header = ["Sekwencja"] + motifs + (["SUMA"] if include_sum else [])
            handle.write(sep.join(header) + "\n")

            for idx, sequence_id in enumerate(seq_ids):
                row_vals = mat[idx, :]
                if mode == "raw":
                    cells = [str(int(value)) for value in row_vals]
                    total = str(int(row_vals.sum()))
                else:
                    cells = [f"{value:.1f}" for value in row_vals]
                    total = f"{row_vals.sum():.1f}"

                output_row = [sequence_id] + cells + ([total] if include_sum else [])
                handle.write(sep.join(output_row) + "\n")
            return

        if layout == "long":
            handle.write(sep.join(["sequence_id", "motif", "value"]) + "\n")
            for seq_idx, sequence_id in enumerate(seq_ids):
                for motif_idx, motif in enumerate(motifs):
                    value = mat[seq_idx, motif_idx]
                    formatted = str(int(value)) if mode == "raw" else f"{value:.1f}"
                    handle.write(sep.join([sequence_id, motif, formatted]) + "\n")
            return

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
    mat = result.matrix(mode)
    motifs = result.motifs
    seq_ids = result.seq_ids

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    def h(text: str) -> str:
        return f"<h2>{text}</h2>\n" if fmt == "html" else f"\n== {text} ==\n"

    def p(text: str) -> str:
        return f"<p>{text}</p>\n" if fmt == "html" else text + "\n"

    out = ""
    if fmt == "html":
        out += "<!doctype html><html><meta charset='utf-8'><body>\n"

    out += h("Raport analizy motywów DNA")
    out += p(f"Data utworzenia wyniku: {result.created_at}")
    if result.fasta_path:
        out += p(f"Plik FASTA: {result.fasta_path}")
    out += p(f"Tryb: {mode} ({'surowe' if mode == 'raw' else 'na 1000 nt'})")
    out += p(f"IUPAC: {'enabled' if result.iupac_enabled else 'disabled'}")
    out += p(f"Dopasowania: {'overlapping' if result.overlapping else 'non-overlapping'}")
    out += p(f"Liczba sekwencji: {len(seq_ids)}")
    out += p(f"Liczba motywów: {len(motifs)}")

    out += h("Podsumowanie")
    total_sum = mat.sum()
    out += p(
        f"Suma wszystkich wartości w macierzy: {total_sum:.1f}"
        if mode != "raw"
        else f"Suma wszystkich zliczeń: {int(total_sum)}"
    )

    if include_table:
        out += h("Tabela wyników (wide)")
        lines = ["\t".join(["Sekwencja"] + motifs + ["SUMA"])]
        limit = min(len(seq_ids), table_max_rows)
        for idx in range(limit):
            sequence_id = seq_ids[idx]
            row = mat[idx, :]
            if mode == "raw":
                cells = [str(int(value)) for value in row]
                total = str(int(row.sum()))
            else:
                cells = [f"{value:.1f}" for value in row]
                total = f"{row.sum():.1f}"
            lines.append("\t".join([sequence_id] + cells + [total]))
        if len(seq_ids) > table_max_rows:
            lines.append(f"... (ucięto do {table_max_rows} wierszy)")

        table_txt = "\n".join(lines)
        out += f"<pre>\n{table_txt}\n</pre>\n" if fmt == "html" else table_txt + "\n"

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

    with open(path, "w", encoding="utf-8") as handle:
        handle.write(out)



def export_figure(fig, path: str, fmt: str = "png", dpi: int = 200) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    kwargs = {"bbox_inches": "tight"}
    if fmt.lower() == "png":
        kwargs["dpi"] = int(dpi)
    fig.savefig(path, **kwargs)



def export_all_figures(figures: Dict[str, object], directory: str, fmt: str = "png", dpi: int = 200) -> int:
    os.makedirs(directory, exist_ok=True)
    count = 0
    for name, fig in figures.items():
        safe_name = "".join(char if char.isalnum() or char in "._-" else "_" for char in name)
        export_figure(fig, os.path.join(directory, f"{safe_name}.{fmt}"), fmt=fmt, dpi=dpi)
        count += 1
    return count


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
        f"Suma wszystkich wartości w macierzy: {mat.sum():.1f}" if mode != "raw" else f"Suma wszystkich zliczeń: {int(mat.sum())}",
        "",
        "Motywy:",
        ", ".join(motifs) if motifs else "-",
        "",
        "Sekwencje:",
        ", ".join(seq_ids[:20]) + (" ..." if len(seq_ids) > 20 else ""),
    ]

    with PdfPages(path) as pdf:
        fig = plt.figure(figsize=(8.27, 11.69))
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

        for name in selected_figure_names:
            fig = figures.get(name)
            if fig is not None:
                pdf.savefig(fig, bbox_inches="tight")
