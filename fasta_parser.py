from __future__ import annotations

from typing import Dict

_ALLOWED_IUPAC = set("ACGTRYSWKMBDHVN")


def load_fasta(path: str, header_id_max_len: int = 8) -> Dict[str, str]:
    """
    Wczytuje plik FASTA i zwraca słownik {ID: sekwencja}.

    Zasady czyszczenia:
    - pełne DNA IUPAC: A C G T R Y S W K M B D H V N
    - U -> T
    - usuwa alignment gapy '-' i '.'
    - '?' -> 'N'
    - usuwa spacje, taby i cyfry
    - inne znaki powodują ValueError
    """
    sequences: Dict[str, str] = {}
    current_id: str | None = None

    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line_num, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                current_id = line[1:].strip()[:header_id_max_len]
                if not current_id:
                    raise ValueError(f"Pusty nagłówek w linii {line_num}")
                if current_id in sequences:
                    raise ValueError(f"Duplikat ID: {current_id}")
                sequences[current_id] = ""
                continue

            if current_id is None:
                raise ValueError("Plik nie zaczyna się od nagłówka FASTA")

            seq_line = _normalize_sequence_line(line)
            bad = {char for char in seq_line if char not in _ALLOWED_IUPAC}
            if bad:
                raise ValueError(f"Niepoprawne znaki {sorted(bad)} w linii {line_num}")

            sequences[current_id] += seq_line

    if not sequences:
        raise ValueError("Nie znaleziono sekwencji")

    return sequences


def _normalize_sequence_line(line: str) -> str:
    seq_line = line.upper()
    seq_line = seq_line.replace("U", "T")
    seq_line = seq_line.replace("-", "").replace(".", "")
    seq_line = seq_line.replace("?", "N")
    seq_line = seq_line.replace(" ", "").replace("\t", "")
    seq_line = "".join(char for char in seq_line if not char.isdigit())
    return seq_line
