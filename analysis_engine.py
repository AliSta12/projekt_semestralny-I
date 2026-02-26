from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional

import numpy as np


@dataclass(frozen=True)
class AnalysisResult:
    seq_ids: List[str]
    motifs: List[str]
    raw_counts: np.ndarray          # shape: (n_seq, n_motifs), dtype=int
    seq_lengths: np.ndarray         # shape: (n_seq,), dtype=int

    fasta_path: Optional[str]
    created_at: str                 # ISO string
    iupac_enabled: bool = True
    overlapping: bool = True

    def matrix(self, mode: str) -> np.ndarray:
        """
        Zwraca macierz wartości dla trybu:
          - 'raw'  -> int
          - 'norm' -> float (na 1000 nt)
        """
        mode = (mode or "raw").lower()
        if mode == "raw":
            return self.raw_counts.astype(float)  # float dla spójności rysowania/formatowania
        elif mode == "norm":
            # (count / length) * 1000, uwaga na długość=0
            lengths = self.seq_lengths.astype(float)
            denom = np.where(lengths > 0, lengths, 1.0)
            return (self.raw_counts.astype(float) / denom[:, None]) * 1000.0
        else:
            raise ValueError(f"Nieznany tryb: {mode}. Dozwolone: raw, norm")


def compute_analysis(
    sequences: Dict[str, str],
    motifs: List[str],
    iupac_count_fn,
    fasta_path: Optional[str] = None,
) -> AnalysisResult:
    """
    Liczy RAW macierz wystąpień (overlapping) dla wszystkich sekwencji x motywów.
    Normalizacja jest liczona później (on demand).
    """
    if not sequences:
        raise ValueError("Brak sekwencji do analizy.")
    if not motifs:
        raise ValueError("Brak motywów do analizy.")

    seq_ids = list(sequences.keys())
    mot_list = list(motifs)

    n_seq = len(seq_ids)
    n_mot = len(mot_list)

    raw = np.zeros((n_seq, n_mot), dtype=int)
    lengths = np.zeros((n_seq,), dtype=int)

    for i, sid in enumerate(seq_ids):
        seq = sequences[sid]
        lengths[i] = len(seq)

        for j, mot in enumerate(mot_list):
            raw[i, j] = int(iupac_count_fn(seq, mot))

    return AnalysisResult(
        seq_ids=seq_ids,
        motifs=mot_list,
        raw_counts=raw,
        seq_lengths=lengths,
        fasta_path=fasta_path,
        created_at=datetime.now().isoformat(timespec="seconds"),
        iupac_enabled=True,
        overlapping=True,
    )