from __future__ import annotations

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
    try:
        return (_IUPAC_MASK[seq_char.upper()] & _IUPAC_MASK[motif_char.upper()]) != 0
    except KeyError:
        return False


def find_positions(sequence: str, motif: str) -> list[int]:
    seq = sequence.upper().replace("U", "T")
    mot = motif.upper().replace("U", "T")

    positions: list[int] = []
    motif_len = len(mot)
    if motif_len == 0 or motif_len > len(seq):
        return positions

    for start in range(len(seq) - motif_len + 1):
        if all(matches_iupac(seq[start + offset], mot[offset]) for offset in range(motif_len)):
            positions.append(start)
    return positions


def count_matches(sequence: str, motif: str) -> int:
    if not sequence or not motif:
        return 0

    seq = sequence.upper().replace("U", "T")
    mot = motif.upper().replace("U", "T")

    try:
        seq_masks = [_IUPAC_MASK[ch] for ch in seq]
        mot_masks = [_IUPAC_MASK[ch] for ch in mot]
    except KeyError as exc:
        raise ValueError(f"Nieznany symbol IUPAC: {exc.args[0]}") from exc

    motif_len = len(mot_masks)
    if motif_len > len(seq_masks):
        return 0

    count = 0
    for start in range(len(seq_masks) - motif_len + 1):
        if all((seq_masks[start + offset] & mot_masks[offset]) != 0 for offset in range(motif_len)):
            count += 1
    return count
