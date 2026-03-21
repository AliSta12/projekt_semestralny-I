from __future__ import annotations

import re
from urllib.parse import urlencode

import requests

NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_USER_AGENT = "DNAApp/1.0 (contact: example@example.com)"


_ACCESSION_PATTERN = re.compile(r"^[A-Z]{1,3}_[0-9]+\.[0-9]+$|^[A-Z]{1,2}[0-9]+\.[0-9]+$")


def is_accession_like(text: str) -> bool:
    text = text.strip()
    if not text:
        return False
    if text.isdigit():
        return True
    return bool(_ACCESSION_PATTERN.match(text))



def http_get(url: str, timeout: int = 30) -> str:
    headers = {"User-Agent": NCBI_USER_AGENT}
    response = requests.get(url, headers=headers, timeout=timeout)
    response.raise_for_status()
    return response.text



def fetch_fasta_by_ids(ids: list[str], db_name: str = "nuccore") -> str:
    if not ids:
        raise ValueError("Brak accession/UID do pobrania.")

    invalid = [item for item in ids if not is_accession_like(item)]
    if invalid:
        raise ValueError("Podaj accession number lub UID, bez wyszukiwania tekstowego.")

    params = {
        "db": db_name,
        "id": ",".join(ids),
        "rettype": "fasta",
        "retmode": "text",
    }
    url = NCBI_EUTILS_BASE + "efetch.fcgi?" + urlencode(params)
    fasta_text = http_get(url)
    if not fasta_text.strip().startswith(">"):
        raise ValueError("NCBI nie zwróciło poprawnego FASTA (sprawdź accession/UID).")
    return fasta_text
