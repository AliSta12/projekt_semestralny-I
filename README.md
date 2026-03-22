# DNA Motif Analyzer

Aplikacja w Pythonie do analizy motywów nukleotydowych w sekwencjach DNA z obsługą kodów IUPAC.

---

## Opis projektu

Program umożliwia:

- wczytywanie sekwencji DNA w formacie FASTA  
- wyszukiwanie motywów nukleotydowych (z obsługą IUPAC)  
- analizę liczby wystąpień motywów  
- wizualizację wyników (heatmapa, wykresy słupkowe)  
- eksport wyników do plików CSV, HTML i PDF  
- pobieranie sekwencji z bazy NCBI (na podstawie accession/UID)  

Aplikacja wykorzystuje podejście macierzowe do reprezentacji wyników oraz integrację z biblioteką Matplotlib do wizualizacji danych.

---

## Funkcjonalności

- obsługa kodów IUPAC (np. N, R, Y)  
- liczenie dopasowań z nakładaniem (overlapping)  
- normalizacja wyników (na 1000 nt)  
- interfejs graficzny (Tkinter)  
- eksport danych i wykresów  

---

## Struktura projektu

- `main.py` – punkt startowy aplikacji  
- `gui_app.py` – interfejs użytkownika (GUI)  
- `analysis_engine.py` – logika analizy motywów  
- `fasta_parser.py` – wczytywanie i walidacja FASTA  
- `iupac.py` – dopasowanie motywów z użyciem kodów IUPAC  
- `export_manager.py` – eksport wyników  
- `export_tab.py` – GUI eksportu  
- `ncbi_client.py` – pobieranie danych z NCBI  

---

## Wymagania

- Python 3.10+

Biblioteki:

```
numpy
matplotlib
requests
```

---

## Instalacja

```bash
pip install -r requirements.txt
```

---

## Uruchomienie

```bash
python main.py
```

---

## Przykładowe użycie

1. Wczytaj plik FASTA lub pobierz dane z NCBI  
2. Dodaj motywy (np. ATG, TATAAA)  
3. Uruchom analizę  
4. Przeglądaj wyniki i wizualizacje  
5. Eksportuj dane do pliku  

---

## Autor

**Alicja Stachura-Matyjewicz**

Projekt wykonany w ramach projektu semestralnego  
na kierunku *Analityk danych biologiczno-medycznych*,  
Uniwersytet Koźmińskiego.
