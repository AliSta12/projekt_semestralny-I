Wersja profesjonalna projektu FASTA

Nowa struktura:
- main.py              -> punkt wejścia
- gui_app.py           -> główna aplikacja Tkinter
- fasta_parser.py      -> wczytywanie i walidacja FASTA
- iupac.py             -> logika dopasowań IUPAC
- ncbi_client.py       -> pobieranie FASTA z NCBI
- analysis_engine.py   -> model wyniku + obliczenia macierzy
- export_manager.py    -> eksport tabel, raportów i wykresów
- export_tab.py        -> zakładka GUI do eksportu

Najważniejsze zmiany:
- rozbicie dużego main.py na moduły o jednej odpowiedzialności
- usunięcie martwego kodu i nieużywanych importów
- uproszczenie clear_results i resetu widoków
- wydzielenie logiki FASTA / IUPAC / NCBI poza GUI
- zachowanie dotychczasowego działania programu
