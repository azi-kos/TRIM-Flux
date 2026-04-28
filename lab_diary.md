# Laboratorijski dnevnik — TRIM-Flux

Diplomska naloga, FKKT/FRI, Univerza v Ljubljani  
Mentor: prof. Miha Moškon  
Avtor: Aljaž Koš

---

## 28.4.2026 — Vzpostavitev repozitorija

### Kaj sem naredil
Naredil fork originalnega TRIM repozitorija (uhlerlab/TRIM) na GitHub pod računom azi-kos/TRIM-Flux. Repozitorij sem klonal lokalno in nastavil dva remote-a: `origin` (moj fork) in `upstream` (original). Ustvaril sem 4 branche: `main`, `replicate-original`, `flux-extension`, `rna-flux`.

### Opažanja in rezultati
- Repozitorij vsebuje originalno kodo: `trim.py` (glavni VAE model), `learn_tcr_embedding.py` (CNN za TCR), `analysis/` (evaluacijske skripte)
- `.mm` fajlov ki jih originalna koda bere ni na GEO — avtorji so objavili `data_processing.py` kot primer kako iz surovih podatkov narediti vhodne `.pkl` fajle
- GEO dataset GSE200996 vsebuje: `.h5` fajle (RNA count matrike), `filtered_contig_*.csv.gz` (TCR sekvence), `*meta.data.txt` (metadata po celici)

---

## 28.4.2026 — Preprocessing podatkov (v teku)

### Kaj sem naredil
Napisal preprocessing notebook (`notebooks/01_preprocessing.ipynb`) ki iz surovih GEO podatkov naredi `.pkl` fajle za TRIM.

### Opažanja in rezultati
_(dopolni po zagonu v Colabu)_

---
